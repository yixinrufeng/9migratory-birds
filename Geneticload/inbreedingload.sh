#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "Usage:"
    echo "  $0 A.hom A.renamed.sorted.vcf.gz outgroup.txt sample.txt output_prefix"
    echo
    echo "Example:"
    echo "  $0 A.hom A.renamed.sorted.vcf.gz outgroup.txt sample.txt A"
    exit 1
fi

HOM="$1"
VCF="$2"
OUTGROUP="$3"
SAMPLES="$4"
PREFIX="$5"

# 重要：
# SnpEff 的 ANN 注释描述的是 ALT 相对于 REF 的功能影响。
# 因此默认只统计 outgroup 推断 REF 为祖先、ALT 为衍生的位点。
#
# 如果你想强行统计 outgroup 推断 ALT 为祖先的位点，可以运行：
# STRICT_REF_ANCESTRAL=0 ./count_roh_derived_snpeff_lof_missense_syn.sh ...
#
# 但通常不建议这样做，因为此时 SnpEff 注释的 ALT 功能影响
# 不等于衍生等位基因的功能影响。
STRICT_REF_ANCESTRAL="${STRICT_REF_ANCESTRAL:-1}"

python3 - "$HOM" "$VCF" "$OUTGROUP" "$SAMPLES" "$PREFIX" "$STRICT_REF_ANCESTRAL" <<'PY'
import sys
import gzip
from collections import defaultdict
from bisect import bisect_right

hom_file, vcf_file, outgroup_file, sample_file, prefix, strict_ref_ancestral = sys.argv[1:7]
strict_ref_ancestral = str(strict_ref_ancestral) == "1"

# ============================================================
# 1. 注释类别定义
# ============================================================

# LoF 定义：
# 这里保留你原脚本中的 LoF 定义。
# 如果你想使用更严格的 LoF，可以删除 inframe_insertion、
# inframe_deletion、splice_region_variant 等较宽松类别。
lof_terms = {
    "transcript_ablation",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_lost",
    "stop_gained",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "splice_region_variant",
    "conservative_inframe_insertion",
    "disruptive_inframe_insertion",
    "conservative_inframe_deletion",
    "disruptive_inframe_deletion"
}

# missense 单独计算，不包含 LoF
missense_terms = {
    "missense_variant"
}

# synonymous 单独计算
syn_terms = {
    "synonymous_variant"
}

# ============================================================
# 2. 基础函数
# ============================================================

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_list(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line)
    return out


def parse_info(info_str):
    info = {}
    for x in info_str.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            info[k] = v
        else:
            info[x] = True
    return info


def get_gt(sample_field, gt_index):
    parts = sample_field.split(":")
    if gt_index >= len(parts):
        return None
    gt = parts[gt_index]
    if gt in (".", "./.", ".|."):
        return None
    return gt.replace("|", "/")


def is_hom_gt(gt):
    if gt is None:
        return False
    a = gt.split("/")
    if len(a) != 2:
        return False
    return a[0] == a[1] and a[0] in ("0", "1")


def infer_ancestral_from_outgroups(fields, outgroup_indices, gt_index):
    """
    规则：
    外群个体中，只要有可用的纯合基因型，就用该纯合 allele 推断祖先状态。

    返回：
      0 = REF 为祖先
      1 = ALT 为祖先
      None = 外群都没有可用纯合基因型
      conflict = 外群中同时出现 0/0 和 1/1，祖先状态冲突
    """
    homo_alleles = []

    for idx in outgroup_indices:
        gt = get_gt(fields[idx], gt_index)
        if not is_hom_gt(gt):
            continue

        allele = gt.split("/")[0]
        homo_alleles.append(int(allele))

    if len(homo_alleles) == 0:
        return None

    if len(set(homo_alleles)) > 1:
        return "conflict"

    return homo_alleles[0]


def is_homozygous_derived(gt, ancestral_code):
    """
    ancestral_code = 0:
      REF 是祖先，ALT 是衍生，纯合衍生 = 1/1

    ancestral_code = 1:
      ALT 是祖先，REF 是衍生，纯合衍生 = 0/0
    """
    if not is_hom_gt(gt):
        return False

    allele = gt.split("/")[0]
    derived_code = str(1 - ancestral_code)

    return allele == derived_code


def get_ann_effects(info, alt):
    """
    从 SnpEff ANN 字段中提取 effect。

    优先使用 ANN 第一列 allele 与当前 ALT 匹配的注释。
    如果没有匹配，则退回使用所有 ANN 注释。
    """
    ann = info.get("ANN", "")
    if ann == "":
        return set(), "no_ann"

    matched = set()
    all_effects = set()

    for item in ann.split(","):
        fs = item.split("|")
        if len(fs) < 2:
            continue

        ann_allele = fs[0]
        effect_field = fs[1]

        effects = set(e for e in effect_field.split("&") if e)
        all_effects.update(effects)

        if ann_allele == alt:
            matched.update(effects)

    if len(matched) > 0:
        return matched, "matched_alt"

    return all_effects, "fallback_all_ann"


def classify_effects(effects):
    """
    将一个位点归为一个互斥类别：

      lof
      missense
      syn
      other

    优先级：
      LoF > missense > synonymous

    也就是说：
      如果一个位点在不同转录本上同时有 LoF 和 missense 注释，
      该位点只归为 LoF，不再计入 missense。

      如果一个位点同时有 missense 和 synonymous 注释，
      该位点只归为 missense，不再计入 synonymous。
    """
    is_lof = len(effects & lof_terms) > 0
    is_missense = len(effects & missense_terms) > 0
    is_syn = len(effects & syn_terms) > 0

    if is_lof:
        return "lof"

    if is_missense:
        return "missense"

    if is_syn:
        return "syn"

    return "other"


def safe_div(a, b):
    if b == 0:
        return "NA"
    return a / b


def fmt(x):
    if x == "NA":
        return "NA"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.10g}"
    return str(x)


# ============================================================
# 3. 读取 ROH 区域
# ============================================================

def load_roh_intervals(hom_file):
    """
    读取 PLINK .hom 文件。
    要求 .hom 文件中包含 IID, CHR, POS1, POS2 四列。

    注意：
    这里不修改染色体名，因此 .hom 与 VCF 中的染色体命名必须一致。
    """
    intervals = defaultdict(lambda: defaultdict(list))

    with open(hom_file) as f:
        header = f.readline().strip().split()

        if len(header) == 0:
            raise RuntimeError("Empty .hom file")

        col = {name: i for i, name in enumerate(header)}

        for needed in ["IID", "CHR", "POS1", "POS2"]:
            if needed not in col:
                raise RuntimeError(f"Cannot find column {needed} in {hom_file}")

        for line in f:
            if not line.strip():
                continue

            fs = line.strip().split()

            sample = fs[col["IID"]]
            chrom = fs[col["CHR"]]
            start = int(float(fs[col["POS1"]]))
            end = int(float(fs[col["POS2"]]))

            if start > end:
                start, end = end, start

            intervals[sample][chrom].append((start, end))

    # 合并每个样本每条染色体上的重叠或相邻 ROH 区间
    merged = defaultdict(dict)

    for sample in intervals:
        for chrom in intervals[sample]:
            arr = sorted(intervals[sample][chrom])
            out = []

            for s, e in arr:
                if not out or s > out[-1][1] + 1:
                    out.append([s, e])
                else:
                    if e > out[-1][1]:
                        out[-1][1] = e

            final = [(s, e) for s, e in out]
            starts = [s for s, e in final]
            merged[sample][chrom] = (final, starts)

    return merged


def in_roh(sample, chrom, pos, roh_intervals):
    if sample not in roh_intervals:
        return False

    if chrom not in roh_intervals[sample]:
        return False

    intervals, starts = roh_intervals[sample][chrom]

    idx = bisect_right(starts, pos) - 1

    if idx < 0:
        return False

    s, e = intervals[idx]

    return s <= pos <= e


# ============================================================
# 4. 主程序
# ============================================================

outgroups = read_list(outgroup_file)
targets = read_list(sample_file)

if len(outgroups) == 0:
    raise RuntimeError("outgroup file is empty")

if len(targets) == 0:
    raise RuntimeError("sample file is empty")

# 如果 sample 文件里包含外群，自动去掉
targets = [x for x in targets if x not in set(outgroups)]

if len(targets) == 0:
    raise RuntimeError("No target samples left after removing outgroups from sample list")

roh_intervals = load_roh_intervals(hom_file)

counts = defaultdict(lambda: defaultdict(int))
stats = defaultdict(int)

seen_header = False

with open_text(vcf_file) as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            seen_header = True

            header = line.split("\t")
            vcf_samples = header[9:]
            sample_to_col = {s: i + 9 for i, s in enumerate(vcf_samples)}

            missing_outgroups = [s for s in outgroups if s not in sample_to_col]
            if missing_outgroups:
                raise RuntimeError(
                    "Outgroup samples not found in VCF: " +
                    ",".join(missing_outgroups)
                )

            missing_targets = [s for s in targets if s not in sample_to_col]
            if missing_targets:
                raise RuntimeError(
                    "Target samples not found in VCF: " +
                    ",".join(missing_targets)
                )

            outgroup_indices = [sample_to_col[s] for s in outgroups]
            target_indices = {s: sample_to_col[s] for s in targets}

            continue

        if not line or line.startswith("#"):
            continue

        if not seen_header:
            raise RuntimeError("VCF header line #CHROM not found before variant records")

        stats["total_vcf_records"] += 1

        fs = line.split("\t")

        if len(fs) < 10:
            stats["skip_malformed_vcf_line"] += 1
            continue

        chrom = fs[0]
        pos = int(fs[1])
        ref = fs[3]
        alt = fs[4]

        # 跳过多等位位点
        if "," in alt:
            stats["skip_multiallelic"] += 1
            continue

        fmt_fields = fs[8].split(":")
        if "GT" not in fmt_fields:
            stats["skip_no_GT"] += 1
            continue

        gt_index = fmt_fields.index("GT")

        ancestral_code = infer_ancestral_from_outgroups(
            fs,
            outgroup_indices,
            gt_index
        )

        if ancestral_code is None:
            stats["skip_no_homozygous_outgroup"] += 1
            continue

        if ancestral_code == "conflict":
            stats["skip_conflicting_outgroups"] += 1
            continue

        if ancestral_code == 0:
            stats["ancestral_REF"] += 1
        elif ancestral_code == 1:
            stats["ancestral_ALT"] += 1

        # 推荐只保留 REF 为祖先的位点
        if strict_ref_ancestral and ancestral_code != 0:
            stats["skip_ancestral_ALT_due_to_STRICT_REF_ANCESTRAL"] += 1
            continue

        info = parse_info(fs[7])
        effects, ann_status = get_ann_effects(info, alt)

        if ann_status == "no_ann":
            stats["skip_no_ANN"] += 1
            continue

        if ann_status == "fallback_all_ann":
            stats["ann_no_alt_match_fallback_all_used"] += 1

        category = classify_effects(effects)

        if category == "other":
            stats["skip_other_annotation"] += 1
            continue

        stats[f"used_{category}_sites"] += 1

        for sample in targets:
            gt = get_gt(fs[target_indices[sample]], gt_index)

            if not is_homozygous_derived(gt, ancestral_code):
                continue

            region = "roh" if in_roh(sample, chrom, pos, roh_intervals) else "nonroh"

            counts[sample][f"{region}_{category}_hom_derived"] += 1


# ============================================================
# 5. 输出结果
# ============================================================

def calc_ratios(c):
    roh_lof_syn = safe_div(
        c["roh_lof_hom_derived"],
        c["roh_syn_hom_derived"]
    )

    nonroh_lof_syn = safe_div(
        c["nonroh_lof_hom_derived"],
        c["nonroh_syn_hom_derived"]
    )

    roh_missense_syn = safe_div(
        c["roh_missense_hom_derived"],
        c["roh_syn_hom_derived"]
    )

    nonroh_missense_syn = safe_div(
        c["nonroh_missense_hom_derived"],
        c["nonroh_syn_hom_derived"]
    )

    if roh_lof_syn == "NA" or nonroh_lof_syn == "NA" or nonroh_lof_syn == 0:
        roh_vs_nonroh_lof_syn = "NA"
    else:
        roh_vs_nonroh_lof_syn = roh_lof_syn / nonroh_lof_syn

    if roh_missense_syn == "NA" or nonroh_missense_syn == "NA" or nonroh_missense_syn == 0:
        roh_vs_nonroh_missense_syn = "NA"
    else:
        roh_vs_nonroh_missense_syn = roh_missense_syn / nonroh_missense_syn

    return {
        "roh_lof_div_syn": roh_lof_syn,
        "nonroh_lof_div_syn": nonroh_lof_syn,
        "roh_vs_nonroh_lof_syn_ratio": roh_vs_nonroh_lof_syn,
        "roh_missense_div_syn": roh_missense_syn,
        "nonroh_missense_div_syn": nonroh_missense_syn,
        "roh_vs_nonroh_missense_syn_ratio": roh_vs_nonroh_missense_syn
    }


columns = [
    "sample",

    "roh_lof_hom_derived",
    "roh_missense_hom_derived",
    "roh_syn_hom_derived",

    "nonroh_lof_hom_derived",
    "nonroh_missense_hom_derived",
    "nonroh_syn_hom_derived",

    "roh_lof_div_syn",
    "nonroh_lof_div_syn",
    "roh_vs_nonroh_lof_syn_ratio",

    "roh_missense_div_syn",
    "nonroh_missense_div_syn",
    "roh_vs_nonroh_missense_syn_ratio"
]

per_sample_out = prefix + ".roh_nonroh.lof_missense_syn.hom_derived_counts.per_sample.tsv"
group_out = prefix + ".roh_nonroh.lof_missense_syn.hom_derived_counts.group_total.tsv"
log_out = prefix + ".roh_nonroh.lof_missense_syn.hom_derived_counts.log"

with open(per_sample_out, "w") as out:
    out.write("\t".join(columns) + "\n")

    for sample in targets:
        c = counts[sample]
        r = calc_ratios(c)

        row = [
            sample,

            c["roh_lof_hom_derived"],
            c["roh_missense_hom_derived"],
            c["roh_syn_hom_derived"],

            c["nonroh_lof_hom_derived"],
            c["nonroh_missense_hom_derived"],
            c["nonroh_syn_hom_derived"],

            r["roh_lof_div_syn"],
            r["nonroh_lof_div_syn"],
            r["roh_vs_nonroh_lof_syn_ratio"],

            r["roh_missense_div_syn"],
            r["nonroh_missense_div_syn"],
            r["roh_vs_nonroh_missense_syn_ratio"]
        ]

        out.write("\t".join(fmt(x) for x in row) + "\n")


group_counts = defaultdict(int)

for sample in targets:
    for k, v in counts[sample].items():
        group_counts[k] += v

group_ratios = calc_ratios(group_counts)

with open(group_out, "w") as out:
    out.write("\t".join(columns) + "\n")

    row = [
        "ALL_TARGETS",

        group_counts["roh_lof_hom_derived"],
        group_counts["roh_missense_hom_derived"],
        group_counts["roh_syn_hom_derived"],

        group_counts["nonroh_lof_hom_derived"],
        group_counts["nonroh_missense_hom_derived"],
        group_counts["nonroh_syn_hom_derived"],

        group_ratios["roh_lof_div_syn"],
        group_ratios["nonroh_lof_div_syn"],
        group_ratios["roh_vs_nonroh_lof_syn_ratio"],

        group_ratios["roh_missense_div_syn"],
        group_ratios["nonroh_missense_div_syn"],
        group_ratios["roh_vs_nonroh_missense_syn_ratio"]
    ]

    out.write("\t".join(fmt(x) for x in row) + "\n")


with open(log_out, "w") as log:
    log.write(f"HOM file: {hom_file}\n")
    log.write(f"VCF file: {vcf_file}\n")
    log.write(f"Outgroup file: {outgroup_file}\n")
    log.write(f"Sample file: {sample_file}\n")
    log.write(f"Number of outgroups: {len(outgroups)}\n")
    log.write(f"Number of target samples: {len(targets)}\n")
    log.write(f"STRICT_REF_ANCESTRAL: {int(strict_ref_ancestral)}\n")
    log.write("\n")

    log.write("Filtering and annotation statistics:\n")
    for k in sorted(stats):
        log.write(f"{k}: {stats[k]}\n")

    log.write("\nCategory definitions:\n")
    log.write("LoF terms:\n")
    log.write(",".join(sorted(lof_terms)) + "\n\n")

    log.write("Missense terms:\n")
    log.write(",".join(sorted(missense_terms)) + "\n\n")

    log.write("Synonymous terms:\n")
    log.write(",".join(sorted(syn_terms)) + "\n\n")

    log.write("Classification priority:\n")
    log.write("LoF > missense > synonymous\n\n")

    log.write("Important note:\n")
    log.write("SnpEff ANN describes ALT relative to REF. STRICT_REF_ANCESTRAL=1 is recommended.\n")
    log.write("With STRICT_REF_ANCESTRAL=1, only sites where REF is inferred as ancestral and ALT as derived are used.\n")

print("Done.")
print("Output files:")
print("  " + per_sample_out)
print("  " + group_out)
print("  " + log_out)
