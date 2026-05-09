#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "Usage:"
    echo "  $0 annotated.biallelic.snp.vcf.gz outgroup.txt target_samples.txt output_prefix"
    echo
    echo "Example:"
    echo "  $0 xx.biallelic.snp.vcf.gz outgroup.txt target_samples.txt result"
    echo
    echo "Default: only sites where outgroup-inferred ancestral allele is REF are used."
    echo "To also keep ancestral-ALT sites, run:"
    echo "  USE_ONLY_REF_ANCESTRAL=0 $0 ..."
    exit 1
fi

VCF="$1"
OUTGROUP="$2"
TARGETS="$3"
PREFIX="$4"

USE_ONLY_REF_ANCESTRAL="${USE_ONLY_REF_ANCESTRAL:-1}"

python3 - "$VCF" "$OUTGROUP" "$TARGETS" "$PREFIX" "$USE_ONLY_REF_ANCESTRAL" <<'PY'
import sys
import gzip
from collections import defaultdict

vcf_file, outgroup_file, target_file, prefix, use_only_ref_ancestral = sys.argv[1:6]
use_only_ref_ancestral = str(use_only_ref_ancestral) == "1"

# LoF definition, following the gorilla paper style and common SnpEff SO terms
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

# Missense means missense_variant only.
# LoF is NOT counted into missense in this modified version.
missense_terms = {
    "missense_variant",
}

syn_terms = {
    "synonymous_variant"
}


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_list(path):
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                vals.append(line)
    return vals


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def get_gt(sample_field, gt_index):
    parts = sample_field.split(":")
    if gt_index >= len(parts):
        return None
    gt = parts[gt_index]
    if gt in (".", "./.", ".|."):
        return None
    return gt.replace("|", "/")


def infer_ancestral_from_outgroups(fields, outgroup_indices, gt_index):
    """
    User rule:
    If either outgroup individual is homozygous, use that homozygous allele as ancestral.
    If both are homozygous but conflict, skip the site.

    Return:
        0 = REF ancestral
        1 = ALT ancestral
        None = no homozygous outgroup
        "conflict" = conflicting homozygous outgroups
    """
    homo = []

    for idx in outgroup_indices:
        gt = get_gt(fields[idx], gt_index)
        if gt is None:
            continue

        alleles = gt.split("/")
        if len(alleles) != 2:
            continue

        a, b = alleles
        if a == b and a in ("0", "1"):
            homo.append(int(a))

    if not homo:
        return None

    unique = set(homo)
    if len(unique) > 1:
        return "conflict"

    return homo[0]


def get_ann_effects(info):
    effects = set()
    ann = info.get("ANN", "")
    if not ann:
        return effects

    for item in ann.split(","):
        fs = item.split("|")
        if len(fs) < 2:
            continue
        for e in fs[1].split("&"):
            if e:
                effects.add(e)

    return effects


def classify_effect(effects):
    """
    Priority for multi-transcript annotations:
      LoF > missense > synonymous > other

    In this modified version:
      LoF and missense are counted separately.
      Missense does NOT include LoF.
    """
    if effects & lof_terms:
        return "lof"
    if effects & missense_terms:
        return "missense"
    if effects & syn_terms:
        return "syn"
    return "other"


def derived_state_for_gt(gt, ancestral_code):
    """
    Return:
      0 = no derived allele
      1 = heterozygous derived
      2 = homozygous derived
      None = missing / non-diploid / unsupported genotype
    """
    alleles = gt.split("/")
    if len(alleles) != 2:
        return None

    if alleles[0] not in ("0", "1") or alleles[1] not in ("0", "1"):
        return None

    derived_code = 1 - ancestral_code
    d = (int(alleles[0]) == derived_code) + (int(alleles[1]) == derived_code)

    return int(d)


def safe_div(a, b):
    if b == 0:
        return "NA"
    return a / b


def calc_metrics(c):
    lof_copies = 2 * c["lof_hom"] + c["lof_het"]
    syn_copies = 2 * c["syn_hom"] + c["syn_het"]
    missense_copies = 2 * c["missense_hom"] + c["missense_het"]

    lof_over_syn = safe_div(lof_copies, syn_copies)
    missense_over_syn = safe_div(missense_copies, syn_copies)

    lof_hom_frac = safe_div(2 * c["lof_hom"], lof_copies)
    syn_hom_frac = safe_div(2 * c["syn_hom"], syn_copies)

    lof_het_frac = safe_div(c["lof_het"], lof_copies)
    syn_het_frac = safe_div(c["syn_het"], syn_copies)

    if lof_hom_frac == "NA" or syn_hom_frac == "NA" or syn_hom_frac == 0:
        hom_frac_ratio = "NA"
    else:
        hom_frac_ratio = lof_hom_frac / syn_hom_frac

    if lof_het_frac == "NA" or syn_het_frac == "NA" or syn_het_frac == 0:
        het_frac_ratio = "NA"
    else:
        het_frac_ratio = lof_het_frac / syn_het_frac

    return {
        "lof_copies": lof_copies,
        "syn_copies": syn_copies,
        "missense_copies": missense_copies,
        "lof_over_syn": lof_over_syn,
        "missense_over_syn": missense_over_syn,
        "lof_hom_frac": lof_hom_frac,
        "syn_hom_frac": syn_hom_frac,
        "lof_hom_frac_over_syn_hom_frac": hom_frac_ratio,
        "lof_het_frac": lof_het_frac,
        "syn_het_frac": syn_het_frac,
        "lof_het_frac_over_syn_het_frac": het_frac_ratio
    }


def fmt(x):
    if x == "NA":
        return "NA"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        # Avoid scientific notation and keep enough precision
        return f"{x:.10f}".rstrip("0").rstrip(".")
    return str(x)


outgroups = read_list(outgroup_file)
targets = read_list(target_file)

if not outgroups:
    sys.exit("ERROR: outgroup.txt is empty.")
if not targets:
    sys.exit("ERROR: target_samples.txt is empty.")

overlap = set(outgroups) & set(targets)
if overlap:
    sys.exit("ERROR: target samples overlap with outgroups: " + ",".join(sorted(overlap)))

# Per-sample counts
counts = {
    s: defaultdict(int) for s in targets
}

stats = defaultdict(int)
header_seen = False

with open_maybe_gzip(vcf_file) as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header_seen = True
            header = line.split("\t")
            samples = header[9:]
            sample_to_col = {s: i + 9 for i, s in enumerate(samples)}

            missing_out = [s for s in outgroups if s not in sample_to_col]
            missing_tar = [s for s in targets if s not in sample_to_col]

            if missing_out:
                sys.exit("ERROR: outgroup samples not found in VCF: " + ",".join(missing_out))
            if missing_tar:
                sys.exit("ERROR: target samples not found in VCF: " + ",".join(missing_tar))

            outgroup_idx = [sample_to_col[s] for s in outgroups]
            target_idx = {s: sample_to_col[s] for s in targets}
            continue

        if not line or line.startswith("#"):
            continue

        if not header_seen:
            sys.exit("ERROR: VCF header line '#CHROM' was not found before variant records.")

        stats["total_records"] += 1
        fields = line.split("\t")

        if len(fields) < 10:
            stats["skip_malformed_record"] += 1
            continue

        ref = fields[3].upper()
        alt = fields[4].upper()

        # Only biallelic SNPs
        if "," in alt or len(ref) != 1 or len(alt) != 1:
            stats["skip_non_biallelic_or_non_snp"] += 1
            continue

        fmt_fields = fields[8].split(":")
        if "GT" not in fmt_fields:
            stats["skip_no_gt"] += 1
            continue
        gt_index = fmt_fields.index("GT")

        ancestral_code = infer_ancestral_from_outgroups(fields, outgroup_idx, gt_index)

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

        # Important:
        # SnpEff ANN describes ALT relative to REF.
        # For derived LoF / missense / synonymous burden, the safest option is to keep
        # only sites where ancestral=REF, so ALT is the derived allele.
        if use_only_ref_ancestral and ancestral_code != 0:
            stats["skip_ancestral_ALT_due_to_strict_mode"] += 1
            continue

        info = parse_info(fields[7])
        effects = get_ann_effects(info)
        if not effects:
            stats["skip_no_ANN"] += 1
            continue

        category = classify_effect(effects)
        if category == "other":
            stats["skip_other_annotation"] += 1
            continue

        stats[f"used_{category}_sites"] += 1

        for s in targets:
            gt = get_gt(fields[target_idx[s]], gt_index)
            if gt is None:
                continue

            dstate = derived_state_for_gt(gt, ancestral_code)
            if dstate is None:
                continue

            if dstate == 0:
                continue

            # dstate = 1: heterozygous derived
            # dstate = 2: homozygous derived
            if category == "lof":
                if dstate == 2:
                    counts[s]["lof_hom"] += 1
                elif dstate == 1:
                    counts[s]["lof_het"] += 1

            elif category == "missense":
                if dstate == 2:
                    counts[s]["missense_hom"] += 1
                elif dstate == 1:
                    counts[s]["missense_het"] += 1

            elif category == "syn":
                if dstate == 2:
                    counts[s]["syn_hom"] += 1
                elif dstate == 1:
                    counts[s]["syn_het"] += 1


# Write per-sample output
out_tsv = prefix + ".derived_burden.per_sample.tsv"
out_group = prefix + ".derived_burden.group_total.tsv"
out_log = prefix + ".derived_burden.log"

columns = [
    "sample",
    "lof_hom", "lof_het", "lof_copies",
    "syn_hom", "syn_het", "syn_copies",
    "missense_hom", "missense_het", "missense_copies",
    "lof_over_syn",
    "missense_over_syn",
    "lof_hom_frac", "syn_hom_frac", "lof_hom_frac_over_syn_hom_frac",
    "lof_het_frac", "syn_het_frac", "lof_het_frac_over_syn_het_frac"
]

with open(out_tsv, "w") as out:
    out.write("\t".join(columns) + "\n")

    for s in targets:
        c = counts[s]
        m = calc_metrics(c)

        row = [
            s,
            c["lof_hom"], c["lof_het"], m["lof_copies"],
            c["syn_hom"], c["syn_het"], m["syn_copies"],
            c["missense_hom"], c["missense_het"], m["missense_copies"],
            m["lof_over_syn"],
            m["missense_over_syn"],
            m["lof_hom_frac"], m["syn_hom_frac"], m["lof_hom_frac_over_syn_hom_frac"],
            m["lof_het_frac"], m["syn_het_frac"], m["lof_het_frac_over_syn_het_frac"]
        ]

        out.write("\t".join(fmt(x) for x in row) + "\n")


# Group total: sum across all target individuals
group_counts = defaultdict(int)
for s in targets:
    for k, v in counts[s].items():
        group_counts[k] += v

group_metrics = calc_metrics(group_counts)

with open(out_group, "w") as out:
    out.write("\t".join(columns) + "\n")

    row = [
        "ALL_TARGETS",
        group_counts["lof_hom"], group_counts["lof_het"], group_metrics["lof_copies"],
        group_counts["syn_hom"], group_counts["syn_het"], group_metrics["syn_copies"],
        group_counts["missense_hom"], group_counts["missense_het"], group_metrics["missense_copies"],
        group_metrics["lof_over_syn"],
        group_metrics["missense_over_syn"],
        group_metrics["lof_hom_frac"], group_metrics["syn_hom_frac"], group_metrics["lof_hom_frac_over_syn_hom_frac"],
        group_metrics["lof_het_frac"], group_metrics["syn_het_frac"], group_metrics["lof_het_frac_over_syn_het_frac"]
    ]

    out.write("\t".join(fmt(x) for x in row) + "\n")


with open(out_log, "w") as log:
    log.write(f"Input VCF: {vcf_file}\n")
    log.write(f"Outgroup samples: {len(outgroups)}\n")
    log.write(f"Target samples: {len(targets)}\n")
    log.write(f"USE_ONLY_REF_ANCESTRAL: {int(use_only_ref_ancestral)}\n")
    log.write("\n")
    for k in sorted(stats):
        log.write(f"{k}: {stats[k]}\n")

print("Done.")
print(f"Per-sample result: {out_tsv}")
print(f"Group-total result: {out_group}")
print(f"Log: {out_log}")
