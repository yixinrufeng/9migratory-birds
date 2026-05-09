#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "Usage:"
    echo "  $0 annotated.biallelic.snp.vcf.gz popX.txt popY.txt outgroup.txt output_prefix"
    echo
    echo "Example:"
    echo "  $0 xx.biallelic.snp.vcf.gz popX.txt popY.txt outgroup.txt result_XY"
    echo
    echo "Optional:"
    echo "  STRICT_DERIVED_ALT=1 $0 ..."
    echo
    echo "STRICT_DERIVED_ALT=1 means only sites where outgroup-inferred ancestral allele is REF are used."
    echo "This is safer when using SnpEff ANN for LoF because ANN describes ALT relative to REF."
    exit 1
fi

VCF="$1"
POPX="$2"
POPY="$3"
OUTGROUP="$4"
PREFIX="$5"

STRICT_DERIVED_ALT="${STRICT_DERIVED_ALT:-0}"

python3 - "$VCF" "$POPX" "$POPY" "$OUTGROUP" "$PREFIX" "$STRICT_DERIVED_ALT" <<'PY'
import sys
import gzip
import math

vcf_file, popx_file, popy_file, outgroup_file, prefix, strict_derived_alt = sys.argv[1:7]
strict_derived_alt = str(strict_derived_alt) == "1"

# LoF definition follows the terms used in the gorilla supplementary methods,
# with extra SnpEff-specific inframe terms included.
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

missense_terms = {
    "missense_variant"
}

synonymous_terms = {
    "synonymous_variant"
}

intergenic_terms = {
    "intergenic_region"
}


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_sample_list(path):
    samples = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                samples.append(line)
    return samples


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
    User-defined rule:
    If either outgroup individual is homozygous at a SNP, define that homozygous allele as ancestral.

    Return:
        0 if REF is ancestral
        1 if ALT is ancestral
        None if no usable homozygous outgroup genotype
        "conflict" if two outgroups are homozygous for different alleles
    """
    homo_alleles = []

    for idx in outgroup_indices:
        gt = get_gt(fields[idx], gt_index)
        if gt is None:
            continue

        alleles = gt.split("/")
        if len(alleles) == 1:
            # haploid call, e.g. "0" or "1"
            if alleles[0] in ("0", "1"):
                homo_alleles.append(int(alleles[0]))
            continue

        if len(alleles) == 2:
            a, b = alleles
            if a == b and a in ("0", "1"):
                homo_alleles.append(int(a))

    if len(homo_alleles) == 0:
        return None

    unique = set(homo_alleles)
    if len(unique) > 1:
        return "conflict"

    return homo_alleles[0]


def derived_count_for_samples(fields, sample_indices, gt_index, ancestral_code):
    """
    Count derived alleles in a population.
    Biallelic SNP only:
        allele 0 = REF
        allele 1 = ALT

    If ancestral_code == 0, derived allele is ALT.
    If ancestral_code == 1, derived allele is REF.
    """
    derived_code = 1 - ancestral_code

    d = 0
    n = 0

    for idx in sample_indices:
        gt = get_gt(fields[idx], gt_index)
        if gt is None:
            continue

        alleles = gt.split("/")
        for a in alleles:
            if a not in ("0", "1"):
                continue
            n += 1
            if int(a) == derived_code:
                d += 1

    return d, n


def get_ann_effects(info):
    effects = set()
    ann = info.get("ANN", "")
    if not ann:
        return effects

    for ann_item in ann.split(","):
        fields = ann_item.split("|")
        if len(fields) < 2:
            continue
        effect_field = fields[1]
        for e in effect_field.split("&"):
            if e:
                effects.add(e)

    return effects


def classify_site(effects):
    """
    Priority:
        LoF > missense > synonymous > intergenic > other

    This mimics a conservative severity-based choice when multiple SnpEff
    transcript annotations are present.
    """
    if effects & lof_terms:
        return "lof"
    if effects & missense_terms:
        return "missense"
    if effects & synonymous_terms:
        return "synonymous"
    if effects and effects.issubset(intergenic_terms):
        return "intergenic"
    return "other"


popx_samples = read_sample_list(popx_file)
popy_samples = read_sample_list(popy_file)
outgroup_samples = read_sample_list(outgroup_file)

if len(outgroup_samples) < 1:
    sys.exit("ERROR: outgroup.txt is empty.")

overlap = set(popx_samples) & set(popy_samples)
if overlap:
    sys.exit("ERROR: popX and popY sample lists overlap: " + ",".join(sorted(overlap)))

if set(outgroup_samples) & set(popx_samples):
    sys.exit("ERROR: outgroup samples overlap with popX samples.")
if set(outgroup_samples) & set(popy_samples):
    sys.exit("ERROR: outgroup samples overlap with popY samples.")

sums = {
    "missense":   {"xy": 0.0, "yx": 0.0, "n": 0},
    "lof":        {"xy": 0.0, "yx": 0.0, "n": 0},
    "synonymous": {"xy": 0.0, "yx": 0.0, "n": 0},
    "intergenic": {"xy": 0.0, "yx": 0.0, "n": 0}
}

stats = {
    "total_records": 0,
    "used_sites": 0,
    "non_biallelic_or_non_snp": 0,
    "no_gt": 0,
    "no_homozygous_outgroup": 0,
    "conflicting_outgroups": 0,
    "ancestral_ref": 0,
    "ancestral_alt": 0,
    "strict_skip_alt_ancestral": 0,
    "no_ann": 0,
    "other_category": 0,
    "no_called_pop": 0
}

with open_maybe_gzip(vcf_file) as f:
    sample_names = None

    for line in f:
        line = line.rstrip("\n")

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header = line.split("\t")
            sample_names = header[9:]
            sample_to_col = {s: i + 9 for i, s in enumerate(sample_names)}

            for label, sample_list in [
                ("popX", popx_samples),
                ("popY", popy_samples),
                ("outgroup", outgroup_samples)
            ]:
                missing = [s for s in sample_list if s not in sample_to_col]
                if missing:
                    sys.exit(f"ERROR: samples in {label} list not found in VCF: " + ",".join(missing))

            popx_idx = [sample_to_col[s] for s in popx_samples]
            popy_idx = [sample_to_col[s] for s in popy_samples]
            outgroup_idx = [sample_to_col[s] for s in outgroup_samples]
            continue

        if not line or line.startswith("#"):
            continue

        stats["total_records"] += 1
        fields = line.split("\t")

        if len(fields) < 10:
            continue

        ref = fields[3].upper()
        alt = fields[4].upper()

        # Restrict to biallelic SNPs
        if "," in alt or len(ref) != 1 or len(alt) != 1:
            stats["non_biallelic_or_non_snp"] += 1
            continue

        fmt = fields[8].split(":")
        if "GT" not in fmt:
            stats["no_gt"] += 1
            continue
        gt_index = fmt.index("GT")

        ancestral_code = infer_ancestral_from_outgroups(fields, outgroup_idx, gt_index)

        if ancestral_code is None:
            stats["no_homozygous_outgroup"] += 1
            continue
        if ancestral_code == "conflict":
            stats["conflicting_outgroups"] += 1
            continue

        if ancestral_code == 0:
            stats["ancestral_ref"] += 1
        elif ancestral_code == 1:
            stats["ancestral_alt"] += 1

        if strict_derived_alt and ancestral_code != 0:
            stats["strict_skip_alt_ancestral"] += 1
            continue

        info = parse_info(fields[7])
        effects = get_ann_effects(info)

        if not effects:
            stats["no_ann"] += 1
            continue

        category = classify_site(effects)

        if category == "other":
            stats["other_category"] += 1
            continue

        dx, nx = derived_count_for_samples(fields, popx_idx, gt_index, ancestral_code)
        dy, ny = derived_count_for_samples(fields, popy_idx, gt_index, ancestral_code)

        if nx == 0 or ny == 0:
            stats["no_called_pop"] += 1
            continue

        fx = dx / nx
        fy = dy / ny

        xy = fx * (1.0 - fy)
        yx = fy * (1.0 - fx)

        sums[category]["xy"] += xy
        sums[category]["yx"] += yx
        sums[category]["n"] += 1
        stats["used_sites"] += 1


neutral_xy = sums["intergenic"]["xy"]
neutral_yx = sums["intergenic"]["yx"]

out_file = prefix + ".Rxy.tsv"
log_file = prefix + ".Rxy.log"

with open(out_file, "w") as out:
    out.write(
        "category\t"
        "n_sites\t"
        "sum_XY_C\t"
        "sum_YX_C\t"
        "sum_XY_intergenic\t"
        "sum_YX_intergenic\t"
        "L_XY\t"
        "L_YX\t"
        "R_X_over_Y\n"
    )

    for cat in ["synonymous", "missense", "lof"]:
        c_xy = sums[cat]["xy"]
        c_yx = sums[cat]["yx"]

        if neutral_xy == 0 or neutral_yx == 0 or c_xy == 0 or c_yx == 0:
            L_xy = float("nan")
            L_yx = float("nan")
            R = float("nan")
        else:
            L_xy = c_xy / neutral_xy
            L_yx = c_yx / neutral_yx
            R = L_xy / L_yx

        out.write(
            f"{cat}\t"
            f"{sums[cat]['n']}\t"
            f"{c_xy:.10g}\t"
            f"{c_yx:.10g}\t"
            f"{neutral_xy:.10g}\t"
            f"{neutral_yx:.10g}\t"
            f"{L_xy:.10g}\t"
            f"{L_yx:.10g}\t"
            f"{R:.10g}\n"
        )

with open(log_file, "w") as log:
    log.write(f"Input VCF: {vcf_file}\n")
    log.write(f"Population X samples: {len(popx_samples)}\n")
    log.write(f"Population Y samples: {len(popy_samples)}\n")
    log.write(f"Outgroup samples: {len(outgroup_samples)}\n")
    log.write(f"STRICT_DERIVED_ALT: {int(strict_derived_alt)}\n")
    log.write("\n")
    for k, v in stats.items():
        log.write(f"{k}: {v}\n")
    log.write("\n")
    for cat in ["intergenic", "synonymous", "missense", "lof"]:
        log.write(
            f"{cat}: n_sites={sums[cat]['n']}, "
            f"sum_XY={sums[cat]['xy']:.10g}, "
            f"sum_YX={sums[cat]['yx']:.10g}\n"
        )

print("Done.")
print(f"Result: {out_file}")
print(f"Log:    {log_file}")

if sums["intergenic"]["xy"] == 0 or sums["intergenic"]["yx"] == 0:
    print("WARNING: intergenic denominator is zero. R values will be NA.")
    print("Check whether your VCF contains intergenic_region annotations.")
