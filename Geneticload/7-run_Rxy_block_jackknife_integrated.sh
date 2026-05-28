#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Integrated Rxy / R_A/B calculation with block jackknife
#
# Usage:
#   bash run_Rxy_block_jackknife_integrated.sh \
#     input.snpeff.ann.vcf.gz \
#     popX.txt \
#     popY.txt \
#     outgroup.txt \
#     output_prefix
#
# Example:
#   STRICT_DERIVED_ALT=1 \
#   REQUIRE_NO_MISSING=1 \
#   N_BLOCKS=100 \
#   BLOCK_MODE=equal_sites \
#   ./run_Rxy_block_jackknife_integrated.sh \
#     all.snpeff.ann.vcf.gz \
#     Spoon-billed.txt \
#     Red-necked.txt \
#     outgroup.txt \
#     Spoon_vs_Red
#
# Optional environment variables:
#   STRICT_DERIVED_ALT=1          default: 1
#   REQUIRE_NO_MISSING=0          default: 0
#   N_BLOCKS=100                  default: 100
#   BLOCK_MODE=equal_sites        equal_sites or bp
#   BLOCK_SIZE_BP=5000000         only used when BLOCK_MODE=bp
#   KEEP_TMP=0                    default: 0
#
# Output:
#   output_prefix.all_sites.Rxy.tsv
#   output_prefix.block_jackknife.Rxy.tsv
#   output_prefix.block_jackknife.summary.Rxy.tsv
#   output_prefix.Rxy.log
# ============================================================

if [ "$#" -lt 5 ]; then
    echo "Usage:"
    echo "  $0 input.snpeff.ann.vcf.gz popX.txt popY.txt outgroup.txt output_prefix"
    echo
    echo "Example:"
    echo "  STRICT_DERIVED_ALT=1 REQUIRE_NO_MISSING=1 N_BLOCKS=100 BLOCK_MODE=equal_sites \\"
    echo "  $0 all.snpeff.ann.vcf.gz Grey-tailed.txt Common_Snipe.txt outgroup.txt Grey_vs_Common"
    echo
    echo "Optional environment variables:"
    echo "  STRICT_DERIVED_ALT=1          Only use sites where ancestral allele is REF; recommended for SnpEff ANN."
    echo "  REQUIRE_NO_MISSING=0          If 1, remove sites with any missing genotype among selected samples."
    echo "  N_BLOCKS=100                  Number of jackknife blocks when BLOCK_MODE=equal_sites."
    echo "  BLOCK_MODE=equal_sites        equal_sites or bp."
    echo "  BLOCK_SIZE_BP=5000000         Physical block size when BLOCK_MODE=bp."
    echo "  KEEP_TMP=0                    If 1, keep temporary files."
    exit 1
fi

VCF="$1"
POPX="$2"
POPY="$3"
OUTGROUP="$4"
PREFIX="$5"

STRICT_DERIVED_ALT="${STRICT_DERIVED_ALT:-1}"
REQUIRE_NO_MISSING="${REQUIRE_NO_MISSING:-0}"
N_BLOCKS="${N_BLOCKS:-100}"
BLOCK_MODE="${BLOCK_MODE:-equal_sites}"
BLOCK_SIZE_BP="${BLOCK_SIZE_BP:-5000000}"
KEEP_TMP="${KEEP_TMP:-0}"

command -v bcftools >/dev/null 2>&1 || {
    echo "ERROR: bcftools not found in PATH." >&2
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 not found in PATH." >&2
    exit 1
}

WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/$(basename "$PREFIX").RxyJK.XXXXXX")

cleanup() {
    if [ "$KEEP_TMP" = "1" ]; then
        echo "Temporary directory kept: $WORKDIR"
    else
        rm -rf "$WORKDIR"
    fi
}
trap cleanup EXIT

SAMPLE_LIST="$WORKDIR/selected_samples.txt"
PRE_VCF="$WORKDIR/input.selected.biallelic.snp.vcf.gz"

cat "$POPX" "$POPY" "$OUTGROUP" \
    | awk 'NF && $1 !~ /^#/ {print $1}' \
    | sort -u > "$SAMPLE_LIST"

echo "[1/3] Preprocessing VCF with bcftools..."

if [ "$REQUIRE_NO_MISSING" = "1" ]; then
    bcftools view -S "$SAMPLE_LIST" -Ou "$VCF" \
        | bcftools norm -m -any -Ou \
        | bcftools view -m2 -M2 -v snps -Ou \
        | bcftools view -e 'F_MISSING>0' -Oz -o "$PRE_VCF"
else
    bcftools view -S "$SAMPLE_LIST" -Ou "$VCF" \
        | bcftools norm -m -any -Ou \
        | bcftools view -m2 -M2 -v snps -Oz -o "$PRE_VCF"
fi

echo "[2/3] Calculating site-level contributions and block jackknife..."

python3 - "$PRE_VCF" "$POPX" "$POPY" "$OUTGROUP" "$PREFIX" \
    "$STRICT_DERIVED_ALT" "$REQUIRE_NO_MISSING" "$N_BLOCKS" \
    "$BLOCK_MODE" "$BLOCK_SIZE_BP" <<'PY'

import sys
import gzip
import math
from array import array

(
    vcf_file,
    popx_file,
    popy_file,
    outgroup_file,
    prefix,
    strict_derived_alt,
    require_no_missing,
    n_blocks_requested,
    block_mode,
    block_size_bp
) = sys.argv[1:11]

strict_derived_alt = str(strict_derived_alt) == "1"
require_no_missing = str(require_no_missing) == "1"
n_blocks_requested = int(n_blocks_requested)
block_mode = str(block_mode)
block_size_bp = int(block_size_bp)

if block_mode not in ("equal_sites", "bp"):
    sys.exit("ERROR: BLOCK_MODE must be either 'equal_sites' or 'bp'.")

# ------------------------------------------------------------
# Annotation categories
# ------------------------------------------------------------

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

cat_to_code = {
    "intergenic": 0,
    "synonymous": 1,
    "missense": 2,
    "lof": 3
}

code_to_cat = {
    0: "intergenic",
    1: "synonymous",
    2: "missense",
    3: "lof"
}

output_categories = ["synonymous", "missense", "lof"]


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

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
    If any outgroup individual is homozygous at a SNP, define that homozygous
    allele as ancestral.

    Return:
        0          REF is ancestral
        1          ALT is ancestral
        None       no usable homozygous outgroup genotype
        conflict   homozygous outgroups support different alleles
    """
    homo_alleles = []

    for idx in outgroup_indices:
        gt = get_gt(fields[idx], gt_index)
        if gt is None:
            continue

        alleles = gt.split("/")

        if len(alleles) == 1:
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
    derived_code = 1 - ancestral_code

    d = 0
    n = 0

    for idx in sample_indices:
        gt = get_gt(fields[idx], gt_index)
        if gt is None:
            continue

        for a in gt.split("/"):
            if a not in ("0", "1"):
                continue

            n += 1

            if int(a) == derived_code:
                d += 1

    return d, n


def get_ann_effects_for_current_alt(info, alt):
    """
    SnpEff ANN format:
        Allele | Annotation | Annotation_Impact | Gene_Name | ...
    """
    effects = set()

    ann = info.get("ANN", "")
    if not ann:
        return effects, "no_ann"

    saw_ann = False
    saw_matching_alt = False

    alt = alt.upper()

    for ann_item in ann.split(","):
        fields = ann_item.split("|")
        if len(fields) < 2:
            continue

        saw_ann = True

        ann_allele = fields[0].upper()
        effect_field = fields[1]

        if ann_allele and ann_allele != alt:
            continue

        saw_matching_alt = True

        for e in effect_field.split("&"):
            if e:
                effects.add(e)

    if not saw_ann:
        return effects, "no_ann"

    if not saw_matching_alt:
        return effects, "no_alt_matching_ann"

    return effects, "ok"


def classify_site(effects):
    """
    Priority:
        LoF > missense > synonymous > intergenic > other
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


def fmt_num(x):
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return "NA"
    return f"{x:.10g}"


def is_finite_number(x):
    return isinstance(x, (float, int)) and math.isfinite(x)


# ------------------------------------------------------------
# Read sample lists
# ------------------------------------------------------------

popx_samples = read_sample_list(popx_file)
popy_samples = read_sample_list(popy_file)
outgroup_samples = read_sample_list(outgroup_file)

if len(popx_samples) == 0:
    sys.exit("ERROR: popX sample list is empty.")

if len(popy_samples) == 0:
    sys.exit("ERROR: popY sample list is empty.")

if len(outgroup_samples) == 0:
    sys.exit("ERROR: outgroup sample list is empty.")

overlap_xy = set(popx_samples) & set(popy_samples)
if overlap_xy:
    sys.exit("ERROR: popX and popY sample lists overlap: " + ",".join(sorted(overlap_xy)))

if set(outgroup_samples) & set(popx_samples):
    sys.exit("ERROR: outgroup samples overlap with popX samples.")

if set(outgroup_samples) & set(popy_samples):
    sys.exit("ERROR: outgroup samples overlap with popY samples.")


# ------------------------------------------------------------
# Store usable site-level contributions
# ------------------------------------------------------------

site_cat = bytearray()
site_xy = array("d")
site_yx = array("d")
site_chrom = []
site_pos = array("I")

stats = {
    "total_records_after_bcftools": 0,
    "non_biallelic_or_non_snp": 0,
    "no_gt": 0,
    "no_homozygous_outgroup": 0,
    "conflicting_outgroups": 0,
    "ancestral_ref": 0,
    "ancestral_alt": 0,
    "strict_skip_alt_ancestral": 0,
    "no_ann": 0,
    "no_alt_matching_ann": 0,
    "other_category": 0,
    "no_called_pop": 0,
    "used_sites": 0
}

category_counts = {
    "intergenic": 0,
    "synonymous": 0,
    "missense": 0,
    "lof": 0
}

category_sum_xy = {
    "intergenic": 0.0,
    "synonymous": 0.0,
    "missense": 0.0,
    "lof": 0.0
}

category_sum_yx = {
    "intergenic": 0.0,
    "synonymous": 0.0,
    "missense": 0.0,
    "lof": 0.0
}


# ------------------------------------------------------------
# Read VCF
# ------------------------------------------------------------

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
                    sys.exit(
                        f"ERROR: samples in {label} list not found in VCF: "
                        + ",".join(missing)
                    )

            popx_idx = [sample_to_col[s] for s in popx_samples]
            popy_idx = [sample_to_col[s] for s in popy_samples]
            outgroup_idx = [sample_to_col[s] for s in outgroup_samples]
            continue

        if not line or line.startswith("#"):
            continue

        stats["total_records_after_bcftools"] += 1

        fields = line.split("\t")

        if len(fields) < 10:
            continue

        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3].upper()
        alt = fields[4].upper()

        if "," in alt or len(ref) != 1 or len(alt) != 1:
            stats["non_biallelic_or_non_snp"] += 1
            continue

        fmt = fields[8].split(":")
        if "GT" not in fmt:
            stats["no_gt"] += 1
            continue

        gt_index = fmt.index("GT")

        ancestral_code = infer_ancestral_from_outgroups(
            fields,
            outgroup_idx,
            gt_index
        )

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
        effects, ann_status = get_ann_effects_for_current_alt(info, alt)

        if ann_status == "no_ann":
            stats["no_ann"] += 1
            continue

        if ann_status == "no_alt_matching_ann":
            stats["no_alt_matching_ann"] += 1
            continue

        category = classify_site(effects)

        if category == "other":
            stats["other_category"] += 1
            continue

        dx, nx = derived_count_for_samples(
            fields,
            popx_idx,
            gt_index,
            ancestral_code
        )

        dy, ny = derived_count_for_samples(
            fields,
            popy_idx,
            gt_index,
            ancestral_code
        )

        if nx == 0 or ny == 0:
            stats["no_called_pop"] += 1
            continue

        fx = dx / nx
        fy = dy / ny

        xy = fx * (1.0 - fy)
        yx = fy * (1.0 - fx)

        code = cat_to_code[category]

        site_cat.append(code)
        site_xy.append(xy)
        site_yx.append(yx)
        site_chrom.append(chrom)
        site_pos.append(pos)

        stats["used_sites"] += 1
        category_counts[category] += 1
        category_sum_xy[category] += xy
        category_sum_yx[category] += yx


n_available = len(site_cat)

if n_available == 0:
    sys.exit("ERROR: no usable sites after filtering.")

if n_available < 2:
    sys.exit("ERROR: fewer than 2 usable sites; block jackknife cannot be performed.")


# ------------------------------------------------------------
# Assign blocks
# ------------------------------------------------------------

site_block = array("I")

if block_mode == "equal_sites":
    actual_n_blocks = min(n_blocks_requested, n_available)

    if actual_n_blocks < 2:
        sys.exit("ERROR: fewer than 2 jackknife blocks.")

    for i in range(n_available):
        b = (i * actual_n_blocks) // n_available
        if b >= actual_n_blocks:
            b = actual_n_blocks - 1
        site_block.append(b)

else:
    block_to_id = {}

    for chrom, pos in zip(site_chrom, site_pos):
        key = (chrom, (pos - 1) // block_size_bp)
        if key not in block_to_id:
            block_to_id[key] = len(block_to_id)
        site_block.append(block_to_id[key])

    actual_n_blocks = len(block_to_id)

    if actual_n_blocks < 2:
        sys.exit("ERROR: fewer than 2 physical blocks. Try smaller BLOCK_SIZE_BP.")


# ------------------------------------------------------------
# Rxy calculation from summed contributions
# ------------------------------------------------------------

def calculate_from_sums(n_by_code, xy_by_code, yx_by_code):
    neutral_xy = xy_by_code[cat_to_code["intergenic"]]
    neutral_yx = yx_by_code[cat_to_code["intergenic"]]
    n_intergenic = n_by_code[cat_to_code["intergenic"]]

    rows = []

    for cat in output_categories:
        c = cat_to_code[cat]

        c_xy = xy_by_code[c]
        c_yx = yx_by_code[c]

        if neutral_xy == 0 or neutral_yx == 0 or c_xy == 0 or c_yx == 0:
            L_xy = float("nan")
            L_yx = float("nan")
            R = float("nan")
        else:
            L_xy = c_xy / neutral_xy
            L_yx = c_yx / neutral_yx
            R = L_xy / L_yx

        rows.append({
            "category": cat,
            "n_sites": n_by_code[c],
            "n_intergenic": n_intergenic,
            "sum_XY_C": c_xy,
            "sum_YX_C": c_yx,
            "sum_XY_intergenic": neutral_xy,
            "sum_YX_intergenic": neutral_yx,
            "L_XY": L_xy,
            "L_YX": L_yx,
            "R_X_over_Y": R
        })

    return rows


def write_rows(out, replicate, rows):
    for row in rows:
        out.write(
            f"{replicate}\t"
            f"{row['category']}\t"
            f"{row['n_sites']}\t"
            f"{row['n_intergenic']}\t"
            f"{fmt_num(row['sum_XY_C'])}\t"
            f"{fmt_num(row['sum_YX_C'])}\t"
            f"{fmt_num(row['sum_XY_intergenic'])}\t"
            f"{fmt_num(row['sum_YX_intergenic'])}\t"
            f"{fmt_num(row['L_XY'])}\t"
            f"{fmt_num(row['L_YX'])}\t"
            f"{fmt_num(row['R_X_over_Y'])}\n"
        )


def write_jk_rows(out, replicate, left_out_block, n_excluded_total, rows):
    for row in rows:
        out.write(
            f"{replicate}\t"
            f"{left_out_block}\t"
            f"{n_excluded_total}\t"
            f"{row['category']}\t"
            f"{row['n_sites']}\t"
            f"{row['n_intergenic']}\t"
            f"{fmt_num(row['sum_XY_C'])}\t"
            f"{fmt_num(row['sum_YX_C'])}\t"
            f"{fmt_num(row['sum_XY_intergenic'])}\t"
            f"{fmt_num(row['sum_YX_intergenic'])}\t"
            f"{fmt_num(row['L_XY'])}\t"
            f"{fmt_num(row['L_YX'])}\t"
            f"{fmt_num(row['R_X_over_Y'])}\n"
        )


header = [
    "replicate",
    "category",
    "n_sites",
    "n_intergenic",
    "sum_XY_C",
    "sum_YX_C",
    "sum_XY_intergenic",
    "sum_YX_intergenic",
    "L_XY",
    "L_YX",
    "R_X_over_Y"
]

jk_header = [
    "replicate",
    "left_out_block",
    "n_excluded_total",
    "category",
    "n_sites",
    "n_intergenic",
    "sum_XY_C",
    "sum_YX_C",
    "sum_XY_intergenic",
    "sum_YX_intergenic",
    "L_XY",
    "L_YX",
    "R_X_over_Y"
]


# ------------------------------------------------------------
# Total sums
# ------------------------------------------------------------

total_n_by_code = [0, 0, 0, 0]
total_xy_by_code = [0.0, 0.0, 0.0, 0.0]
total_yx_by_code = [0.0, 0.0, 0.0, 0.0]

for code in range(4):
    cat = code_to_cat[code]
    total_n_by_code[code] = category_counts[cat]
    total_xy_by_code[code] = category_sum_xy[cat]
    total_yx_by_code[code] = category_sum_yx[cat]


# ------------------------------------------------------------
# Per-block sums
# ------------------------------------------------------------

block_n = [[0, 0, 0, 0] for _ in range(actual_n_blocks)]
block_xy = [[0.0, 0.0, 0.0, 0.0] for _ in range(actual_n_blocks)]
block_yx = [[0.0, 0.0, 0.0, 0.0] for _ in range(actual_n_blocks)]

for i in range(n_available):
    b = site_block[i]
    c = site_cat[i]

    block_n[b][c] += 1
    block_xy[b][c] += site_xy[i]
    block_yx[b][c] += site_yx[i]


# ------------------------------------------------------------
# Output all-site result
# ------------------------------------------------------------

all_out_file = prefix + ".all_sites.Rxy.tsv"
jk_out_file = prefix + ".block_jackknife.Rxy.tsv"
jk_summary_file = prefix + ".block_jackknife.summary.Rxy.tsv"
log_file = prefix + ".Rxy.log"

all_rows = calculate_from_sums(
    total_n_by_code,
    total_xy_by_code,
    total_yx_by_code
)

R_all_by_cat = {
    row["category"]: row["R_X_over_Y"]
    for row in all_rows
}

with open(all_out_file, "w") as out:
    out.write("\t".join(header) + "\n")
    write_rows(out, "all", all_rows)


# ------------------------------------------------------------
# Output leave-one-block-out jackknife replicates
# ------------------------------------------------------------

jk_values_by_cat = {cat: [] for cat in output_categories}

with open(jk_out_file, "w") as out:
    out.write("\t".join(jk_header) + "\n")

    for b in range(actual_n_blocks):
        leave_n = [
            total_n_by_code[c] - block_n[b][c]
            for c in range(4)
        ]

        leave_xy = [
            total_xy_by_code[c] - block_xy[b][c]
            for c in range(4)
        ]

        leave_yx = [
            total_yx_by_code[c] - block_yx[b][c]
            for c in range(4)
        ]

        n_excluded_total = sum(block_n[b])

        rows = calculate_from_sums(
            leave_n,
            leave_xy,
            leave_yx
        )

        for row in rows:
            if is_finite_number(row["R_X_over_Y"]):
                jk_values_by_cat[row["category"]].append(row["R_X_over_Y"])

        write_jk_rows(
            out,
            b + 1,
            b + 1,
            n_excluded_total,
            rows
        )


# ------------------------------------------------------------
# Jackknife summary
# ------------------------------------------------------------

def mean(values):
    return sum(values) / len(values) if values else float("nan")


with open(jk_summary_file, "w") as out:
    out.write(
        "category\t"
        "R_all\t"
        "n_blocks_total\t"
        "n_blocks_used\t"
        "jackknife_mean\t"
        "jackknife_se\t"
        "lower_1SE\t"
        "upper_1SE\t"
        "lower_2SE\t"
        "upper_2SE\t"
        "lower_95_normal\t"
        "upper_95_normal\n"
    )

    for cat in output_categories:
        values = jk_values_by_cat[cat]
        values = [v for v in values if is_finite_number(v)]

        theta_all = R_all_by_cat.get(cat, float("nan"))

        if len(values) < 2 or not is_finite_number(theta_all):
            jk_mean = float("nan")
            jk_se = float("nan")
        else:
            b_eff = len(values)
            jk_mean = mean(values)
            jk_var = (b_eff - 1) / b_eff * sum(
                (v - jk_mean) ** 2 for v in values
            )
            jk_se = math.sqrt(jk_var)

        out.write(
            f"{cat}\t"
            f"{fmt_num(theta_all)}\t"
            f"{actual_n_blocks}\t"
            f"{len(values)}\t"
            f"{fmt_num(jk_mean)}\t"
            f"{fmt_num(jk_se)}\t"
            f"{fmt_num(theta_all - jk_se if is_finite_number(jk_se) else float('nan'))}\t"
            f"{fmt_num(theta_all + jk_se if is_finite_number(jk_se) else float('nan'))}\t"
            f"{fmt_num(theta_all - 2 * jk_se if is_finite_number(jk_se) else float('nan'))}\t"
            f"{fmt_num(theta_all + 2 * jk_se if is_finite_number(jk_se) else float('nan'))}\t"
            f"{fmt_num(theta_all - 1.96 * jk_se if is_finite_number(jk_se) else float('nan'))}\t"
            f"{fmt_num(theta_all + 1.96 * jk_se if is_finite_number(jk_se) else float('nan'))}\n"
        )


# ------------------------------------------------------------
# Log file
# ------------------------------------------------------------

with open(log_file, "w") as log:
    log.write(f"Input VCF after bcftools preprocessing: {vcf_file}\n")
    log.write(f"Population X file: {popx_file}\n")
    log.write(f"Population Y file: {popy_file}\n")
    log.write(f"Outgroup file: {outgroup_file}\n")
    log.write(f"Population X samples: {len(popx_samples)}\n")
    log.write(f"Population Y samples: {len(popy_samples)}\n")
    log.write(f"Outgroup samples: {len(outgroup_samples)}\n")
    log.write("\n")

    log.write("Options\n")
    log.write(f"STRICT_DERIVED_ALT: {int(strict_derived_alt)}\n")
    log.write(f"REQUIRE_NO_MISSING: {int(require_no_missing)}\n")
    log.write(f"N_BLOCKS requested: {n_blocks_requested}\n")
    log.write(f"N_BLOCKS used: {actual_n_blocks}\n")
    log.write(f"BLOCK_MODE: {block_mode}\n")
    log.write(f"BLOCK_SIZE_BP: {block_size_bp}\n")
    log.write("\n")

    log.write("Filtering statistics\n")
    for k, v in stats.items():
        log.write(f"{k}: {v}\n")
    log.write("\n")

    log.write("Usable site counts and sums\n")
    for cat in ["intergenic", "synonymous", "missense", "lof"]:
        log.write(
            f"{cat}: "
            f"n_sites={category_counts[cat]}, "
            f"sum_XY={category_sum_xy[cat]:.10g}, "
            f"sum_YX={category_sum_yx[cat]:.10g}\n"
        )

    log.write("\n")

    log.write("Block summary\n")
    for b in range(actual_n_blocks):
        log.write(
            f"block_{b + 1}: "
            f"n_total={sum(block_n[b])}, "
            f"n_intergenic={block_n[b][cat_to_code['intergenic']]}, "
            f"n_synonymous={block_n[b][cat_to_code['synonymous']]}, "
            f"n_missense={block_n[b][cat_to_code['missense']]}, "
            f"n_lof={block_n[b][cat_to_code['lof']]}\n"
        )

    log.write("\n")

    if strict_derived_alt:
        log.write(
            "Note: STRICT_DERIVED_ALT=1 was used. "
            "Only sites where outgroup-inferred ancestral allele is REF were retained. "
            "This is recommended because SnpEff ANN describes ALT relative to REF.\n"
        )
    else:
        log.write(
            "Warning: STRICT_DERIVED_ALT=0 was used. "
            "When ancestral allele is ALT, the derived allele is REF, but SnpEff ANN still describes ALT. "
            "This can make LoF/missense/synonymous categories difficult to interpret as derived effects.\n"
        )

    if block_mode == "equal_sites":
        log.write(
            "\nNote: BLOCK_MODE=equal_sites was used. "
            "Usable sites were divided into consecutive blocks according to VCF order. "
            "For genomic block jackknife, the VCF should be sorted by chromosome and position.\n"
        )
    else:
        log.write(
            "\nNote: BLOCK_MODE=bp was used. "
            "Blocks were defined by chromosome and physical windows of BLOCK_SIZE_BP.\n"
        )


print("Done.")
print(f"All-site result:             {all_out_file}")
print(f"Block jackknife result:      {jk_out_file}")
print(f"Block jackknife summary:     {jk_summary_file}")
print(f"Log:                         {log_file}")
print(f"Usable sites:                {n_available}")
print(f"Jackknife blocks:            {actual_n_blocks}")

if category_counts["intergenic"] == 0:
    print("WARNING: no usable intergenic sites. R values will be NA.")

if category_sum_xy["intergenic"] == 0 or category_sum_yx["intergenic"] == 0:
    print("WARNING: intergenic denominator is zero. R values may be NA.")

PY

echo "[3/3] Finished."
