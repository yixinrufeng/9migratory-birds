#!/usr/bin/env python3
import argparse
import gzip
import random
import math
from collections import defaultdict
from statistics import median

# LoF categories used in the mountain gorilla Science 2015 Fig. 4B-style analysis
LOF_TERMS = {
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

SYN_TERMS = {"synonymous_variant"}


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_sample_list(path):
    samples = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                samples.append(line.split()[0])
    return samples


def get_info_value(info, key):
    prefix = key + "="
    for item in info.split(";"):
        if item.startswith(prefix):
            return item[len(prefix):]
    return None


def classify_snpeff_ann(info, alt):
    """
    Classify variant according to snpEff ANN field.
    Only the ALT allele is interpreted here.
    Therefore, this script only uses sites where outgroup supports REF as ancestral.
    """
    ann = get_info_value(info, "ANN")
    if ann is None:
        return None

    entries = ann.split(",")
    matched_terms = []

    # Prefer ANN records whose allele field matches ALT.
    for e in entries:
        fields = e.split("|")
        if len(fields) < 2:
            continue
        ann_allele = fields[0]
        consequence = fields[1]

        if ann_allele == alt:
            matched_terms.extend(consequence.split("&"))

    # Fallback for indels or cases where snpEff allele string does not exactly match ALT.
    if not matched_terms:
        for e in entries:
            fields = e.split("|")
            if len(fields) >= 2:
                matched_terms.extend(fields[1].split("&"))

    terms = set(matched_terms)

    # If a site is LoF in any transcript, classify as LoF.
    if terms & LOF_TERMS:
        return "LoF"

    if terms & SYN_TERMS:
        return "SYN"

    return None


def parse_gt(sample_field, gt_index):
    fields = sample_field.split(":")
    if gt_index >= len(fields):
        return []

    gt = fields[gt_index]

    if gt in [".", "./.", ".|."]:
        return []

    alleles = gt.replace("|", "/").split("/")

    if "." in alleles:
        return []

    # Only biallelic 0/1 genotypes are accepted.
    if any(a not in ["0", "1"] for a in alleles):
        return []

    return alleles


def af_bin(af, bins):
    b = int(math.floor(af * bins))
    if b < 0:
        b = 0
    if b >= bins:
        b = bins - 1
    return b


def main():
    parser = argparse.ArgumentParser(
        description="Figure 4B-like analysis from one snpEff-annotated VCF with target and outgroup samples."
    )
    parser.add_argument("--vcf", required=True, help="snpEff annotated VCF, .vcf or .vcf.gz")
    parser.add_argument("--target-list", required=True, help="list.txt: target population sample names")
    parser.add_argument("--outgroup-list", required=True, help="outgroup.txt: outgroup sample names")
    parser.add_argument("--prefix", default="fig4b", help="output prefix")
    parser.add_argument("--nrep", type=int, default=10000, help="number of synonymous resampling replicates")
    parser.add_argument("--bins", type=int, default=20, help="AF bins used when exact dAC/AN matching fails")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--outgroup-max-missing", type=float, default=0.2,
                        help="maximum missing allele fraction allowed in outgroup; default 0.2")
    parser.add_argument("--target-max-missing", type=float, default=0.5,
                        help="maximum missing allele fraction allowed in target population; default 0.5")
    parser.add_argument("--include-nonpass", action="store_true",
                        help="include non-PASS variants; default only PASS or '.' variants are used")
    args = parser.parse_args()

    random.seed(args.seed)

    target_samples = read_sample_list(args.target_list)
    outgroup_samples = read_sample_list(args.outgroup_list)

    target_set = set(target_samples)
    outgroup_set = set(outgroup_samples)

    records = {"LoF": [], "SYN": []}

    n_total = 0
    n_biallelic = 0
    n_pass = 0
    n_aa_ref = 0
    n_annotated = 0
    n_target_segregating = 0

    samples = None
    target_indices = []
    outgroup_indices = []

    with open_text(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]

                sample_to_idx = {s: i for i, s in enumerate(samples)}

                missing_target = [s for s in target_samples if s not in sample_to_idx]
                missing_outgroup = [s for s in outgroup_samples if s not in sample_to_idx]

                if missing_target:
                    raise ValueError(
                        "These target samples are not in VCF header: " + ",".join(missing_target)
                    )
                if missing_outgroup:
                    raise ValueError(
                        "These outgroup samples are not in VCF header: " + ",".join(missing_outgroup)
                    )

                target_indices = [sample_to_idx[s] for s in target_samples]
                outgroup_indices = [sample_to_idx[s] for s in outgroup_samples]

                if len(target_indices) == 0:
                    raise ValueError("No target samples found.")
                if len(outgroup_indices) == 0:
                    raise ValueError("No outgroup samples found.")

                continue

            if not line.strip():
                continue

            n_total += 1

            a = line.rstrip("\n").split("\t")
            if len(a) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info, fmt = a[:9]
            genotype_fields = a[9:]

            # Only biallelic variants.
            if "," in alt:
                continue
            n_biallelic += 1

            if not args.include_nonpass:
                if filt not in ["PASS", "."]:
                    continue
            n_pass += 1

            category = classify_snpeff_ann(info, alt)
            if category is None:
                continue
            n_annotated += 1

            fmt_keys = fmt.split(":")
            if "GT" not in fmt_keys:
                continue
            gt_index = fmt_keys.index("GT")

            # ------------------------------------------------------------
            # 1. Outgroup polarization:
            #    keep only sites where called outgroup alleles are all REF.
            #    Thus REF = ancestral, ALT = derived.
            # ------------------------------------------------------------
            og_AN = 0
            og_AC_alt = 0
            og_possible_AN = 2 * len(outgroup_indices)

            for idx in outgroup_indices:
                alleles = parse_gt(genotype_fields[idx], gt_index)
                if not alleles:
                    continue
                og_AN += len(alleles)
                og_AC_alt += alleles.count("1")

            if og_possible_AN == 0:
                continue

            og_missing_frac = 1 - (og_AN / og_possible_AN)

            if og_missing_frac > args.outgroup_max_missing:
                continue

            # Require at least one called outgroup allele.
            if og_AN == 0:
                continue

            # Keep only sites where all called outgroup alleles are REF.
            # Then ALT can be treated as the derived allele.
            if og_AC_alt != 0:
                continue

            n_aa_ref += 1

            # ------------------------------------------------------------
            # 2. Target population:
            #    count derived allele frequency and derived homozygotes.
            # ------------------------------------------------------------
            target_AN = 0
            target_dAC = 0
            has_hom_derived = False
            target_possible_AN = 2 * len(target_indices)

            for idx in target_indices:
                alleles = parse_gt(genotype_fields[idx], gt_index)
                if not alleles:
                    continue

                target_AN += len(alleles)
                target_dAC += alleles.count("1")

                # Diploid derived homozygote: 1/1 or 1|1
                if len(alleles) >= 2 and all(x == "1" for x in alleles):
                    has_hom_derived = True

            if target_possible_AN == 0:
                continue

            target_missing_frac = 1 - (target_AN / target_possible_AN)

            if target_missing_frac > args.target_max_missing:
                continue

            # Use only segregating sites in the target population.
            if not (0 < target_dAC < target_AN):
                continue

            n_target_segregating += 1

            af = target_dAC / target_AN

            records[category].append({
                "chrom": chrom,
                "pos": int(pos),
                "dAC": target_dAC,
                "AN": target_AN,
                "af": af,
                "hom": 1 if has_hom_derived else 0
            })

    lof_sites = records["LoF"]
    syn_sites = records["SYN"]

    # ------------------------------------------------------------
    # Output site-level records.
    # ------------------------------------------------------------
    with open(args.prefix + ".site_records.tsv", "w") as out:
        out.write("category\tchrom\tpos\tdAC\tAN\tAF\thas_derived_homozygote\n")
        for cat in ["LoF", "SYN"]:
            for r in records[cat]:
                out.write(
                    f"{cat}\t{r['chrom']}\t{r['pos']}\t{r['dAC']}\t{r['AN']}\t"
                    f"{r['af']:.8f}\t{r['hom']}\n"
                )

    # ------------------------------------------------------------
    # If no LoF or synonymous sites, stop with summary.
    # ------------------------------------------------------------
    with open(args.prefix + ".filtering_counts.tsv", "w") as out:
        out.write("step\tcount\n")
        out.write(f"all_variant_records\t{n_total}\n")
        out.write(f"biallelic_records\t{n_biallelic}\n")
        out.write(f"pass_or_dot_records\t{n_pass}\n")
        out.write(f"LoF_or_syn_annotated_records\t{n_annotated}\n")
        out.write(f"outgroup_REF_ancestral_records\t{n_aa_ref}\n")
        out.write(f"target_segregating_records\t{n_target_segregating}\n")
        out.write(f"final_LoF_sites\t{len(lof_sites)}\n")
        out.write(f"final_synonymous_sites\t{len(syn_sites)}\n")

    if len(lof_sites) == 0 or len(syn_sites) == 0:
        with open(args.prefix + ".summary.tsv", "w") as out:
            out.write(
                "n_lof_sites\tn_syn_sites\tlof_nhom\tmedian_syn_nhom\t"
                "scaled_lof\tp_value\tnote\n"
            )
            out.write(
                f"{len(lof_sites)}\t{len(syn_sites)}\tNA\tNA\tNA\tNA\t"
                "No LoF sites or no synonymous sites after filtering.\n"
            )
        return

    # Observed LoF n_hom: number of LoF sites with at least one derived homozygote.
    lof_nhom = sum(x["hom"] for x in lof_sites)

    # ------------------------------------------------------------
    # Synonymous resampling matched by derived allele count / allele number.
    # If exact dAC/AN matching is not available, use AF-bin matching.
    # ------------------------------------------------------------
    syn_by_exact = defaultdict(list)
    syn_by_bin = defaultdict(list)

    for s in syn_sites:
        syn_by_exact[(s["dAC"], s["AN"])].append(s)
        syn_by_bin[af_bin(s["af"], args.bins)].append(s)

    available_bins = sorted(syn_by_bin.keys())

    syn_nhom_reps = []
    fallback_count = 0

    for rep in range(args.nrep):
        nhom = 0

        for l in lof_sites:
            pool = syn_by_exact.get((l["dAC"], l["AN"]), [])

            # Fallback to AF-bin matching if no exact dAC/AN synonymous site exists.
            if len(pool) == 0:
                fallback_count += 1
                b = af_bin(l["af"], args.bins)
                if b not in syn_by_bin:
                    b = min(available_bins, key=lambda x: abs(x - b))
                pool = syn_by_bin[b]

            chosen = random.choice(pool)
            nhom += chosen["hom"]

        syn_nhom_reps.append(nhom)

    med_syn = median(syn_nhom_reps)

    if med_syn == 0:
        scaled_lof = float("nan")
    else:
        scaled_lof = lof_nhom / med_syn

    # Same direction as the paper: proportion of synonymous resamples with n_hom < observed LoF n_hom.
    p_value = sum(1 for x in syn_nhom_reps if x < lof_nhom) / args.nrep

    with open(args.prefix + ".syn_sampling.tsv", "w") as out:
        out.write("rep\tsyn_nhom\tscaled_syn_nhom\n")
        for i, nhom in enumerate(syn_nhom_reps, start=1):
            scaled = nhom / med_syn if med_syn != 0 else float("nan")
            out.write(f"{i}\t{nhom}\t{scaled}\n")

    with open(args.prefix + ".summary.tsv", "w") as out:
        out.write(
            "n_lof_sites\tn_syn_sites\tlof_nhom\tmedian_syn_nhom\t"
            "scaled_lof\tp_value\tfallback_count\tnrep\n"
        )
        out.write(
            f"{len(lof_sites)}\t{len(syn_sites)}\t{lof_nhom}\t{med_syn}\t"
            f"{scaled_lof}\t{p_value}\t{fallback_count}\t{args.nrep}\n"
        )


if __name__ == "__main__":
    main()
