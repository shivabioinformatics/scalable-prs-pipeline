#!/usr/bin/env python3
"""
calculate_prs.py
----------------
This is the heart of the whole pipeline. Everything before this (QC, normalization,
imputation) was just getting the data ready. Now I actually compute the score.

The math is straightforward:
  PRS_j = sum over all SNPs i of (dosage_ij * beta_i)

Where:
  - dosage is how many copies of the effect allele a person carries (0, 1, or 2)
  - beta is the effect size from the GWAS summary statistics

A higher PRS means higher genetic risk. I can then bin people into percentiles
and see who falls in the high-risk tail of the distribution.

I also handle:
  - SNP matching between the VCF and the GWAS stats (they won't always overlap perfectly)
  - Allele flipping — sometimes my VCF has the ref/alt swapped relative to the GWAS
  - P-value filtering — I can threshold which variants to include based on significance
  - LD clumping awareness — noted in comments, would use PLINK in production
"""

import sys
import os
import argparse
import math
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Polygenic Risk Scores")
    parser.add_argument("--vcf", required=True, help="Input VCF (post-QC, post-imputation)")
    parser.add_argument("--gwas", required=True, help="GWAS summary statistics TSV")
    parser.add_argument("--output", required=True, help="Output PRS scores TSV")
    parser.add_argument("--p-threshold", type=float, default=0.05,
                        help="P-value threshold for including variants (default: 0.05)")
    parser.add_argument("--variant-log", default=None,
                        help="Optional: log which variants were used/skipped")
    return parser.parse_args()


def load_gwas_stats(filepath, p_threshold):
    """
    Load the GWAS summary statistics and filter by p-value.

    In a real PRS, you'd also do LD clumping here (using PLINK --clump)
    to avoid counting correlated SNPs twice. For this demo I'm just
    using a p-value cutoff, which is the simplest approach. More
    sophisticated methods like LDpred2 or PRS-CS adjust betas
    directly for LD structure.
    """
    gwas = {}
    total = 0
    kept = 0
    skipped_p = 0

    with open(filepath, "r") as f:
        header = f.readline().strip().split("\t")

        # figure out which columns are which — different GWAS files use different names
        col_map = {}
        for i, col in enumerate(header):
            col_map[col.upper()] = i

        # I need at minimum: SNP ID, effect allele, beta, and p-value
        snp_col = col_map.get("SNP", col_map.get("RSID", col_map.get("VARIANT_ID", 0)))
        ea_col = col_map.get("EFFECT_ALLELE", col_map.get("A1", col_map.get("ALT", 3)))
        oa_col = col_map.get("OTHER_ALLELE", col_map.get("A2", col_map.get("REF", 4)))
        beta_col = col_map.get("BETA", col_map.get("EFFECT", col_map.get("LOG_OR", 5)))
        p_col = col_map.get("P", col_map.get("PVALUE", col_map.get("P_VALUE", 7)))

        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < max(snp_col, ea_col, beta_col, p_col) + 1:
                continue

            total += 1
            rsid = fields[snp_col]
            effect_allele = fields[ea_col]
            other_allele = fields[oa_col]

            try:
                beta = float(fields[beta_col])
                pval = float(fields[p_col])
            except ValueError:
                continue

            # filter on p-value threshold
            # in practice with a real GWAS you'd see most variants fail this filter
            # and only the genome-wide significant (or suggestive) ones make it through
            if pval > p_threshold:
                skipped_p += 1
                continue

            kept += 1
            gwas[rsid] = {
                "effect_allele": effect_allele,
                "other_allele": other_allele,
                "beta": beta,
                "pval": pval
            }

    print(f"  GWAS stats: {total} total, {kept} pass p < {p_threshold}, {skipped_p} filtered out")
    return gwas


def parse_vcf_genotypes(vcf_path):
    """
    Read the VCF and extract genotype dosages for each sample at each variant.

    Dosage = number of alt alleles (0, 1, or 2).
    For imputed data you'd want the dosage field (DS) instead of hard calls,
    but for this demo I'm using hard genotype calls and converting them.
    """
    sample_names = []
    genotype_data = {}  # rsid -> {sample: dosage}
    variant_alleles = {}  # rsid -> (ref, alt)

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                sample_names = fields[9:]
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            rsid = fields[2]
            ref = fields[3]
            alt = fields[4]
            fmt = fields[8].split(":")
            genotypes = fields[9:]

            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            variant_alleles[rsid] = (ref, alt)
            dosages = {}

            for i, geno_str in enumerate(genotypes):
                parts = geno_str.split(":")
                gt = parts[gt_idx] if gt_idx < len(parts) else "./."

                if gt in ("./.", ".|."):
                    # missing genotype — I'll impute with the mean dosage later
                    dosages[sample_names[i]] = None
                else:
                    alleles = gt.replace("|", "/").split("/")
                    # count alt alleles
                    dosage = sum(1 for a in alleles if a.isdigit() and int(a) > 0)
                    dosages[sample_names[i]] = dosage

            genotype_data[rsid] = dosages

    return sample_names, genotype_data, variant_alleles


def calculate_scores(sample_names, genotype_data, variant_alleles, gwas_stats):
    """
    The actual PRS calculation. For each person, I walk through every variant
    that passed QC and is in the GWAS stats, multiply their dosage by the
    effect size (beta), and sum it all up.

    I also handle allele flipping: if the effect allele in the GWAS is my
    REF allele (not ALT), I need to flip the dosage (2 - dosage) so the
    math comes out right. This happens more often than you'd think because
    different studies encode ref/alt differently.
    """
    # find overlapping variants between my VCF and the GWAS stats
    common_variants = set(genotype_data.keys()) & set(gwas_stats.keys())
    print(f"  Variants in VCF: {len(genotype_data)}")
    print(f"  Variants in GWAS (post-filter): {len(gwas_stats)}")
    print(f"  Overlapping variants: {len(common_variants)}")

    if len(common_variants) == 0:
        print("  WARNING: No overlapping variants found! Something is probably wrong.")
        print("  Check that your VCF rsIDs match the GWAS summary stats rsIDs.")

    # for missing genotypes, I'll use the mean dosage across all called samples
    # this is a simple mean imputation — not great but acceptable for a demo
    variant_mean_dosage = {}
    for rsid in common_variants:
        dosages = [d for d in genotype_data[rsid].values() if d is not None]
        variant_mean_dosage[rsid] = sum(dosages) / len(dosages) if dosages else 0.0

    # calculate PRS for each sample
    scores = OrderedDict()
    variant_contributions = []  # for logging which variants mattered most

    for sample in sample_names:
        prs = 0.0
        n_used = 0
        n_missing = 0

        for rsid in sorted(common_variants):
            gwas = gwas_stats[rsid]
            beta = gwas["beta"]
            effect_allele = gwas["effect_allele"]
            ref, alt = variant_alleles[rsid]

            # get the dosage (or impute if missing)
            dosage = genotype_data[rsid].get(sample)
            if dosage is None:
                dosage = variant_mean_dosage[rsid]
                n_missing += 1

            # handle allele flipping
            # my VCF counts ALT alleles. if the GWAS effect allele is my REF,
            # I need to flip: dosage becomes (2 - dosage)
            if effect_allele == ref:
                dosage = 2.0 - dosage
            elif effect_allele != alt:
                # alleles don't match at all — skip this variant
                # this could be a strand issue or a triallelic site
                continue

            prs += dosage * beta
            n_used += 1

        scores[sample] = {"prs": prs, "variants_used": n_used, "variants_imputed": n_missing}

    return scores, common_variants


def write_output(filepath, scores, sample_names):
    """
    Write out the final PRS scores. Each line is one patient with their
    raw score, the z-score (standardized), and their percentile ranking.

    The z-score is important for clinical interpretation — a raw PRS of 0.47
    doesn't mean anything by itself. But a z-score of +2.1 means that person
    is 2.1 standard deviations above the population mean, which puts them
    in roughly the top 2% of genetic risk.
    """
    # compute z-scores and percentiles
    all_scores = [scores[s]["prs"] for s in sample_names]
    mean_prs = sum(all_scores) / len(all_scores)
    variance = sum((s - mean_prs) ** 2 for s in all_scores) / len(all_scores)
    std_prs = math.sqrt(variance) if variance > 0 else 1.0

    # rank for percentiles
    sorted_scores = sorted(all_scores)

    with open(filepath, "w") as f:
        f.write("SAMPLE\tPRS_RAW\tPRS_ZSCORE\tPERCENTILE\tRISK_CATEGORY\tVARIANTS_USED\tVARIANTS_IMPUTED\n")

        for sample in sample_names:
            info = scores[sample]
            raw = info["prs"]
            zscore = (raw - mean_prs) / std_prs
            # simple percentile: what fraction of scores are below this one
            rank = sum(1 for s in sorted_scores if s <= raw)
            percentile = (rank / len(sorted_scores)) * 100

            # risk categories based on percentile
            # these thresholds are commonly used in PRS literature
            if percentile >= 95:
                category = "HIGH_RISK"
            elif percentile >= 80:
                category = "ELEVATED"
            elif percentile >= 20:
                category = "AVERAGE"
            else:
                category = "LOW_RISK"

            f.write(f"{sample}\t{raw:.6f}\t{zscore:.4f}\t{percentile:.1f}\t{category}\t{info['variants_used']}\t{info['variants_imputed']}\n")

    print(f"\n  PRS Results Summary:")
    print(f"    Mean PRS: {mean_prs:.6f}")
    print(f"    Std PRS:  {std_prs:.6f}")
    print(f"    Range:    [{min(all_scores):.6f}, {max(all_scores):.6f}]")


def main():
    args = parse_args()

    print("=" * 60)
    print("Polygenic Risk Score Calculation")
    print("=" * 60)

    # step 1: load GWAS summary stats with p-value filter
    print(f"\n1. Loading GWAS summary statistics (p < {args.p_threshold})...")
    gwas_stats = load_gwas_stats(args.gwas, args.p_threshold)

    # step 2: parse VCF genotypes
    print(f"\n2. Parsing VCF genotypes from {args.vcf}...")
    sample_names, genotype_data, variant_alleles = parse_vcf_genotypes(args.vcf)
    print(f"  {len(sample_names)} samples, {len(genotype_data)} variants")

    # step 3: compute PRS
    print(f"\n3. Computing PRS (weighted sum)...")
    scores, used_variants = calculate_scores(sample_names, genotype_data, variant_alleles, gwas_stats)

    # step 4: write results
    print(f"\n4. Writing results to {args.output}...")
    write_output(args.output, scores, sample_names)

    print(f"\nDone! Scores written to {args.output}")


if __name__ == "__main__":
    main()
