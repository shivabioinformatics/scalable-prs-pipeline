#!/usr/bin/env python3
"""
generate_test_data.py
---------------------
I need simulated data to demo this pipeline without needing real patient genomes.
This script creates:
  - A fake multi-sample VCF with ~500 SNPs across 20 samples
  - GWAS summary statistics with effect sizes (betas) for those SNPs
  - A tiny reference panel VCF for the imputation step
  - A placeholder chain file for the liftover step

None of this is real patient data — it's all randomly generated so I can
show the pipeline runs end-to-end without needing IRB approval or huge files.
"""

import os
import random
import math

# so results are reproducible every time I regenerate
random.seed(42)

# ----- configuration -----
NUM_VARIANTS = 500        # how many SNPs to simulate
NUM_SAMPLES = 20          # how many "patients"
NUM_REF_SAMPLES = 50      # reference panel is bigger so imputation has something to work with
CHROMOSOMES = list(range(1, 23))  # autosomes only, keeping it simple
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def make_output_dir():
    """just make sure the data directory exists before I start writing files"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def generate_variant_positions():
    """
    Pick random positions across chromosomes.
    I'm spacing them out so they look semi-realistic — not just random noise
    clustered in one spot.
    """
    variants = []
    for i in range(NUM_VARIANTS):
        chrom = random.choice(CHROMOSOMES)
        # spread positions across a realistic range (1 to 250 million)
        pos = random.randint(10000, 250000000)
        ref_allele = random.choice(["A", "C", "G", "T"])
        # alt allele has to be different from the ref — obviously
        alt_options = [a for a in ["A", "C", "G", "T"] if a != ref_allele]
        alt_allele = random.choice(alt_options)
        rsid = f"rs{100000 + i}"
        variants.append((chrom, pos, rsid, ref_allele, alt_allele))

    # sort by chromosome then position so the VCF is properly ordered
    variants.sort(key=lambda x: (x[0], x[1]))
    return variants


def generate_genotype():
    """
    Give me a random genotype: 0/0, 0/1, or 1/1
    I'm weighting towards 0/0 (homozygous ref) since most variants in a
    real genome are going to be reference. This makes the PRS distribution
    look more realistic.
    """
    r = random.random()
    if r < 0.60:
        return "0/0"
    elif r < 0.85:
        return "0/1"
    else:
        return "1/1"


def generate_genotype_with_missing():
    """
    Same as above but ~5% of calls are missing (./.)
    Real sequencing data always has some missing calls — this gives
    the imputation step something to actually do.
    """
    if random.random() < 0.05:
        return "./."
    return generate_genotype()


def write_vcf(filepath, variants, sample_names, include_missing=False):
    """
    Write a valid VCF file. I kept the header minimal but correct —
    enough for downstream tools to parse it without complaining.
    """
    geno_func = generate_genotype_with_missing if include_missing else generate_genotype

    with open(filepath, "w") as f:
        # VCF header — the bare minimum that bcftools and beagle expect
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=generate_test_data.py\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')

        # contig lines — tools like bcftools need these
        for c in CHROMOSOMES:
            f.write(f"##contig=<ID=chr{c},length=250000000>\n")

        # column header
        header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        header_cols.extend(sample_names)
        f.write("\t".join(header_cols) + "\n")

        # variant lines
        for chrom, pos, rsid, ref, alt in variants:
            # calculate allele frequency from the genotypes i'm about to generate
            genotypes = []
            alt_count = 0
            total_alleles = 0
            for _ in sample_names:
                gt = geno_func()
                genotypes.append(gt)
                if gt != "./.":
                    alleles = gt.split("/")
                    alt_count += alleles.count("1")
                    total_alleles += 2

            af = alt_count / total_alleles if total_alleles > 0 else 0.0

            info = f"AF={af:.4f}"
            fmt = "GT:DP"

            # add depth values — random but in a realistic range
            geno_fields = []
            for gt in genotypes:
                dp = random.randint(15, 60) if gt != "./." else 0
                geno_fields.append(f"{gt}:{dp}")

            line_parts = [f"chr{chrom}", str(pos), rsid, ref, alt, ".", "PASS", info, fmt]
            line_parts.extend(geno_fields)
            f.write("\t".join(line_parts) + "\n")


def write_gwas_summary_stats(filepath, variants):
    """
    Create fake GWAS summary statistics.
    In a real project these would come from a published GWAS — like the
    CARDIoGRAMplusC4D consortium for coronary artery disease.

    Each variant gets a beta (effect size) and a p-value.
    Most variants have tiny effects and high p-values (not significant).
    A handful have larger effects and low p-values — these are the ones
    that actually drive the PRS.
    """
    with open(filepath, "w") as f:
        f.write("SNP\tCHR\tPOS\tEFFECT_ALLELE\tOTHER_ALLELE\tBETA\tSE\tP\tAF\n")

        for chrom, pos, rsid, ref, alt in variants:
            # most variants: small effect, bad p-value
            # ~10% of variants: moderate effect, decent p-value
            # this mimics what a real GWAS looks like
            if random.random() < 0.10:
                # "significant" hit — these drive the PRS
                beta = random.gauss(0, 0.15)
                se = abs(beta) / random.uniform(2.5, 5.0)
                p = random.uniform(1e-8, 5e-4)
            else:
                # noise — these barely contribute
                beta = random.gauss(0, 0.02)
                se = abs(beta) / random.uniform(0.5, 1.5) if abs(beta) > 0.001 else 0.01
                p = random.uniform(0.01, 1.0)

            af = random.uniform(0.05, 0.50)
            f.write(f"{rsid}\t{chrom}\t{pos}\t{alt}\t{ref}\t{beta:.6f}\t{se:.6f}\t{p:.2e}\t{af:.4f}\n")


def write_chain_file(filepath):
    """
    A chain file tells CrossMap how to convert coordinates between genome builds.
    This is a placeholder — real chain files are huge and come from UCSC.
    For the demo, I just need the file to exist so the pipeline doesn't crash.
    """
    with open(filepath, "w") as f:
        f.write("# Placeholder chain file for hg19 -> hg38 liftover\n")
        f.write("# In production, download from:\n")
        f.write("# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz\n")
        f.write("#\n")
        # write a minimal valid chain entry so the file isn't completely empty
        for c in CHROMOSOMES:
            f.write(f"chain 1000 chr{c} 250000000 + 0 250000000 chr{c} 250000000 + 0 250000000 {c}\n")
            f.write("250000000\n\n")


def main():
    print("=" * 60)
    print("Generating simulated test data for PRS pipeline")
    print("=" * 60)

    make_output_dir()
    variants = generate_variant_positions()

    # 1. Patient VCF — this is what a clinic would send us
    sample_names = [f"PATIENT_{i:03d}" for i in range(1, NUM_SAMPLES + 1)]
    vcf_path = os.path.join(OUTPUT_DIR, "sample_input.vcf")
    write_vcf(vcf_path, variants, sample_names, include_missing=True)
    print(f"  [OK] Patient VCF: {vcf_path}")
    print(f"       {NUM_VARIANTS} variants x {NUM_SAMPLES} samples, ~5% missing calls")

    # 2. GWAS summary statistics — the "risk model"
    gwas_path = os.path.join(OUTPUT_DIR, "gwas_summary_stats.tsv")
    write_gwas_summary_stats(gwas_path, variants)
    print(f"  [OK] GWAS summary stats: {gwas_path}")
    print(f"       {NUM_VARIANTS} variants with effect sizes and p-values")

    # 3. Reference panel — used by Beagle for imputation
    ref_names = [f"REF_{i:03d}" for i in range(1, NUM_REF_SAMPLES + 1)]
    ref_path = os.path.join(OUTPUT_DIR, "reference_panel.vcf")
    write_vcf(ref_path, variants, ref_names, include_missing=False)
    print(f"  [OK] Reference panel: {ref_path}")
    print(f"       {NUM_VARIANTS} variants x {NUM_REF_SAMPLES} reference samples")

    # 4. Chain file — for coordinate liftover
    chain_path = os.path.join(OUTPUT_DIR, "chain_file.chain")
    write_chain_file(chain_path)
    print(f"  [OK] Chain file: {chain_path}")

    print()
    print("All test data generated successfully!")
    print(f"Files are in: {OUTPUT_DIR}")
    print()


if __name__ == "__main__":
    main()
