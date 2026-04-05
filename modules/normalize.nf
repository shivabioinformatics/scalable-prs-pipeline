/*
 * modules/normalize.nf — Genome Build Normalization
 *
 * this is the step a lot of people skip and then wonder why their
 * results look weird. the problem: my patient VCF might be on hg19
 * but the GWAS summary stats are on hg38 (or vice versa). if I don't
 * convert them to the same build, the positions won't line up and
 * I'll miss most of the variants.
 *
 * in production: CrossMap for coordinate liftover + bcftools norm
 * for the demo: I simulate the normalization (data is already on same build)
 *
 * I also want to talk about batch normalization here since the feedback
 * mentioned it. Batch normalization in PRS context means adjusting the
 * raw scores so they're comparable across different genotyping batches
 * or ancestry groups. A raw PRS of 0.5 means different things for
 * European vs African ancestry populations because allele frequencies
 * differ. Z-score normalization within ancestry groups fixes this.
 */

process NORMALIZE_VCF {
    tag "Normalize"
    publishDir "${params.outdir}/normalized", mode: 'copy'

    input:
    path vcf
    path chain_file

    output:
    path "normalized.vcf",  emit: normalized_vcf

    script:
    /*
     * Steps:
     * 1. CrossMap liftover — converts coordinates from one build to another
     *    using the chain file (maps old positions to new positions)
     * 2. bcftools norm — left-aligns indels and splits multiallelic sites
     *    so every variant has exactly one ref and one alt allele
     *
     * For the demo, the simulated data is already on the right build,
     * so I just copy it through. The commands below show what it would
     * look like in production.
     *
     * Production commands (commented out):
     *   CrossMap.py vcf chain_file input.vcf hg38.fa normalized_raw.vcf
     *   bcftools norm -f hg38.fa -m -both normalized_raw.vcf > normalized.vcf
     *
     * Batch normalization of PRS scores (done in the scoring step):
     *   z_score = (raw_prs - population_mean) / population_std
     *   This makes scores comparable across ancestry groups
     */
    """
    echo "Normalizing VCF to ${params.genome_build}..."
    echo "Chain file: ${chain_file}"

    # for the demo, the data is already normalized
    # in production, this would run CrossMap + bcftools norm
    cp ${vcf} normalized.vcf

    echo "Normalization complete"
    echo "Note: In production, this step would run:"
    echo "  1. CrossMap.py vcf <chain> <input> <ref.fa> <output>"
    echo "  2. bcftools norm -f <ref.fa> -m -both <input> > <output>"
    """
}
