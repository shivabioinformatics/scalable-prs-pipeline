/*
 * modules/impute.nf — Genotype Imputation
 *
 * real sequencing data never covers every single base. there are always
 * gaps — maybe the read depth was too low at some positions, or the
 * variant caller wasn't confident enough to make a call. but my PRS
 * model expects genotypes at specific positions.
 *
 * imputation fills in those gaps by using patterns of linkage
 * disequilibrium (LD) from a reference panel. basically, if you know
 * someone has allele A at position 1000 and allele G at position 3000,
 * and in the reference panel those two alleles are almost always paired
 * with allele T at position 2000, then you can infer position 2000
 * is probably T even if you didn't sequence it directly.
 *
 * Beagle is the standard tool for this. it's fast and accurate.
 * in production this is the most compute-heavy step — that's why
 * I have 8 CPUs and 16GB allocated, and why I use spot instances on AWS.
 */

process IMPUTE_GENOTYPES {
    tag "Impute"
    publishDir "${params.outdir}/imputed", mode: 'copy'

    input:
    path vcf
    path ref_panel

    output:
    path "imputed.vcf",     emit: imputed_vcf

    script:
    /*
     * Beagle takes:
     *   - gt: the target genotypes (my patient VCF with missing calls)
     *   - ref: the reference panel (1000 Genomes, etc.)
     *   - out: output prefix
     *
     * After imputation, I filter on DR2 (dosage R-squared) > 0.3
     * DR2 measures imputation quality — below 0.3 means the imputation
     * wasn't confident and the genotype is probably wrong.
     *
     * Production command:
     *   java -jar beagle.jar gt=input.vcf ref=panel.vcf out=imputed
     *   bcftools filter -i 'INFO/DR2 > 0.3' imputed.vcf.gz > imputed_filtered.vcf
     *
     * For the demo, I do a simple mean-genotype imputation of missing calls
     * using a python one-liner. Not as good as Beagle, but shows the concept.
     */
    """
    echo "Imputing missing genotypes..."
    echo "Reference panel: ${ref_panel}"
    echo "Note: In production, this runs Beagle with a 1000 Genomes reference panel"

    # demo imputation: replace ./. with 0/0 (simplest possible imputation)
    # in production Beagle does this properly using LD patterns
    sed 's|\\./\\.|0/0|g' ${vcf} > imputed.vcf

    # count how many sites were imputed
    TOTAL=\$(grep -v "^#" ${vcf} | wc -l)
    MISSING=\$(grep -v "^#" ${vcf} | grep -c "\\./\\." || true)
    echo "Imputed \${MISSING} missing calls across \${TOTAL} variant sites"
    echo ""
    echo "Production workflow:"
    echo "  1. java -jar beagle.jar gt=<input> ref=<panel> out=<output>"
    echo "  2. Filter: bcftools filter -i 'INFO/DR2 > 0.3' <output>"
    """
}
