/*
 * modules/qc.nf — Quality Control for the input VCF
 *
 * this is step 1 of everything. if my input data is bad, the whole
 * PRS is meaningless. I'm checking call rates, allele frequencies,
 * and read depth. variants that fail get dropped.
 *
 * in production I'd use bcftools stats + bcftools filter here.
 * for the demo I wrote a python script that does the same checks
 * so it runs without needing bcftools installed.
 */

process QC_VCF {
    tag "QC"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path vcf

    output:
    path "filtered.vcf",    emit: filtered_vcf
    path "qc_stats.tsv",    emit: qc_stats

    script:
    /*
     * the qc_filter.py script checks three things per variant:
     *   1. call rate >= 0.95 (at least 95% of samples have this variant called)
     *   2. MAF >= 0.01 (minor allele freq — super rare variants aren't useful for PRS)
     *   3. mean depth >= 10 (low depth = low confidence genotypes)
     *
     * it also checks per-sample call rate and flags bad samples
     */
    """
    python3 "${projectDir}/bin/qc_filter.py" \\
        --vcf ${vcf} \\
        --out-vcf filtered.vcf \\
        --out-stats qc_stats.tsv \\
        --min-call-rate ${params.min_call_rate} \\
        --min-maf ${params.min_maf} \\
        --min-depth ${params.min_depth} \\
        --min-sample-call-rate ${params.min_sample_call_rate}
    """
}
