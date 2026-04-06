/*
 * modules/score.nf — PRS Calculation
 *
 * this is the actual scoring step. takes the cleaned, imputed genotypes
 * and the GWAS summary statistics, and computes a weighted sum:
 *
 *   PRS = sum(dosage_i * beta_i) for all variants i
 *
 * the calculate_prs.py script handles:
 *   - matching variants between VCF and GWAS stats
 *   - allele flipping (when ref/alt are swapped)
 *   - p-value filtering (only use variants below the threshold)
 *   - z-score normalization of the final scores
 *   - risk categorization (low / average / elevated / high)
 */

process CALCULATE_PRS {
    tag "Score"
    publishDir "${params.outdir}/scores", mode: 'copy'

    input:
    path vcf
    path gwas_stats

    output:
    path "prs_scores.tsv",  emit: prs_scores

    script:
    """
    python3 "${projectDir}/bin/calculate_prs.py" \\
        --vcf ${vcf} \\
        --gwas ${gwas_stats} \\
        --output prs_scores.tsv \\
        --p-threshold ${params.p_threshold}
    """
}
