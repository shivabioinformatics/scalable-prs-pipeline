/*
 * modules/visualize.nf — Generate Figures
 *
 * five plots that tell the story of the analysis:
 *   1. QC summary — did the data pass quality checks?
 *   2. PRS distribution — bell curve with risk thresholds
 *   3. Risk stratification — who's in each risk group?
 *   4. Variant effects — Manhattan-style plot of effect sizes
 *   5. Score comparison — ranked lollipop chart per patient
 *
 * these are the kinds of figures you'd put in a clinical report
 * or a paper. matplotlib makes them publication-quality.
 */

process VISUALIZE_PRS {
    tag "Visualize"
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path prs_scores
    path qc_stats
    path gwas_stats

    output:
    path "*.png",           emit: figures

    script:
    """
    python3 "${projectDir}/bin/visualize_prs.py" \\
        --scores ${prs_scores} \\
        --qc-stats ${qc_stats} \\
        --gwas ${gwas_stats} \\
        --outdir .
    """
}
