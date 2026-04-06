/*
 * modules/database.nf — Store Results in SQLite
 *
 * flat files are fine for a demo, but in a clinical setting you need
 * a real database. I'm using SQLite here (single file, no server),
 * but the same schema works on PostgreSQL (Amazon RDS) or you could
 * put the variant weights in DynamoDB for fast key-value lookups.
 *
 * the database has three tables:
 *   - samples: patient scores (what a clinician queries)
 *   - variants: GWAS weights (the risk model)
 *   - qc_metrics: per-sample quality stats
 */

process STORE_RESULTS {
    tag "Database"
    publishDir "${params.outdir}/database", mode: 'copy'

    input:
    path prs_scores
    path gwas_stats
    path qc_stats

    output:
    path "prs_results.db",  emit: database

    script:
    """
    python3 "${projectDir}/bin/store_results.py" \\
        --scores ${prs_scores} \\
        --gwas ${gwas_stats} \\
        --qc-stats ${qc_stats} \\
        --db prs_results.db
    """
}
