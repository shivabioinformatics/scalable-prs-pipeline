#!/usr/bin/env nextflow

/*
 * main.nf — Entry point for the PRS pipeline
 *
 * gotta have DSL2 turned on for the modular setup to work.
 * everything is split into separate modules and a sub-workflow
 * so each step is independently testable and swappable.
 *
 * the basic flow is:
 *   VCF in -> QC -> normalize -> impute -> score -> store -> visualize
 *
 * I'm using channels to pass data between steps — that's the Nextflow way
 * of keeping things parallel and reproducible.
 */

nextflow.enable.dsl = 2

// ---- pull in my modules ----
include { QC_VCF }              from './modules/qc'
include { NORMALIZE_VCF }       from './modules/normalize'
include { IMPUTE_GENOTYPES }    from './modules/impute'
include { CALCULATE_PRS }       from './modules/score'
include { STORE_RESULTS }       from './modules/database'
include { VISUALIZE_PRS }       from './modules/visualize'

// ---- pull in the sub-workflow that chains everything ----
include { PRS_WORKFLOW }        from './subworkflows/prs_workflow'

/*
 * print a banner so I know the pipeline actually started
 * and which parameters it's using
 */
log.info """
============================================================
  POLYGENIC RISK SCORE (PRS) PIPELINE
============================================================
  Input VCF       : ${params.input_vcf}
  GWAS Stats      : ${params.gwas_stats}
  Reference Panel : ${params.ref_panel}
  Chain File      : ${params.chain_file}
  Output Dir      : ${params.outdir}
  P-value Cutoff  : ${params.p_threshold}
  Min Call Rate   : ${params.min_call_rate}
  Min MAF         : ${params.min_maf}
  Genome Build    : ${params.genome_build}
============================================================
"""

/*
 * the main workflow — just calls the sub-workflow
 * I keep it separate so I could add more sub-workflows later
 * (like a different disease model or a validation track)
 */
workflow {
    // create input channels from the parameters
    vcf_ch    = Channel.fromPath(params.input_vcf, checkIfExists: true)
    gwas_ch   = Channel.fromPath(params.gwas_stats, checkIfExists: true)
    ref_ch    = Channel.fromPath(params.ref_panel, checkIfExists: true)
    chain_ch  = Channel.fromPath(params.chain_file, checkIfExists: true)

    // run the full PRS pipeline
    PRS_WORKFLOW(vcf_ch, gwas_ch, ref_ch, chain_ch)
}

/*
 * what happens when the pipeline finishes
 * just a nice message so I know it completed and where to find results
 */
workflow.onComplete {
    log.info """
    ============================================================
    Pipeline Complete!
    Results are in: ${params.outdir}
    ============================================================
    """.stripIndent()
}
