/*
 * subworkflows/prs_workflow.nf — Full PRS Pipeline Chain
 *
 * this is where all the modules get wired together in order.
 * each step's output becomes the next step's input via channels.
 *
 * the flow:
 *   VCF -> QC -> Normalize -> Impute -> Score -> Store -> Visualize
 *
 * I could've put this all in main.nf, but splitting it into a
 * sub-workflow means I can potentially add other workflows later
 * (like a validation workflow, or a different disease model) without
 * cluttering the main entry point.
 */

include { QC_VCF }           from '../modules/qc'
include { NORMALIZE_VCF }    from '../modules/normalize'
include { IMPUTE_GENOTYPES } from '../modules/impute'
include { CALCULATE_PRS }    from '../modules/score'
include { STORE_RESULTS }    from '../modules/database'
include { VISUALIZE_PRS }    from '../modules/visualize'

workflow PRS_WORKFLOW {
    take:
    vcf_ch          // patient VCF file
    gwas_ch         // GWAS summary statistics
    ref_panel_ch    // reference panel for imputation
    chain_ch        // chain file for liftover

    main:
    // step 1: quality control
    // filter out bad variants and flag bad samples
    QC_VCF(vcf_ch)

    // step 2: normalize to common genome build
    // make sure coordinates match between my data and the GWAS stats
    NORMALIZE_VCF(QC_VCF.out.filtered_vcf, chain_ch)

    // step 3: impute missing genotypes
    // fill in gaps using LD patterns from the reference panel
    IMPUTE_GENOTYPES(NORMALIZE_VCF.out.normalized_vcf, ref_panel_ch)

    // step 4: calculate PRS
    // the weighted sum — dosage * beta for each variant
    CALCULATE_PRS(IMPUTE_GENOTYPES.out.imputed_vcf, gwas_ch)

    // step 5: store in database
    // put scores, variants, and QC metrics into SQLite
    STORE_RESULTS(
        CALCULATE_PRS.out.prs_scores,
        gwas_ch,
        QC_VCF.out.qc_stats
    )

    // step 6: generate visualizations
    // 5 plots for the results report
    VISUALIZE_PRS(
        CALCULATE_PRS.out.prs_scores,
        QC_VCF.out.qc_stats,
        gwas_ch
    )

    emit:
    scores   = CALCULATE_PRS.out.prs_scores
    database = STORE_RESULTS.out.database
    figures  = VISUALIZE_PRS.out.figures
}
