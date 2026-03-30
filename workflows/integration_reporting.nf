/*
 * Phase 6: Integration & Reporting
 *
 * Merges all annotations, generates master CSVs, HTML report, and figures.
 */

include { INTEGRATE_ANNOTATIONS } from '../modules/local/integrate_annotations'
include { GENERATE_REPORT }       from '../modules/local/generate_report'
include { GENERATE_FIGURES }      from '../modules/local/generate_figures'

workflow INTEGRATION_REPORTING {
    take:
    ch_substrates            // [sample_id, substrates_filtered.tsv]
    ch_substrates_unfiltered // [sample_id, substrates_all.tsv]
    ch_annotations           // [sample_id, tool_name, annotations.csv] (mixed)
    ch_enrichment            // enrichment_results.csv

    main:
    // Group annotations by sample_id
    ch_annotations_grouped = ch_annotations
        .groupTuple(by: 0)
        .map { sample_id, tool_names, files -> tuple(sample_id, files) }

    // Join substrates with their annotations
    ch_integrate_input = ch_substrates
        .join(ch_substrates_unfiltered)
        .join(ch_annotations_grouped, remainder: true)

    INTEGRATE_ANNOTATIONS(ch_integrate_input)

    // Collect all integrated CSVs for reporting
    ch_all_integrated = INTEGRATE_ANNOTATIONS.out.master_csv.collect()

    GENERATE_REPORT(ch_all_integrated, ch_enrichment)
    GENERATE_FIGURES(ch_all_integrated)

    emit:
    master_csv = INTEGRATE_ANNOTATIONS.out.master_csv
    report     = GENERATE_REPORT.out.report
    figures    = GENERATE_FIGURES.out.figures
}
