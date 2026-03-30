/*
 * Phase 4: Substrate Identification & Filtering
 *
 * Per-component proximity analysis, T5SS special handling,
 * enrichment testing, and system-level filtering.
 */

include { PROXIMITY_ANALYSIS }  from '../modules/local/proximity_analysis'
include { T5SS_HANDLER }        from '../modules/local/t5ss_handler'
include { ENRICHMENT_TESTING }  from '../modules/local/enrichment_testing'
include { SYSTEM_FILTERING }    from '../modules/local/system_filtering'

workflow SUBSTRATE_IDENTIFICATION {
    take:
    ch_gene_order     // [sample_id, gene_order.tsv]
    ch_ss_components  // [sample_id, ss_components.tsv]
    ch_predictions    // [sample_id, predictions_validated.tsv]
    ch_valid_systems  // [sample_id, valid_systems.tsv]

    main:
    // Join all inputs by sample_id
    ch_proximity_input = ch_gene_order
        .join(ch_ss_components)
        .join(ch_predictions)

    // Per-component proximity analysis (NOT system-boundary)
    PROXIMITY_ANALYSIS(ch_proximity_input)

    // T5SS special handling (self-secreting autotransporters)
    ch_t5ss_input = ch_ss_components
        .join(ch_predictions)

    T5SS_HANDLER(ch_t5ss_input)

    // Merge proximity substrates with T5SS substrates, then filter
    ch_filter_input = PROXIMITY_ANALYSIS.out.substrates
        .join(T5SS_HANDLER.out.substrates)
        .join(ch_valid_systems)
        .join(ch_predictions)

    SYSTEM_FILTERING(ch_filter_input)

    // Collect all filtered substrates for enrichment testing
    ch_all_substrates = SYSTEM_FILTERING.out.filtered.collect()

    ENRICHMENT_TESTING(ch_all_substrates)

    emit:
    substrates            = SYSTEM_FILTERING.out.filtered     // [sample_id, substrates_filtered.tsv]
    substrates_unfiltered = SYSTEM_FILTERING.out.unfiltered   // [sample_id, substrates_all.tsv]
    t5ss_domains          = T5SS_HANDLER.out.domains          // [sample_id, t5ss_domains.tsv]
    enrichment            = ENRICHMENT_TESTING.out.results    // enrichment_results.csv
}
