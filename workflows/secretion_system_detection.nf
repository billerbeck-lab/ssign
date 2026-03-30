/*
 * Phase 2: Secretion System Detection
 *
 * Runs MacSyFinder v2 with TXSScan models to identify secretion systems,
 * then validates systems by wholeness threshold.
 */

include { MACSYFINDER }      from '../modules/local/macsyfinder'
include { VALIDATE_SYSTEMS } from '../modules/local/validate_systems'

workflow SECRETION_SYSTEM_DETECTION {
    take:
    ch_proteins   // [sample_id, proteins.faa]
    ch_gene_info  // [sample_id, gene_info.tsv]

    main:
    MACSYFINDER(ch_proteins)

    // Join MacSyFinder results with gene_info for component mapping
    ch_msf_with_info = MACSYFINDER.out.results
        .join(ch_gene_info)

    VALIDATE_SYSTEMS(ch_msf_with_info)

    emit:
    ss_components = VALIDATE_SYSTEMS.out.components  // [sample_id, ss_components.tsv]
    valid_systems = VALIDATE_SYSTEMS.out.systems     // [sample_id, valid_systems.tsv]
}
