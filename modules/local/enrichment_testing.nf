/*
 * ENRICHMENT_TESTING — statistical enrichment of SS types near extracellular proteins
 */

process ENRICHMENT_TESTING {
    tag "enrichment"
    label 'process_low'

    input:
    path(substrate_files)  // collected from all genomes

    output:
    path("enrichment_results.csv"), emit: results
    path("enrichment_summary.txt"), emit: summary

    script:
    """
    enrichment_testing.py \\
        --substrate-files ${substrate_files} \\
        --out-csv enrichment_results.csv \\
        --out-summary enrichment_summary.txt
    """
}
