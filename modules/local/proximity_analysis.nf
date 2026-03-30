/*
 * PROXIMITY_ANALYSIS — find extracellular proteins near SS components
 *
 * Uses PER-COMPONENT proximity (not system-boundary).
 * Never spans contigs.
 */

process PROXIMITY_ANALYSIS {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(gene_order), path(ss_components), path(predictions)

    output:
    tuple val(sample_id), path("${sample_id}_proximity_substrates.tsv"), emit: substrates

    script:
    """
    proximity_analysis.py \\
        --gene-order ${gene_order} \\
        --ss-components ${ss_components} \\
        --predictions ${predictions} \\
        --sample ${sample_id} \\
        --window ${params.proximity_window} \\
        --conf-threshold ${params.conf_threshold} \\
        --output ${sample_id}_proximity_substrates.tsv
    """
}
