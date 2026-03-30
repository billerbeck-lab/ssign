/*
 * CROSS_VALIDATE — merge and cross-validate prediction tool outputs
 *
 * Flags unreliable predictions (e.g., DeepSecE T3SS without MacSyFinder T3SS).
 */

process CROSS_VALIDATE {
    tag "$sample_id"
    label 'process_low'
    // container: ssign-base

    input:
    tuple val(sample_id), path(deeplocpro), path(deepsece), path(valid_systems)

    output:
    tuple val(sample_id), path("${sample_id}_predictions_validated.tsv"), emit: validated

    script:
    """
    cross_validate_predictions.py \\
        --deeplocpro ${deeplocpro} \\
        --deepsece ${deepsece} \\
        --valid-systems ${valid_systems} \\
        --sample ${sample_id} \\
        --conf-threshold ${params.conf_threshold} \\
        --output ${sample_id}_predictions_validated.tsv
    """
}
