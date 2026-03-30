/*
 * DEEPSECE — predict secretion system substrate type (MIT license)
 */

process DEEPSECE {
    tag "$sample_id"
    label 'process_high'
    // container: ssign-base (DeepSecE is pip-installable)

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_deepsece.tsv"), emit: predictions

    script:
    """
    run_deepsece.py \\
        --input ${proteins} \\
        --sample ${sample_id} \\
        --output ${sample_id}_deepsece.tsv
    """
}
