/*
 * SIGNALP — predict signal peptides (DTU license required)
 */

process SIGNALP {
    tag "$sample_id"
    label 'process_medium'
    // User must provide SignalP via params.signalp_path

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_signalp.tsv"), emit: predictions

    when:
    !params.skip_signalp && params.signalp_path

    script:
    """
    run_signalp.py \\
        --input ${proteins} \\
        --sample ${sample_id} \\
        --signalp-path "${params.signalp_path}" \\
        --output ${sample_id}_signalp.tsv
    """
}
