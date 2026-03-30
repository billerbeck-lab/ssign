/*
 * DEEPLOCPRO — predict subcellular localization (DTU license required)
 */

process DEEPLOCPRO {
    tag "$sample_id"
    label 'process_high'
    // User must provide DeepLocPro via params.deeplocpro_path

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_deeplocpro.tsv"), emit: predictions

    script:
    """
    run_deeplocpro.py \\
        --input ${proteins} \\
        --sample ${sample_id} \\
        --deeplocpro-path "${params.deeplocpro_path}" \\
        --conf-threshold ${params.conf_threshold} \\
        --output ${sample_id}_deeplocpro.tsv
    """
}
