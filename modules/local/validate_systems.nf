/*
 * VALIDATE_SYSTEMS — filter MacSyFinder results by wholeness and map components
 */

process VALIDATE_SYSTEMS {
    tag "$sample_id"
    label 'process_low'
    // container: ssign-base

    input:
    tuple val(sample_id), path(msf_dir), path(gene_info)

    output:
    tuple val(sample_id), path("${sample_id}_ss_components.tsv"), emit: components
    tuple val(sample_id), path("${sample_id}_valid_systems.tsv"), emit: systems

    script:
    """
    validate_macsyfinder_systems.py \\
        --msf-dir ${msf_dir} \\
        --gene-info ${gene_info} \\
        --sample ${sample_id} \\
        --wholeness-threshold ${params.wholeness_threshold} \\
        --excluded-systems "${params.excluded_systems}" \\
        --out-components ${sample_id}_ss_components.tsv \\
        --out-systems ${sample_id}_valid_systems.tsv
    """
}
