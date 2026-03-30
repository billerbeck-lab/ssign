/*
 * INTERPROSCAN — domain annotation (local or EBI REST API)
 */

process INTERPROSCAN_LOCAL {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('interproscan'), path("${sample_id}_interproscan.csv"), emit: annotations

    script:
    """
    run_interproscan.py \\
        --mode local \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --db "${params.interproscan_db}" \\
        --output ${sample_id}_interproscan.csv
    """
}

process INTERPROSCAN_REMOTE {
    tag "$sample_id"
    label 'process_network'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('interproscan'), path("${sample_id}_interproscan.csv"), emit: annotations

    script:
    """
    run_interproscan.py \\
        --mode remote \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --output ${sample_id}_interproscan.csv
    """
}
