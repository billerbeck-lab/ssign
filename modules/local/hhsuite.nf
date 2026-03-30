/*
 * HHSUITE — HH-suite annotation (local hhblits/hhsearch or MPI Toolkit API)
 *
 * CRITICAL: Remote mode uses "alignment" parameter, NOT "sequence"!
 */

process HHSUITE_LOCAL {
    tag "$sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('hhsuite'), path("${sample_id}_hhsuite.csv"), emit: annotations

    script:
    """
    run_hhsuite.py \\
        --mode local \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --pfam-db "${params.hhsuite_pfam_db}" \\
        --pdb70-db "${params.hhsuite_pdb70_db}" \\
        --uniclust-db "${params.hhsuite_uniclust_db}" \\
        --output ${sample_id}_hhsuite.csv
    """
}

process HHSUITE_REMOTE {
    tag "$sample_id"
    label 'process_network'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('hhsuite'), path("${sample_id}_hhsuite.csv"), emit: annotations

    script:
    """
    run_hhsuite.py \\
        --mode remote \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --output ${sample_id}_hhsuite.csv
    """
}
