/*
 * BLASTP — BLASTp annotation in local or remote mode
 *
 * Bug fix preserved: hit_desc.split(" >")[0] before filtering
 * (NCBI concatenates multiple hit descriptions with ">")
 */

process BLASTP_LOCAL {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('blastp'), path("${sample_id}_blastp.csv"), emit: annotations

    script:
    def exclude_taxid = params.blastp_exclude_taxid ? "-negative_taxids ${params.blastp_exclude_taxid}" : ""
    """
    run_blastp.py \\
        --mode local \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --db "${params.blastp_db}" \\
        --evalue ${params.blastp_evalue} \\
        --min-pident ${params.blastp_min_pident} \\
        --min-qcov ${params.blastp_min_qcov} \\
        ${exclude_taxid ? "--exclude-taxid ${params.blastp_exclude_taxid}" : ""} \\
        --output ${sample_id}_blastp.csv
    """
}

process BLASTP_REMOTE {
    tag "$sample_id"
    label 'process_network'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('blastp'), path("${sample_id}_blastp.csv"), emit: annotations

    script:
    """
    run_blastp.py \\
        --mode remote \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --evalue ${params.blastp_evalue} \\
        --min-pident ${params.blastp_min_pident} \\
        --min-qcov ${params.blastp_min_qcov} \\
        ${params.blastp_exclude_taxid ? "--exclude-taxid ${params.blastp_exclude_taxid}" : ""} \\
        --output ${sample_id}_blastp.csv
    """
}
