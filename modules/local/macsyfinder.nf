/*
 * MACSYFINDER — detect secretion systems using TXSScan models
 */

process MACSYFINDER {
    tag "$sample_id"
    label 'process_medium'
    // container: ssign-annotation

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("macsyfinder_out"), emit: results

    script:
    """
    macsyfinder \\
        --sequence-db ${proteins} \\
        --db-type ordered_replicon \\
        --models TXSScan all \\
        --out-dir macsyfinder_out \\
        --worker ${task.cpus}
    """
}
