/*
 * FOLDSEEK — structural homology search
 *
 * Uses QUAL-01 strategy: batch search first, per-protein fallback
 * with relaxed E-value on retry.
 * Uses qtmscore (query-normalized), NOT alntmscore.
 */

process FOLDSEEK_LOCAL {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('foldseek'), path("${sample_id}_foldseek.csv"), emit: annotations

    script:
    """
    run_foldseek.py \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --db "${params.foldseek_db}" \\
        --evalue ${params.foldseek_evalue} \\
        --tmscore ${params.foldseek_tmscore} \\
        --output ${sample_id}_foldseek.csv
    """
}
