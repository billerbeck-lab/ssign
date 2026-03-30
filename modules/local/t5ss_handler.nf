/*
 * T5SS_HANDLER — handle self-secreting autotransporters
 *
 * T5aSS components ARE substrates (self-secreting).
 * Uses PF03797 as ground truth, NOT DeepLocPro threshold.
 * Classifies into: Classical AT, Minimal passenger, Barrel-only, OMP/Porin.
 */

process T5SS_HANDLER {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(ss_components), path(predictions)

    output:
    tuple val(sample_id), path("${sample_id}_t5ss_substrates.tsv"), emit: substrates
    tuple val(sample_id), path("${sample_id}_t5ss_domains.tsv"),    emit: domains

    script:
    """
    t5ss_handler.py \\
        --ss-components ${ss_components} \\
        --predictions ${predictions} \\
        --sample ${sample_id} \\
        --out-substrates ${sample_id}_t5ss_substrates.tsv \\
        --out-domains ${sample_id}_t5ss_domains.tsv
    """
}
