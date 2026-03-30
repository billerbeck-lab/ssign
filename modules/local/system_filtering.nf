/*
 * SYSTEM_FILTERING — apply exclusion filters and DSE cross-genome validation
 *
 * Critical bug fix: DSE-only substrates require predicted SS type to exist
 * in the same genome (dse_type_in_genome check).
 */

process SYSTEM_FILTERING {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(proximity_subs), path(t5ss_subs),
          path(valid_systems), path(predictions)

    output:
    tuple val(sample_id), path("${sample_id}_substrates_filtered.tsv"), emit: filtered
    tuple val(sample_id), path("${sample_id}_substrates_all.tsv"),      emit: unfiltered

    script:
    """
    system_filtering.py \\
        --proximity-substrates ${proximity_subs} \\
        --t5ss-substrates ${t5ss_subs} \\
        --valid-systems ${valid_systems} \\
        --predictions ${predictions} \\
        --sample ${sample_id} \\
        --excluded-systems "${params.excluded_systems}" \\
        --required-fraction-correct ${params.required_fraction_correct} \\
        --out-filtered ${sample_id}_substrates_filtered.tsv \\
        --out-all ${sample_id}_substrates_all.tsv
    """
}
