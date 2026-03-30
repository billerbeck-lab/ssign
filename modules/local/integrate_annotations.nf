/*
 * INTEGRATE_ANNOTATIONS — merge all annotation sources into master CSV
 */

process INTEGRATE_ANNOTATIONS {
    tag "$sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(substrates_filtered), path(substrates_all),
          path(annotation_files)

    output:
    tuple val(sample_id), path("${sample_id}_master.csv"), emit: master_csv

    script:
    """
    integrate_annotations.py \\
        --substrates-filtered ${substrates_filtered} \\
        --substrates-all ${substrates_all} \\
        --annotations ${annotation_files} \\
        --sample ${sample_id} \\
        --output ${sample_id}_master.csv
    """
}
