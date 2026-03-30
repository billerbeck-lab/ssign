/*
 * DETECT_FORMAT — identify input file type (GenBank, GFF3, or raw contigs)
 *
 * Accepts one or two files. Stages them as path() so they're available
 * in containerized environments. Outputs the format as an env variable
 * and passes files through.
 */

process DETECT_FORMAT {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(input_files)

    output:
    tuple val(sample_id), path(input_files), env(FORMAT), emit: result

    script:
    """
    FORMAT=\$(detect_input_format.py "${input_files[0]}")
    """
}
