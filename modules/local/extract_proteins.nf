/*
 * EXTRACT_PROTEINS — extract protein sequences and gene info from annotated files
 *
 * Accepts GenBank or GFF3(+FASTA). Files are staged via path().
 */

process EXTRACT_PROTEINS {
    tag "$sample_id"
    label 'process_low'
    // container: ssign-base

    input:
    tuple val(sample_id), path(input_files)

    output:
    tuple val(sample_id), path("${sample_id}_proteins.faa"), emit: proteins
    tuple val(sample_id), path("${sample_id}_gene_info.tsv"), emit: gene_info

    script:
    if (input_files instanceof List && input_files.size() > 1) {
        // Two files: GFF3 + FASTA
        """
        extract_proteins.py \\
            --input ${input_files[0]} \\
            --fasta ${input_files[1]} \\
            --sample ${sample_id} \\
            --out-proteins ${sample_id}_proteins.faa \\
            --out-gene-info ${sample_id}_gene_info.tsv
        """
    } else {
        // Single file: GenBank
        def the_file = input_files instanceof List ? input_files[0] : input_files
        """
        extract_proteins.py \\
            --input ${the_file} \\
            --sample ${sample_id} \\
            --out-proteins ${sample_id}_proteins.faa \\
            --out-gene-info ${sample_id}_gene_info.tsv
        """
    }
}
