/*
 * EXTRACT_GENE_ORDER — produce contig-sorted gene order for proximity analysis
 */

process EXTRACT_GENE_ORDER {
    tag "$sample_id"
    label 'process_low'
    // container: ssign-base

    input:
    tuple val(sample_id), path(gene_info)

    output:
    tuple val(sample_id), path("${sample_id}_gene_order.tsv"), emit: order

    script:
    """
    extract_gene_order.py \\
        --gene-info ${gene_info} \\
        --output ${sample_id}_gene_order.tsv
    """
}
