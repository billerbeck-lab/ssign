/*
 * PRODIGAL — predict protein-coding genes from unannotated contigs
 */

process PRODIGAL {
    tag "$sample_id"
    label 'process_medium'
    // container: ssign-annotation

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_proteins.faa"), emit: proteins
    tuple val(sample_id), path("${sample_id}_gene_info.tsv"), emit: gene_info

    script:
    """
    # Run Prodigal in metagenomic mode
    prodigal \\
        -i ${contigs} \\
        -a ${sample_id}_proteins_raw.faa \\
        -o ${sample_id}_genes.gff \\
        -f gff \\
        -p meta

    # Convert Prodigal output to our standard gene_info format
    prodigal_to_gene_info.py \\
        --proteins ${sample_id}_proteins_raw.faa \\
        --gff ${sample_id}_genes.gff \\
        --sample ${sample_id} \\
        --out-proteins ${sample_id}_proteins.faa \\
        --out-gene-info ${sample_id}_gene_info.tsv
    """
}
