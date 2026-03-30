/*
 * Phase 1: Input Processing
 *
 * Detects input format (FASTA contigs, GenBank, GFF3), runs ORF prediction
 * if needed, extracts proteins and gene order information.
 */

include { DETECT_FORMAT }      from '../modules/local/detect_format'
include { PRODIGAL }           from '../modules/local/prodigal'
include { EXTRACT_PROTEINS }   from '../modules/local/extract_proteins'
include { EXTRACT_GENE_ORDER } from '../modules/local/extract_gene_order'

workflow INPUT_PROCESSING {
    take:
    ch_input  // channel: [sample_id, [file(s)]]

    main:
    // Detect input format for each sample
    DETECT_FORMAT(ch_input)

    // Branch by detected format
    DETECT_FORMAT.out.result
        .branch {
            genbank: it[2] == 'genbank'
            gff3:    it[2] == 'gff3'
            contigs: it[2] == 'fasta_contigs'
        }
        .set { ch_branched }

    // Raw contigs → Prodigal for ORF prediction
    // Extract first file from the list
    ch_contigs = ch_branched.contigs.map { sample, files, fmt ->
        tuple(sample, files[0])
    }
    PRODIGAL(ch_contigs)

    // GenBank/GFF3 → extract proteins directly
    // Pass files as-is (extract_proteins.py handles both single and paired)
    ch_annotated = ch_branched.genbank
        .map { sample, files, fmt -> tuple(sample, files) }
        .mix(
            ch_branched.gff3.map { sample, files, fmt -> tuple(sample, files) }
        )
    EXTRACT_PROTEINS(ch_annotated)

    // Merge protein outputs from both paths
    ch_proteins = EXTRACT_PROTEINS.out.proteins
        .mix(PRODIGAL.out.proteins)

    ch_gene_info = EXTRACT_PROTEINS.out.gene_info
        .mix(PRODIGAL.out.gene_info)

    // Extract gene order for proximity analysis
    EXTRACT_GENE_ORDER(ch_gene_info)

    emit:
    proteins    = ch_proteins                    // [sample_id, proteins.faa]
    gene_info   = ch_gene_info                   // [sample_id, gene_info.tsv]
    gene_order  = EXTRACT_GENE_ORDER.out.order   // [sample_id, gene_order.tsv]
}
