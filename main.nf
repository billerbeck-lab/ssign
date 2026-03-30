#!/usr/bin/env nextflow
/*
 * ssign — Secretion-system Identification for Gram Negatives
 *
 * A Nextflow DSL2 pipeline for identifying secretion system substrates
 * in gram-negative bacterial genomes.
 *
 * Pipeline phases:
 *   1. Input Processing     — format detection, ORF prediction, protein extraction
 *   2. SS Detection         — MacSyFinder v2 with TXSScan models
 *   3. Secreted Prediction  — DeepLocPro, DeepSecE, SignalP 6.0
 *   4. Substrate ID         — proximity analysis, T5SS handling, enrichment, filtering
 *   5. Annotation (opt.)    — BLASTp, HH-suite, Foldseek, InterProScan, ProtParam, etc.
 *   6. Integration          — merge annotations, generate report + figures
 *
 * Usage:
 *   nextflow run ssign --input genome.gbff --outdir results
 *   nextflow run ssign --input samplesheet.csv --outdir results -profile docker
 *
 * License: GPL-3.0
 */

nextflow.enable.dsl = 2

/*
 * ──────────────────────────────────────────────────────────────────────
 *  VALIDATE PARAMETERS
 * ──────────────────────────────────────────────────────────────────────
 */

if (!params.input) {
    error "ERROR: --input is required. Provide a genome file or samplesheet CSV."
}

/*
 * ──────────────────────────────────────────────────────────────────────
 *  LOG HEADER
 * ──────────────────────────────────────────────────────────────────────
 */

log.info """\
    ╔═══════════════════════════════════════════╗
    ║  ssign — SS Identification for Gram (-)   ║
    ╚═══════════════════════════════════════════╝

    input          : ${params.input}
    outdir         : ${params.outdir}
    excluded SS    : ${params.excluded_systems}
    skip annotation: ${params.skip_annotation}
    profile        : ${workflow.profile}
    ─────────────────────────────────────────────
    """.stripIndent()

/*
 * ──────────────────────────────────────────────────────────────────────
 *  INCLUDE SUBWORKFLOWS
 * ──────────────────────────────────────────────────────────────────────
 */

include { INPUT_PROCESSING }            from './workflows/input_processing'
include { SECRETION_SYSTEM_DETECTION }  from './workflows/secretion_system_detection'
include { SECRETED_PROTEIN_PREDICTION } from './workflows/secreted_protein_prediction'
include { SUBSTRATE_IDENTIFICATION }    from './workflows/substrate_identification'
include { OPTIONAL_ANNOTATION }         from './workflows/optional_annotation'
include { INTEGRATION_REPORTING }       from './workflows/integration_reporting'

/*
 * ──────────────────────────────────────────────────────────────────────
 *  HELPER: parse input — single file or samplesheet CSV
 * ──────────────────────────────────────────────────────────────────────
 */

def parse_input(input_path) {
    def input_file = file(input_path, checkIfExists: true)

    // If it's a CSV samplesheet, parse rows
    if (input_file.name.endsWith('.csv')) {
        return Channel
            .fromPath(input_path)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def sample = row.sample
                def files = []
                if (row.input_1) files << file(row.input_1, checkIfExists: true)
                if (row.containsKey('input_2') && row.input_2) {
                    files << file(row.input_2, checkIfExists: true)
                }
                return tuple(sample, files)
            }
    }

    // Single file — derive sample name from filename
    def sample_name = input_file.baseName
        .replaceAll(/_genomic$/, '')
        .replaceAll(/\.(gbff|gbk|gb|gff3?|fasta|fna|fa)$/, '')
    return Channel.of(tuple(sample_name, [input_file]))
}

/*
 * ──────────────────────────────────────────────────────────────────────
 *  MAIN WORKFLOW
 * ──────────────────────────────────────────────────────────────────────
 */

workflow {
    // Parse input into channel of [sample_id, [file(s)]]
    ch_input = parse_input(params.input)

    // Phase 1: Input Processing
    INPUT_PROCESSING(ch_input)

    // Phase 2: Secretion System Detection
    SECRETION_SYSTEM_DETECTION(
        INPUT_PROCESSING.out.proteins,
        INPUT_PROCESSING.out.gene_info
    )

    // Phase 3: Secreted Protein Prediction
    SECRETED_PROTEIN_PREDICTION(
        INPUT_PROCESSING.out.proteins,
        SECRETION_SYSTEM_DETECTION.out.valid_systems
    )

    // Phase 4: Substrate Identification & Filtering
    SUBSTRATE_IDENTIFICATION(
        INPUT_PROCESSING.out.gene_order,
        SECRETION_SYSTEM_DETECTION.out.ss_components,
        SECRETED_PROTEIN_PREDICTION.out.predictions,
        SECRETION_SYSTEM_DETECTION.out.valid_systems
    )

    // Phase 5: Optional Annotation
    if (!params.skip_annotation) {
        OPTIONAL_ANNOTATION(
            SUBSTRATE_IDENTIFICATION.out.substrates,
            INPUT_PROCESSING.out.proteins
        )

        // Phase 6: Integration & Reporting
        INTEGRATION_REPORTING(
            SUBSTRATE_IDENTIFICATION.out.substrates,
            SUBSTRATE_IDENTIFICATION.out.substrates_unfiltered,
            OPTIONAL_ANNOTATION.out.annotations,
            SUBSTRATE_IDENTIFICATION.out.enrichment
        )
    } else {
        // Phase 6 without annotations
        INTEGRATION_REPORTING(
            SUBSTRATE_IDENTIFICATION.out.substrates,
            SUBSTRATE_IDENTIFICATION.out.substrates_unfiltered,
            Channel.empty(),
            SUBSTRATE_IDENTIFICATION.out.enrichment
        )
    }
}

/*
 * ──────────────────────────────────────────────────────────────────────
 *  COMPLETION HANDLER
 * ──────────────────────────────────────────────────────────────────────
 */

workflow.onComplete {
    log.info """\
        ─────────────────────────────────────────────
        ssign pipeline complete!
        Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration : ${workflow.duration}
        Output   : ${params.outdir}
        ─────────────────────────────────────────────
        """.stripIndent()
}
