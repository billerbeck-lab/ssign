/*
 * Phase 3: Secreted Protein Prediction
 *
 * Runs DeepLocPro, DeepSecE, and optionally SignalP to predict
 * which proteins are secreted and by which pathway.
 * Then cross-validates predictions across tools.
 */

include { DEEPLOCPRO }       from '../modules/local/deeplocpro'
include { DEEPSECE }         from '../modules/local/deepsece'
include { SIGNALP }          from '../modules/local/signalp'
include { CROSS_VALIDATE }   from '../modules/local/cross_validate'

workflow SECRETED_PROTEIN_PREDICTION {
    take:
    ch_proteins       // [sample_id, proteins.faa]
    ch_valid_systems  // [sample_id, valid_systems.tsv]

    main:
    // Run prediction tools in parallel
    DEEPLOCPRO(ch_proteins)
    DEEPSECE(ch_proteins)

    // SignalP is optional (DTU license required)
    if (!params.skip_signalp && params.signalp_path) {
        SIGNALP(ch_proteins)
        ch_signalp = SIGNALP.out.predictions
    } else {
        ch_signalp = Channel.empty()
    }

    // Cross-validate: merge predictions, flag unreliable DSE T3SS
    ch_for_validation = DEEPLOCPRO.out.predictions
        .join(DEEPSECE.out.predictions)
        .join(ch_valid_systems)

    CROSS_VALIDATE(ch_for_validation)

    emit:
    predictions = CROSS_VALIDATE.out.validated  // [sample_id, predictions_validated.tsv]
    deeplocpro  = DEEPLOCPRO.out.predictions    // [sample_id, deeplocpro.tsv]
    deepsece    = DEEPSECE.out.predictions       // [sample_id, deepsece.tsv]
}
