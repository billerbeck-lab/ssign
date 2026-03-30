/*
 * Phase 5: Optional Annotation
 *
 * All tools are independently skippable.
 * Each database tool supports local and remote modes.
 */

include { BLASTP_LOCAL }     from '../modules/local/blastp'
include { BLASTP_REMOTE }    from '../modules/local/blastp'
include { HHSUITE_LOCAL }    from '../modules/local/hhsuite'
include { HHSUITE_REMOTE }   from '../modules/local/hhsuite'
include { FOLDSEEK_LOCAL }   from '../modules/local/foldseek'
include { INTERPROSCAN_LOCAL }  from '../modules/local/interproscan'
include { INTERPROSCAN_REMOTE } from '../modules/local/interproscan'
include { PROTPARAM }        from '../modules/local/protparam'

workflow OPTIONAL_ANNOTATION {
    take:
    ch_substrates  // [sample_id, substrates_filtered.tsv]
    ch_proteins    // [sample_id, proteins.faa]

    main:
    // Collect all annotation outputs into a single channel per sample
    ch_annotations = Channel.empty()

    // --- BLASTp ---
    if (!params.skip_blastp) {
        if (params.blastp_mode == 'local') {
            BLASTP_LOCAL(ch_substrates.join(ch_proteins))
            ch_annotations = ch_annotations.mix(BLASTP_LOCAL.out.annotations)
        } else {
            BLASTP_REMOTE(ch_substrates.join(ch_proteins))
            ch_annotations = ch_annotations.mix(BLASTP_REMOTE.out.annotations)
        }
    }

    // --- HH-suite ---
    if (!params.skip_hhsuite) {
        if (params.hhsuite_mode == 'local') {
            HHSUITE_LOCAL(ch_substrates.join(ch_proteins))
            ch_annotations = ch_annotations.mix(HHSUITE_LOCAL.out.annotations)
        } else {
            HHSUITE_REMOTE(ch_substrates.join(ch_proteins))
            ch_annotations = ch_annotations.mix(HHSUITE_REMOTE.out.annotations)
        }
    }

    // --- Foldseek ---
    if (!params.skip_foldseek) {
        FOLDSEEK_LOCAL(ch_substrates.join(ch_proteins))
        ch_annotations = ch_annotations.mix(FOLDSEEK_LOCAL.out.annotations)
    }

    // --- InterProScan ---
    if (!params.skip_interproscan) {
        if (params.interproscan_mode == 'local') {
            INTERPROSCAN_LOCAL(ch_substrates.join(ch_proteins))
            ch_annotations = ch_annotations.mix(INTERPROSCAN_LOCAL.out.annotations)
        } else {
            INTERPROSCAN_REMOTE(ch_substrates.join(ch_proteins))
            ch_annotations = ch_annotations.mix(INTERPROSCAN_REMOTE.out.annotations)
        }
    }

    // --- ProtParam (always local, lightweight) ---
    if (!params.skip_protparam) {
        PROTPARAM(ch_substrates.join(ch_proteins))
        ch_annotations = ch_annotations.mix(PROTPARAM.out.annotations)
    }

    emit:
    annotations = ch_annotations  // [sample_id, tool_name, annotations.csv]
}
