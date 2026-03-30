/*
 * PROTPARAM — BioPython physicochemical property computation
 *
 * Computes: MW, pI, GRAVY, instability index, aromaticity, charge at pH7
 */

process PROTPARAM {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('protparam'), path("${sample_id}_protparam.csv"), emit: annotations

    script:
    """
    compute_protparam.py \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --output ${sample_id}_protparam.csv
    """
}
