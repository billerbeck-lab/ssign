/*
 * GENERATE_FIGURES — publication-quality visualization
 */

process GENERATE_FIGURES {
    label 'process_medium'

    input:
    path(master_csvs)

    output:
    path("figures/"), emit: figures

    script:
    """
    mkdir -p figures
    generate_figures.py \\
        --master-csvs ${master_csvs} \\
        --outdir figures \\
        --dpi ${params.dpi}
    """
}
