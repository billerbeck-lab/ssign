/*
 * GENERATE_REPORT — create HTML and text summary reports
 */

process GENERATE_REPORT {
    label 'process_low'

    input:
    path(master_csvs)
    path(enrichment)

    output:
    path("ssign_report.html"), emit: report
    path("ssign_report.txt"),  emit: report_txt

    script:
    """
    generate_report.py \\
        --master-csvs ${master_csvs} \\
        --enrichment ${enrichment} \\
        --out-html ssign_report.html \\
        --out-txt ssign_report.txt
    """
}
