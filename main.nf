log.info """\
    SINGLE CELL NANOPORE READS SIMULATOR - N F   P I P E L I N E
    ===================================
    matrix                        : ${params.matrix} 
    transcriptome                 : ${params.transcriptome}
    matrix rownames               : ${params.features}
    GTF                           : ${params.gtf}
    filtered out barcodes         : ${params.filtered_out_bc}
    simulate cell types           : ${params.sim_celltypes}
    ref cell type annotation      : ${params.cell_types_annotation}
    error model                   : ${params.error_model}
    Qscore model                  : ${params.qscore_model}
    outdir                        : ${params.outdir}
    amplification rate            : ${params.amp}
    """
    .stripIndent()

include { COUNT_SIMULATOR } from './modules/modules.nf'
include { TEMPLATE_MAKER } from './modules/modules.nf'
include { ERRORS_SIMULATOR } from './modules/modules.nf'
include { GROUND_TRUTH } from './modules/modules.nf'
include { QC } from './modules/modules.nf'

workflow {
    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)
    
    matrix_ch = params.matrix != null ? file(params.matrix, checkIfExists: true) : 
                                             file("no_matrix_counts", type: "file")

    barcodes_ch = params.filtered_out_bc != null ? file(params.filtered_out_bc, checkIfExists: true) : 
                                            file("no_barcode_counts", type: "file")

    if (barcodes_ch.name == "no_barcode_counts" && matrix_ch.name == "no_matrix_counts") {
        log.error("Please provide one of the parameters: 'matrix counts' or 'barcodes counts'.")
        System.exit(1) 
    }
    
    gtf_ch = params.features != "transcript_id" ? Channel.fromPath(params.gtf, checkIfExists: true) : 
                                             file("no_gtf", type: "file")

    cell_types_ch = params.sim_celltypes == true ? Channel.fromPath(params.cell_types_annotation, checkIfExists: true) : 
                                             file("no_cell_types", type: "file")

    error_model_ch = params.error_model != null ? Channel.fromPath(params.error_model, checkIfExists: true) : 
                                             file("no_error_model", type: "file")

    qscore_model_ch = params.qscore_model != null ? Channel.fromPath(params.qscore_model, checkIfExists: true) : 
                                             file("no_qscore_model", type: "file")

    if (params.sim_celltypes) {  
        counts_ch = COUNT_SIMULATOR(matrix_ch, cell_types_ch)
        template_ch = TEMPLATE_MAKER(counts_ch, transcriptome_ch, barcodes_ch, gtf_ch)

    } else { 
        template_ch = TEMPLATE_MAKER(matrix_ch, transcriptome_ch, barcodes_ch, gtf_ch)
    }
    gr_truth_ch = GROUND_TRUTH(template_ch)
    error_ch    = ERRORS_SIMULATOR(template_ch, error_model_ch, qscore_model_ch)
    qc_ch       = QC(error_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/${params.projetctName}.html\n" : "Oops .. something went wrong" )
}