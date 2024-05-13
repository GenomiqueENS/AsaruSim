log.info """\
    ASARU - SINGLE CELL NANOPORE READS SIMULATOR P I P E L I N E
    ===================================
    matrix                        : ${params.matrix} 
    transcriptome                 : ${params.transcriptome}
    matrix rownames               : ${params.features}
    barcodes counts distribution  : ${params.bc_counts}
    simulate cell types           : ${params.sim_celltypes}
    cell type annotation          : ${params.cell_types_annotation}
    GTF                           : ${params.gtf}
    pre-trained error model       : ${params.trained_model}
    identity model                : ${params.badread_identity}
    error model                   : ${params.error_model}
    Qscore model                  : ${params.qscore_model}
    build erro model              : ${params.build_model}
    FASTQ model                   : ${params.model_fastq}
    reference genome              : ${params.ref_genome}
    amplification rate            : ${params.amp}
    outdir                        : ${params.outdir}
    """
    .stripIndent()

include { SUBSAMPLE } from './modules/errorModel.nf'
include { ALIGNMENT } from './modules/errorModel.nf'
include { ERROR_MODLING } from './modules/errorModel.nf'
include { IDENTITY } from './modules/errorModel.nf'

include { COUNT_SIMULATOR } from './modules/modules.nf'
include { TEMPLATE_MAKER } from './modules/modules.nf'
include { ERRORS_SIMULATOR } from './modules/modules.nf'
include { GROUND_TRUTH } from './modules/modules.nf'
include { QC } from './modules/modules.nf'

workflow {
    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)
    
    matrix_ch = params.matrix != null ? file(params.matrix, checkIfExists: true) : 
                                             file("no_matrix_counts", type: "file")

    barcodes_ch = params.bc_counts != null ? file(params.bc_counts, checkIfExists: true) : 
                                            file("no_barcode_counts", type: "file")

    if (barcodes_ch.name == "no_barcode_counts" && matrix_ch.name == "no_matrix_counts") {
        log.error("Please provide one of the parameters: 'matrix counts' or 'barcodes counts'.")
        System.exit(1) 
    }
    
    gtf_ch = params.features != "transcript_id" ? Channel.fromPath(params.gtf, checkIfExists: true) : 
                                             file("no_gtf", type: "file")

    cell_types_ch = params.sim_celltypes == true ? Channel.fromPath(params.cell_types_annotation, checkIfExists: true) : 
                                             file("no_cell_types", type: "file")

    error_model_ch = params.error_model != null ? file(params.error_model) : 
                                             file("no_error_model", type: "file")

    qscore_model_ch = params.qscore_model != null ? file(params.qscore_model) : 
                                             file("no_qscore_model", type: "file")
    
    if (params.build_model) {
        fastq_ch = Channel.fromPath(params.model_fastq, checkIfExists: true)
        genome_ch = Channel.fromPath(params.ref_genome, checkIfExists: true)
        sub_fastq_ch = SUBSAMPLE(fastq_ch)
        paf_ch = ALIGNMENT(sub_fastq_ch, genome_ch)
        ERROR_MODLING(sub_fastq_ch, genome_ch, paf_ch)
        error_model_ch = ERROR_MODLING.out.error_model_ch
        qscore_model_ch = ERROR_MODLING.out.qscore_model_ch
    }

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