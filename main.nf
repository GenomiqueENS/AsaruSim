params.help = false

help_message = """
    Usage:
        nextflow run main.nf --matrix <path> [options]

    Description:
        ASARU - SINGLE CELL NANOPORE READS SIMULATOR PIPELINE
        Simulates nanopore reads from single-cell data.

    Required parameters:
        --matrix <path>                Path to the expression matrix.
        --transcriptome <path>         Path to the transcriptome file.
        --features <value>              Features type provided in the count matrix (The matrix row names)
                                       (transcript_id, gene_id or gene_name)[default: transcript_id]

    Optional parameters:
        --barcodes_counts <value>      Distribution of barcode counts.
        --full_length <boolean>        Indicates if transcripts are full length [default: false].
        --simulate_cell_types <boolean>Simulate cell types [default: false].
        --cell_type_annotation <path>  Path to cell type annotation.
        --gtf <path>                   Path to the GTF file.
        --trained_model <path>         Pre-trained error model.
        --identity_model <value>       Identity model for Badread.
        --error_model <path>           Custom error model.
        --qscore_model <path>          Q score model.
        --build_model <boolean>        Build an error model [default: false].
        --fastq_model <path>           Path to the FASTQ model.
        --ref_genome <path>            Path to the reference genome.
        --umi_duplication <value>      UMI duplication.
        --pcr_cycles <int>             Number of PCR amplification cycles.
        --pcr_dup_rate <value>         PCR duplication rate.
        --pcr_error_rate <value>       PCR error rate.
        --pcr_total_reads <int>        Total number of PCR reads.
        --truncation_model <path>      Path to truncation probabilities .CSV file [default: bin/models/truncation_default_model.csv].
        --outdir <path>                Output directory.

    Example:
        nextflow run main.nf --matrix 'path/to/matrix.csv' --transcriptome 'path/to/transcriptome.fa' --outdir 'results'

    """

if (params.help) {
    println help_message
    System.exit(0)
}

log.info """\
    ASARU - SINGLE CELL NANOPORE READS SIMULATOR P I P E L I N E
    ===================================
    matrix                        : ${params.matrix} 
    transcriptome                 : ${params.transcriptome}
    matrix rownames               : ${params.features}
    barcodes counts distribution  : ${params.bc_counts}
    full length transcripts       : ${params.full_length}
    novel transcripts             : ${params.novel_transcripts}
    simulate cell types           : ${params.sim_celltypes}
    cell type annotation          : ${params.cell_types_annotation}
    GTF                           : ${params.gtf}
    pre-trained error model       : ${params.trained_model}
    identity model                : ${params.badread_identity}
    error model                   : ${params.error_model}
    Qscore model                  : ${params.qscore_model}
    build new mode                : ${params.build_model}
    FASTQ model                   : ${params.fastq_model}
    Truncation model              : ${params.truncation_model}
    reference genome              : ${params.ref_genome}
    UMI duplication               : ${params.umi_duplication}
    PCR amplification cycles      : ${params.pcr_cycles}
    PCR duplication rate          : ${params.pcr_dup_rate}
    PCR error rate                : ${params.pcr_error_rate}
    Poly dT mean length           : ${params.dT_LENGTH}
    Adapter sequence              : ${params.ADPTER_SEQ}
    TSO sequence                  : ${params.TSO_SEQ}
    outdir                        : ${params.outdir}
    """
    .stripIndent()

logo_ch = channel.fromPath('./images/asarusim_v2.png')
config_params_ch = Channel.value(params)
                          .map { p -> p.collect { key, value -> "--conf_params $key:$value" }}

workflow_params_ch = Channel.from(workflow).map { p -> p.collect { value -> "--work_params '$value'" }}

include { SUBSAMPLE } from './modules/errorModel.nf'
include { ALIGNMENT } from './modules/errorModel.nf'
include { ERROR_MODLING } from './modules/errorModel.nf'
include { IDENTITY_ESTIMATION } from './modules/errorModel.nf'
include { LENGTH_ESTIMATION } from './modules/errorModel.nf'
include { BATCH_SIZE_FOR_BADREAD } from './modules/errorModel.nf'
include { TRUNCATION_ESTIMATION } from './modules/truncationModel.nf'
include { GFF_TO_FASTA } from './modules/gff_to_fasta.nf'

include { COUNT_SIMULATOR } from './modules/modules.nf'
include { TEMPLATE_MAKER } from './modules/modules.nf'
include { PCR_SIMULATOR } from './modules/PCR.nf'
include { ERRORS_SIMULATOR } from './modules/modules.nf'
include { GROUND_TRUTH } from './modules/modules.nf'
include { QC } from './modules/modules.nf'

def validSPARSIMOptions = ['Bacher', 'Camp', 'Chu', 'Engel', 'Horning', 'Tung', 'Zheng', 'Macosko', 'Saunders']
validSPARSIMOptions = validSPARSIMOptions.collect { it + '_param_preset' }

workflow {
    //transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)

    
    if (params.matrix in validSPARSIMOptions){
        matrix_ch = file(params.matrix, type: "file")

    } else {
        matrix_ch =  params.matrix != null ? file(params.matrix, checkIfExists: true) : 
                                             file("no_matrix_counts", type: "file")
    }

    barcodes_ch = params.bc_counts != null ? file(params.bc_counts, checkIfExists: true) : 
                                            file("no_barcode_counts", type: "file")

    if (barcodes_ch.name == "no_barcode_counts" && matrix_ch.name == "no_matrix_counts" ) {
        log.error("Please provide one of the parameters: 'matrix counts' or 'barcodes counts'.")
        System.exit(1) 
    }

    transcriptome_ch =  params.transcriptome != null ? file(params.transcriptome, checkIfExists: true) : 
                                             file("no_transcriptome", type: "file")
    
    gtf_ch = (params.features != "transcript_id" || params.transcriptome == null) ? Channel.fromPath(params.gtf, checkIfExists: true) : 
                                             file("no_gtf", type: "file")

    cell_types_ch = params.sim_celltypes == true && (params.matrix !in validSPARSIMOptions)? Channel.fromPath(params.cell_types_annotation, checkIfExists: true) : 
                                             file("no_cell_types", type: "file")

    error_model_ch = params.error_model != null ? file(params.error_model) : 
                                             file("no_error_model", type: "file")

    qscore_model_ch = params.qscore_model != null ? file(params.qscore_model) : 
                                             file("no_qscore_model", type: "file")
    
    identity_ch = params.badread_identity != null ? channel.from(params.badread_identity) :
                                                    channel.from("96,2,98")
    length_dist_ch = params.length_dist != null ? channel.from(params.length_dist) :
                                                    channel.from("0.37,0.0,824.94")
    truncation_ch = params.truncation_model != null ? file(params.truncation_model) :
                                                    file("bin/models/truncation_default_model.csv", type: "file")
    
    if (params.transcriptome == null) {
        if (params.gtf == null){
            error("Please provide an annotation file in GTF format! e.g '--gtf annoation.gtf' or transcriptome in FASTA format") 
        }else{
            genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
            transcriptome_ch = GFF_TO_FASTA(gtf_ch, genome_ch)
        }
    }

    if (params.build_model) {
        fastq_ch = Channel.fromPath(params.fastq_model, checkIfExists: true)
        sub_fastq_ch = SUBSAMPLE(fastq_ch)
        paf_ch = ALIGNMENT(sub_fastq_ch, transcriptome_ch)
        ERROR_MODLING(sub_fastq_ch, transcriptome_ch, paf_ch)
        identity_ch = IDENTITY_ESTIMATION(paf_ch)
        truncation_ch = TRUNCATION_ESTIMATION(paf_ch)
        if (params.features != "transcript_id") {
            length_dist_ch = LENGTH_ESTIMATION(fastq_ch)
        }
        error_model_ch = ERROR_MODLING.out.error_model_ch
        qscore_model_ch = ERROR_MODLING.out.qscore_model_ch
    }

    if (params.sim_celltypes) {  
        counts_ch = COUNT_SIMULATOR(matrix_ch, cell_types_ch)
        TEMPLATE_MAKER(counts_ch, transcriptome_ch, barcodes_ch, gtf_ch, length_dist_ch, truncation_ch)

    } else { 
        TEMPLATE_MAKER(matrix_ch, transcriptome_ch, barcodes_ch, gtf_ch, length_dist_ch, truncation_ch)
    }

    template_fa_ch = TEMPLATE_MAKER.out.template
    template_log_ch = TEMPLATE_MAKER.out.logfile

    if (params.pcr_cycles > 0) {
        template_fa_ch = PCR_SIMULATOR(template_fa_ch, template_log_ch)
    }

    gr_truth_ch = GROUND_TRUTH(template_fa_ch)
    batch_size_ch  = BATCH_SIZE_FOR_BADREAD(template_fa_ch)
    error_ch    = ERRORS_SIMULATOR(template_fa_ch, error_model_ch, qscore_model_ch, identity_ch) //, batch_size_ch.toInteger())
    qc_ch       = QC(error_ch, config_params_ch, workflow_params_ch, logo_ch, template_log_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/${params.projetctName}.html\n" : "Oops .. something went wrong" )
}