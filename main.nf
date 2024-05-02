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

process COUNT_SIMULATOR {
    publishDir params.outdir, mode:'copy'
    
    input:
    path matrix
    path cell_types_annotation

    output:
    path "simulated_counts.csv"
    
    script:
    """
    Rscript $projectDir/bin/counts_simulator.R $matrix \
    $cell_types_annotation \
    simulated_counts.csv
    """
}

process TEMPLATE_MAKER {
    publishDir params.outdir, mode:'copy'
    
    input:
    path matrix
    path transcriptome
    path barcodes
    path gtf

    output:
    path "template.fa"
    
    script:
        def gtf = params.features != "transcript_id" ? "--features $params.features --gtf $gtf" : ""
        def unfiltered = barcodes.name != "no_barcode_counts" ? "--unfilteredBC $barcodes" : ""
    """
    python3.11 $projectDir/bin/template_maker.py  --matrix ${matrix} \
    --transcriptome $transcriptome \
    $unfiltered \
    $gtf \
    --thread 1 \
    --adapter $params.ADPTER_SEQ \
    --TSO $params.TSO_SEQ \
    --len_dT $params.dT_LENGTH \
    -o template.fa
    """
}

process ERRORS_SIMULATOR {
    publishDir params.outdir, mode:'copy'

    input:
    path template
    path error_model
    path qscore_model

    output:
    path "simulated.fastq"

    script:
        def error_model_arg = error_model.name != "no_error_model" ? "--badread-error-model $error_model" : ""
        def qscore_model_arg = qscore_model.name != "no_qscore_model" ? "--badread-qscore-model $qscore_model" : ""
    """
    python3.11 $projectDir/bin/error_simulator.py --thread $params.threads \
    -t $template \
    --badread-identity $params.badread_identity \
    $error_model_arg \
    $qscore_model_arg \
    -o simulated.fastq
    """
}

process GROUND_TRUTH {
    publishDir params.outdir, mode:'copy'
    
    input:
    path template

    output:
    path "ground_truth.tsv"

    script:
    """
    seqkit seq -n $template | sed 's/_/\t/g'  > ground_truth.tsv

    """
}

process QC {
    publishDir params.outdir, mode:'copy'
    
    input:
    path fastq

    output:
    path "${params.projetctName}.html"

    script:
    """
    python3.11 $projectDir/bin/QC.py -q $fastq -n $params.projetctName
    """
}

workflow {
    matrix_ch = Channel.fromPath(params.matrix, checkIfExists: true)

    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)

    barcodes_ch = params.filtered_out_bc != null ? Channel.fromPath(params.filtered_out_bc, checkIfExists: true) : 
                                             file("no_barcode_counts", type: "file")

    gtf_ch = params.features != "transcript_id" ? Channel.fromPath(params.gtf, checkIfExists: true) : 
                                             file("no_gtf", type: "file")

    cell_types_ch = params.sim_celltypes == true ? Channel.fromPath(params.cell_types_annotation, checkIfExists: true) : 
                                             file("no_cell_types", type: "file")

    error_model_ch = params.error_model != null ? Channel.fromPath(params.error_model, checkIfExists: true) : 
                                             file("no_error_model", type: "file")

    qscore_model_ch = params.qscore_model != null ? Channel.fromPath(params.qscore_model, checkIfExists: true) : 
                                             file("no_qscore_model", type: "file")

    if (params.sim_celltypes == true) {  
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