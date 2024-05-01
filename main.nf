log.info """\
    SINGLE CELL NANOPORE READS SIMULATOR - N F   P I P E L I N E
    ===================================
    matrix                   : ${params.matrix} 
    transcriptome            : ${params.transcriptome}
    barcodes                 : ${params.barcodes}
    read length distribution : ${params.ref_distribution}
    GTF                      : ${params.gtf} 
    amplification rate       : ${params.amp}
    outdir                   : ${params.outdir}
    """
    .stripIndent()

process READ_GENERATOR {
    publishDir params.outdir, mode:'copy'
    
    input:
    path transcriptome
    path barcodes
    path ref_dist

    output:
    path "template.fa"
    
    script:
        def dist_arg = ref_dist.name != "random_dist" ? "-f $ref_dist" : ""
    """
    python3.11 $projectDir/bin/perfect_read_generator.py  -i ${barcodes} \
    -r $transcriptome \
    $dist_arg \
    --amp-rate $params.amp \
    --thread $params.threads \
    -o template.fa
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
        def gtf = params.features != "transcript_id" ? "--gtf $gtf" : ""
    """
    python3.11 $projectDir/bin/template_maker.py  --matrix ${matrix} \
    --unfilteredBC ${barcodes} \
    --transcriptome $transcriptome \
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
    seqkit seq -n $template | sed 's/_/\t/3' > ground_truth.tsv

    """
}

//| sed 's/_/\t/g' > ground_truth.tsv
 //   sed -i -e '1i\BC\tUMI\tindex\tfeature' ground_truth.tsv

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
    barcodes_ch = Channel.fromPath(params.barcodes, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf)

    dist_ch = params.ref_distribution != null ? Channel.fromPath(params.ref_distribution, checkIfExists: true) : 
                                             file("random_dist", type: "file")
    error_model_ch = params.error_model != null ? Channel.fromPath(params.error_model, checkIfExists: true) : 
                                             file("no_error_model", type: "file")
    qscore_model_ch = params.qscore_model != null ? Channel.fromPath(params.qscore_model, checkIfExists: true) : 
                                             file("no_qscore_model", type: "file")
    
    template_ch = TEMPLATE_MAKER(matrix_ch, transcriptome_ch, barcodes_ch, dist_ch)
    quant_ch    = ERRORS_SIMULATOR(template_ch, error_model_ch, qscore_model_ch)
    gr_truth_ch = GROUND_TRUTH(template_ch)
    qc_ch       = QC(quant_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/${params.projetctName}.html\n" : "Oops .. something went wrong" )
}