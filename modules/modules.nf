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
        def unfiltered = barcodes.name != "no_barcode_counts"? "--unfilteredBC $barcodes" : ""
    """
    python3.11 $projectDir/bin/template_maker.py  --matrix ${matrix} \
    --transcriptome $transcriptome \
    $unfiltered \
    $gtf \
    --thread 1 \
    --amp $params.amp \
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
    def error_model_arg = ""
    def qscore_model_arg = ""
    if (params.trained_model) {
        error_model_arg = "--badread-error-model ${params.trained_model}"
        qscore_model_arg = "--badread-qscore-model ${params.trained_model}"
    } else {
        error_model_arg = error_model.name != "no_error_model" ? "--badread-error-model $error_model" : ""
        qscore_model_arg = qscore_model.name != "no_qscore_model" ? "--badread-qscore-model $qscore_model" : ""
    }
    """
    python3.11 $projectDir/bin/SLSim/error_simulator.py --thread $params.threads \
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