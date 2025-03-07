process SUBSAMPLE {
    input:
    path reads_fastq
    
    output:
    path "sub_reads.fq.gz"

    
    script:
    """
    seqkit head -n 100000 $reads_fastq -o sub_reads.fq.gz
    """
}

process ALIGNMENT {
    input:
    path reads_fastq
    path ref_transcriptome
    

    output:
    path "alignments.paf.gz"

    
    script:
    """
    /minimap2-2.28_x64-linux/minimap2 -t $params.threads -c -x map-ont $ref_transcriptome $reads_fastq | gzip > alignments.paf.gz
    """
}

process ERROR_MODLING {
    input:
    path reads_fastq
    path ref_transcriptome
    path alignment
    
    output:
    path "new_error_model",  emit: error_model_ch
    path "new_qscore_model",  emit: qscore_model_ch

    script:
    """
    python3.11 /bin/Badread/badread-runner.py error_model --reference $ref_transcriptome --reads $reads_fastq --alignment $alignment > new_error_model
    python3.11 /bin/Badread/badread-runner.py qscore_model --reference $ref_transcriptome --reads $reads_fastq --alignment $alignment > new_qscore_model
    """
}


process IDENTITY_ESTIMATION {
    input:
    path alignment
    
    output:
    stdout emit: beta_params

    script:
    def identity_model = params.identity_model == null ? "gap_excluded_identity" : params.identity_model
    """
    python3.11 $projectDir/bin/params_estimator/estimator.py -r identity -m $identity_model -p $alignment
    """
}

process LENGTH_ESTIMATION {
    input:
    path reads_fastq
    
    output:
    stdout emit: length_params

    script:
    """
    python3.11 $projectDir/bin/params_estimator/estimator.py -r length -f $reads_fastq -t $params.threads 
    """
}

process BATCH_SIZE_FOR_BADREAD {
    input:
    path template
    
    output:
    stdout emit: batch_size

    script:
    """
    wc -l $template | cut -d " " -f1
    """
}