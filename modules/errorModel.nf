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
    path ref_genome
    

    output:
    path "alignments.paf.gz"

    
    script:
    """
    /minimap2-2.28_x64-linux/minimap2 -t $params.threads -c -x map-ont $ref_genome $reads_fastq | gzip > alignments.paf.gz
    """
}

process ERROR_MODLING {
    input:
    path reads_fastq
    path ref_genome
    path alignment
    
    output:
    path "new_error_model",  emit: error_model_ch
    path "new_qscore_model",  emit: qscore_model_ch

    script:
    """
    python3.11 /bin/Badread/badread-runner.py error_model --reference $ref_genome --reads $reads_fastq --alignment $alignment > new_error_model
    python3.11 /bin/Badread/badread-runner.py qscore_model --reference $ref_genome --reads $reads_fastq --alignment $alignment > new_qscore_model
    """
}

process IDENTITY_ESTIMATION {
    input:
    path alignment
    
    output:
    stdout emit: beta_params

    script:
    """
    python3.11 $projectDir/bin/params_estimator/estimator.py -r identity -p $alignment 
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
