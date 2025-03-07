process COUNT_SIMULATOR {
    publishDir params.outdir, mode:'copy'
    cache false
    
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
    publishDir params.outdir, mode:'copy', pattern: "log.csv"
    
    input:
    path matrix
    path transcriptome
    path barcodes
    path ref_gtf
    val length_dist
    path truncation
    path intron_r
    path ref_genome

    output:
    path "template.fa", emit: template
    path "log.csv", emit: logfile
    
    script:
        def matrix = matrix.name != "no_matrix_counts"? "--matrix $matrix" : ""
        def gtf = params.features != "transcript_id" ? "--features $params.features --gtf $ref_gtf" : ""
        def unfiltered = barcodes.name != "no_barcode_counts"? "--unfilteredBC $barcodes" : ""
        def full_length = params.full_length? "--full_length" : "--truncation_model $truncation"
        def intron_retention = params.intron_retention? "--intron_retention --IR_model $intron_r --gtf $ref_gtf" : ""
        def genome = ref_genome.name != "no_genome"? "--ref_genome $ref_genome" : ""
    """
    python3.11 $projectDir/bin/AsaruSim.py template_maker  $matrix \
    --transcriptome $transcriptome \
    $unfiltered \
    $gtf \
    --thread 1 \
    --amp $params.umi_duplication \
    --adapter $params.ADPTER_SEQ \
    --TSO $params.TSO_SEQ \
    --len_dT $params.dT_LENGTH \
    $full_length \
    $intron_retention \
    $genome \
    --unspliced_ratio $params.unspliced_ratio \
    --length_dist ${length_dist.trim()} \
    -o template.fa \
    --log log.csv
    """
}

process ERRORS_SIMULATOR {
    publishDir params.outdir, mode:'copy', pattern: '*.fa.gz'

    input:
    path template
    path error_model
    path qscore_model
    val identity
    val batch_size

    output:
    path "simulated.fastq", emit: fastq
    path "template.fa.gz", optional: true
    path "amplified_reads.fa.gz", optional: true

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
    python3.11 $projectDir/bin/AsaruSim.py call_badread  -t $template \
        --badread-identity ${identity.trim()} \
        $error_model_arg \
        $qscore_model_arg \
        --batch-size ${Math.min((batch_size/(2*params.threads)).toInteger(), 100000)} \
        --thread $params.threads \
        -o simulated.fastq &&
        gzip -f $template
    """
}

process GROUND_TRUTH {
    publishDir params.outdir, mode:'copy'
    
    input:
    path template

    output:
    path "ground_truth.tsv.gz"

    script:
    """
    seqkit seq -n $template | sed 's/_/\t/g' | gzip > ground_truth.tsv.gz

    """
}

process QC {
    publishDir params.outdir, mode:'copy'
    
    input:
    path fastq
    val conf_params
    val work_params
    path logo
    path log_csv

    output:
    path "${params.projetctName}.html"
    path "simulated.fastq.gz"

    script:
    """
    bedread_version=\$(echo "-e \$(badread --version | head -n 1)")
    nextflow_version=\$(echo "-e \$(Rscript -e "ip <- installed.packages()[, c('Package', 'Version')]; 
                                    cat(paste(ip[,1], ip[,2], sep=': '), sep='\n')" | grep SPARSim)")

    python3.11 $projectDir/bin/QC.py -q $fastq --stats $log_csv -n $params.projetctName \
            ${conf_params.join(' ').replaceAll(/[\(\)\$]/, "")} \
            ${work_params.join('').replaceAll(/[\(\)\$]/, "")} \
            "\$bedread_version" \
            "\$nextflow_version" \
    && gzip -f $fastq
    """
}