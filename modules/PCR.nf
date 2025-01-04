process PCR_SIMULATOR {
    publishDir params.outdir, mode:'copy', pattern: 'template.fa.gz'
    
    input:
    path fasta
    path template_log

    output:
    path "amplified_reads.fa", emit: fasta
    path "template.fa.gz", optional: true
    
    script:
        def totalNumber = params.pcr_total_reads == null? "" : "--totalNumber $params.pcr_total_reads"

    """
    python3.11 $projectDir/bin/AsaruSim.py PCR -f ${fasta} \
    --cycles $params.pcr_cycles \
    --dup $params.pcr_dup_rate \
    --error $params.pcr_error_rate \
    --thread $params.threads \
    $totalNumber \
    --seed $params.seed \
    --maker_log ${template_log} \
    --out amplified_reads.fa &&
    gzip -f ${fasta}
    """
}