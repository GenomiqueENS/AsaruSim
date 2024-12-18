process PCR_SIMULATOR {
    publishDir params.outdir, mode:'copy'
    
    input:
    path fasta
    path template_log

    output:
    path "amplified_reads.fa" 
    
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
    --out amplified_reads.fa
    """
}