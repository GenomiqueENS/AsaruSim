process PCR_SIMULATOR {
    publishDir params.outdir, mode:'copy'
    
    input:
    path fasta

    output:
    path "amplified_reads.fa" 
    
    script:
        def totalNamber = params.pcr_total_reads == null? "" : "--totalNamber $params.pcr_total_reads"

    """
    python3.11 $projectDir/bin/PCR.py -f ${fasta} \
    --cycles $params.pcr_cycles \
    --dup $params.pcr_dup_rate \
    --error $params.pcr_error_rate \
    --thread $params.threads \
    $totalNamber \
    --seed 2024 \
    --out amplified_reads.fa
    """
}