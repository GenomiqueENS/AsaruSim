process PCR_SIMULATOR {
    publishDir params.outdir, mode:'copy'
    
    input:
    path fasta

    output:
    path "amplified_reads.fa" 
    
    script:
        def totalNamber = params.totalNamber == false? "" : "--totalNamber $params.totalNamber"

    """
    python3.11 $projectDir/bin/PCR.py -f ${fasta} \
    --cycles $params.cycles \
    --dup $params.dup_rate \
    --error $params.error_rate \
    --thread $params.threads \
    $totalNamber \
    --seed 2024 \
    --out amplified_reads.fa
    """
}