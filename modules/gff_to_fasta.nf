process GFF_TO_FASTA {
    input:
    path gff
    path fasta
    
    output:
    path "transcriptome_assembly.fa",  emit: transcriptome_assembly

    script:
    """
    gffread $gff -g $fasta -w transcriptome_assembly.fa
    """
}