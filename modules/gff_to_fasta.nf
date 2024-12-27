process GFF_TO_FASTA {
    input:
    path gff
    path ref_genome
    
    output:
    path "transcriptome_assembly.fa",  emit: transcriptome_assembly

    script:
    def genome = ref_genome.name != "no_genome"? "--ref_genome $ref_genome" : ""

    """
    gffread $gff -g $genome -w transcriptome_assembly.fa
    """
}