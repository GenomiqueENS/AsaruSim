#!/bin/bash

# Check if main.nf exists
if [ ! -f "main.nf" ]; then
    echo "ERROR: cannot found main.nf." >&2
    exit 1
fi

# Check if wget is installed
if [ ! -x "$(command -v wget)" ]; then
    echo "ERROR: wget is not installed/not in the PATH." >&2
    exit 1
fi

# Download dataset
if [ ! -f "sub_pbmc_datasets.tar.gz" ]; then
    wget https://zenodo.org/records/12731409/files/sub_pbmc_datasets.tar.gz
fi

# Untar dataset
if [ ! -d "dataset" ]; then
    tar -xzf sub_pbmc_datasets.tar.gz
fi

# Download reference transcriptome
if [ ! -f "dataset/Homo_sapiens.GRCh38.cdna.all.fa" ]; then
    wget -O - https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz | gzip -d > dataset/Homo_sapiens.GRCh38.cdna.all.fa
fi

# Download GTF
if [ ! -f "dataset/GRCh38-2020-A-genes.gtf" ]; then
    wget -O - https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz | gzip -d > dataset/GRCh38-2020-A-genes.gtf
fi

# Download reference genome
if [ ! -f "dataset/GRCh38-2020-A-genome.fa" ]; then
    wget -O - https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz | gzip -d > dataset/GRCh38-2020-A-genome.fa
fi

# Check if nextflow is in the path
if [ ! -x "$(command -v nextflow)" ]; then
    echo "ERROR: nextflow is not installed/not in the PATH." >&2
    exit 1
fi

echo "* Execute basic workflow"
nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                     --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                     --features gene_name \
                     --gtf dataset/GRCh38-2020-A-genes.gtf