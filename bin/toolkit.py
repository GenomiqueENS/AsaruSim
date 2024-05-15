# Author: Ali Hamraoui
# Date: 25-04-2024

import random
from scipy import stats

def fetch_transcript_id_by_gene(gene, transcripts_index):
    transcript_ids = transcripts_index.get(gene, [])
    if transcript_ids:
        return transcript_ids[0]
    else:
        return None


def cut_sequence(sequence, length_distribution):
    """
    Function to truncate sequences to a random length based on a given log-normal distribution.
    Inputs:
    sequences: List of sequences (byte strings).
    length_distribution: Tuple (shape, loc, scale) defining the parameters of the log-normal distribution.
    Returns: List of truncated sequences.
    """
    shape, loc, scale = length_distribution
    cut_length = int(stats.lognorm.rvs(shape, loc, scale))-100
    cut_length = min(cut_length, len(sequence))
    return sequence[:cut_length]


def parse_gtf(gtf_file, index_by, protein_coding=False):
    
    transcripts = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            columns = line.strip().split("\t")
            if len(columns) > 8:
                Attribute = columns[8].replace('"', '')
                feature = columns[2]
                if feature != 'transcript': 
                    continue
                    
                attributes = dict(item.strip().split(' ',1) for item in Attribute.strip().split(';') if item.strip())
                transcript_id = attributes['transcript_id']
                gene_id = attributes['gene_id']
                gene_name = attributes['gene_name']
                transcript_type = attributes.get('transcript_biotype', attributes.get('transcript_type', ''))
                
                index_key = None
                if index_by == 'gene_id':
                    index_key = gene_id
                elif index_by == 'gene_name':
                    index_key = gene_name
                
                if index_key:
                    if index_key not in transcripts:
                        transcripts[index_key] = []
                    
                    if not protein_coding:
                        transcripts[index_key].append(transcript_id)
                    elif protein_coding and transcript_type == 'protein_coding':
                        transcripts[index_key].append(transcript_id)

    return transcripts


def custom_key_function(header):
    """Extrait l'identifiant sans la version."""
    return header.split('.')[0]


def generate_random_umi(length=12):
    """Génère une séquence UMI aléatoire de longueur spécifiée."""
    bases = 'ATCG'
    return ''.join(random.choice(bases) for _ in range(length))


def cDNA_pool(trns, transcriptome):
    if trns in transcriptome:
        sequence = transcriptome[trns]
        return sequence[:]