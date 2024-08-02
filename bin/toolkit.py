# Author: Ali Hamraoui
# Date: 25-04-2024
import multiprocessing as mp
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import random
from scipy import stats
import numpy as np

def adjust_lambda_cut(base_lambda, seq_len):
    return base_lambda / seq_len


def cut_sequence_proba(sequence, base_lambda=1):
        lambda_param = adjust_lambda_cut(base_lambda, len(sequence))
        prob = selection_probas(range(len(sequence)), lambda_param)
        cut_length = np.random.choice(range(len(sequence)), p=prob[::-1])
        return sequence[-cut_length:]


def adjust_lambda(base_lambda, lengths):
    std_dev = np.std(lengths)
    return base_lambda / std_dev


def selection_probas(lengths, lambda_param):
    probabilities = np.exp(-lambda_param * np.array(lengths)**2)
    return probabilities / probabilities.sum()


def logN_probas(transcript_lengths, shape = 0.39146, scale = 830):
    probabilities = [stats.lognorm.pdf(length, shape, scale=scale) for length in transcript_lengths]
    return probabilities/sum(probabilities)


def fetch_transcript_id_by_gene(gene, transcripts_index, base_lambda=0.001):
    transcript_ids = transcripts_index.get(gene, [])
    
    if transcript_ids:
        trans_len_dict = {trns: length for trns, length in transcript_ids}
        
        lengths = list(trans_len_dict.values())
        transcripts = list(trans_len_dict.keys())

        probabilities = logN_probas(lengths)

        if np.any(np.isnan(probabilities)):
            return random.choice(transcripts)
        
        selected_index = np.random.choice(range(len(lengths)), p=probabilities)
        
        return transcripts[selected_index]
        
    else:
        return None


def fetch_transcript_id_by_gene_old(gene, transcripts_index):
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
    return sequence[-cut_length:]


def parse_gtf(gtf_file, index_by, protein_coding=False):

    if index_by == 'gene_id' or index_by == 'gene_name':
        pass
    else:
        logging.error("\033[91mError: Please specify the feature type in the provided matrix (transcript_id, gene_id or gene_name) using '--features' argument.\033[0m")
        sys.exit(1)
        
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
                transcript_type = attributes.get('transcript_biotype', attributes.get('transcript_type', ''))

                if index_by in attributes:
                    index_key = attributes[index_by]
                else:
                    continue
                
                start = int(columns[3]) 
                end = int(columns[4])  
                length = end - start 

                if index_key not in transcripts:
                        transcripts[index_key] = []
                    
                if not protein_coding:
                    transcripts[index_key].append((transcript_id, length))
                    
                elif protein_coding and transcript_type == 'protein_coding':
                    transcripts[index_key].append((transcript_id, length))
    
    return transcripts


def parse_gtf_small(gtf_file, index_by, protein_coding=False):

    if index_by == 'gene_id' or index_by == 'gene_name':
        pass
    else:
        logging.error("\033[91mError: Please specify the feature type in the provided matrix (transcript_id, gene_id or gene_name) using '--features' argument.\033[0m")
        sys.exit(1)
        
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
                transcript_type = attributes.get('transcript_biotype', attributes.get('transcript_type', ''))

                index_key = attributes[index_by]
                
                start = int(columns[3]) 
                end = int(columns[4])  
                current_length = end - start 

                if index_key in transcripts:
                    if len(transcripts[index_key]) > 1:
                        _, length = transcripts[index_key]
                        
                        if current_length < length:
                            if not protein_coding:
                                transcripts[index_key] = [transcript_id, length]
                                
                            elif protein_coding and transcript_type == 'protein_coding':
                                transcripts[index_key] = [transcript_id, length]

                else:
                    transcripts[index_key] = [transcript_id, current_length]
   
    transcripts = {key:[transcript_id[0]]  for key, transcript_id in transcripts.items()}
    
    return transcripts


def parse_gtf_old(gtf_file, index_by, protein_coding=True):
    
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


def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,pbar = True, pbar_update = 500,  *arg, **kwargs):
    executor = ProcessPoolExecutor(n_process)

    max_queue = n_process * 2
    if pbar:
        pbar = tqdm(unit = 'read', desc='Processed')

    futures = {}
    n_job_in_queue = 0
    while True:
        while n_job_in_queue < max_queue:
            i = next(iterator, None)
            if not i:
                break
            futures[executor.submit(func, i, *arg, **kwargs)] = None
            n_job_in_queue += 1

        job = next(as_completed(futures), None)

        if job is None:
            break
        else:
            n_job_in_queue -= 1
            pbar.update(pbar_update)
            yield job
            del futures[job]