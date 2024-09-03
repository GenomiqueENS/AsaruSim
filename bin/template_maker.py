# Author: Ali Hamraoui
# Date: 25-04-2024

import random
import glob
import os
import argparse
import logging
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from pyfaidx import Fasta
from toolkit import *

YELLOW = "\033[93m"
GRAY = "\033[90m"
RESET = "\033[0m"

logging.basicConfig(level=logging.INFO, format=GRAY+ '%(message)s' +RESET)

def setup_template_parameters(parent_parser):
    description="Nanopore template maker."
    if parent_parser:
        parser = argparse.ArgumentParser(parents=[parent_parser], 
                                        description=description)
    else :
        parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-r','--transcriptome', type=str, help="Path to the transcriptome FASTA file.")
    parser.add_argument('-m','--matrix', type=str, help="Path to the SPARSim matrix CSV file.")
    parser.add_argument('-g','--gtf', type=str, help="Path to the annotation GTF file.")
    parser.add_argument('-b','--unfilteredBC', type=str, help="Path to the unfiltered barcode counts CSV file.")

    parser.add_argument('-o','--outFasta', type=str, default="template.fa", help="Output file name for the generated template sequences.")
    parser.add_argument('-t','--threads', type=int, default=1, help="Number of threads to use.")
    parser.add_argument('-f','--features', type=str, default="transcript_id", help="Feature rownames in input matrix.")
    parser.add_argument('--on_disk', action='store_true', help="Whether to use Faidx to index the FASTA or load the FASTA in memory.")
    parser.add_argument('-c','--amp', type=int, default=1, help="amplification rate.")

    parser.add_argument('--full_length', action='store_true', help="Simulate a full length transcripts.")
    parser.add_argument('--length_dist', type=str, default="0.37,0.0,824.94", help="amplification rate.")

    parser.add_argument('--adapter', type=str, default="ATGCGTAGTCAGTCATGATC", help="Adapter sequence.")
    parser.add_argument('--TSO', type=str, default="ATGCGTAGTCAGTCATGATC", help="TSO sequence.")
    parser.add_argument('--len_dT', type=str, default=15, help="Poly-dT sequence.")
    parser.add_argument('--log', type=str, help="Path to the log file CSV.")
    return parser


class Transcriptome:
    def __init__(self, fasta_file):
        self.on_disk = True
        self.transcripts = self.load(fasta_file)

    def load(self, fasta_file):
        if self.on_disk:
            filtered_transcripts = Fasta(fasta_file, key_function=custom_key_function)
        else : 
            transcripts = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
            filtered_transcripts = {key.split(".")[0]: value.seq for key, value in transcripts.items()}
        return filtered_transcripts

    def get_sequence(self, transcript_id):
        if self.on_disk:
            if transcript_id in self.transcripts:
                sequence = self.transcripts[transcript_id]
                return sequence[:]
            else :
                return ""
        else:
            return str(self.transcripts.get(transcript_id, ""))


class TemplateGenerator:
    def __init__(self, transcriptome, adapter, len_dT, TSO, outFasta, threads, amp, full_length, length_dist):
        self.COMPLEMENT_MAP  = str.maketrans('ATCG', 'TAGC')
        self.transcriptome = transcriptome
        self.full_length = full_length
        self.length_dist = length_dist
        self.adapter = adapter
        self.dT = "T"*int(len_dT)
        self.TSO = TSO
        self.amp = amp
        self.outFasta = outFasta
        self.threads = threads
        self.counter = 0
        self.unfound = 0

    def generate_random_umi(self, length=12):
        bases = 'ATCG'
        return ''.join(random.choice(bases) for _ in range(length))

    def complement(self, sequence):
        return sequence.translate(self.COMPLEMENT_MAP)
    
    def make_random_template(self, cell, umi_counts):
        if self.threads == 1:
            filename = self.outFasta
            mode = 'a'
        else:
            filename = f"tmp_cells/{cell}.fa"
            mode = 'w'
            
        max_len = 4000
        min_len = 3000
        p_cut = 0.88
        p_long = 0.9999
        
        with open(filename, mode) as fasta:
            idx = 1
            trns_ids = random.choices(list(self.transcriptome.transcripts.keys()), k=umi_counts)
            for idx, trns in enumerate(trns_ids):
                cDNA = self.transcriptome.get_sequence(trns)
                if cDNA:
                    if len(cDNA) > min_len:
                            cDNA = cut_sequence_proba(cDNA)
                                
                    umi = self.generate_random_umi()
                    seq = f"{self.adapter}{cell}{umi}{self.dT}{cDNA}{self.TSO}"
                    if random.randint(0, 1) == 1:
                        seq = self.complement(seq)[::-1]
                    for j in range(np.random.poisson(self.amp)):
                        fasta.write(f">{cell}_{umi}_{trns}_{idx}\n"
                                    f"{seq}\n")
                        idx += 1


    def make_template(self, cell, umi_counts, transcripts_index, features):
        if self.threads == 1:
            filename = self.outFasta
            mode = 'a'
        else:
            filename = f"tmp_cells/{cell}.fa"
            mode = 'w'
            
        transcript_id = True if features=="transcript_id" else False

        max_len = 2000
        min_len = 1500
        p_cut = 0.8
        p_long = 0.999
        
        with open(filename, mode) as fasta:
            idx = 0
            unfound = 0
            for trns, count in umi_counts.items():
                count=int(count)
                if count > 1:
                    trns = trns if transcript_id else fetch_transcript_id_by_gene(trns, transcripts_index)
                    cDNA = self.transcriptome.get_sequence(trns) if trns else None
                    if cDNA:   
                        for i in range(count):
                            umi = self.generate_random_umi()
                            seq = f"{self.adapter}{cell}{umi}{self.dT}{cDNA}{self.TSO}"
                            if random.randint(0, 1) == 1:
                                seq = self.complement(seq)[::-1]
                            for j in range(np.random.poisson(self.amp)):
                                fasta.write(f">{cell}_{umi}_{trns}_{idx}\n"
                                            f"{seq}\n")
                                
                                idx += 1
                    else : unfound += 1
            self.counter += idx
            self.unfound += unfound


def template_maker(args):
    
    logging.info("_____________________________")
    logging.info("Output file name : %s" , args.outFasta)
    logging.info("Threads: %d", args.threads)
    logging.info("Features to simulate: %s", args.features)
    logging.info("Path to the transcriptome file: %s", args.transcriptome)
    logging.info("Path to the matrix file: %s", args.matrix)
    logging.info("Path to the unfiltered barcode counts file: %s", args.unfilteredBC)
    logging.info("Adapter sequence: %s", args.adapter)
    logging.info("TSO sequence: %s", args.TSO)
    logging.info("Poly-dT sequence: %s", args.len_dT)
    logging.info("Load the FASTA in memory?: %s", args.on_disk)
    logging.info("_______________________________")
    
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        handler.setFormatter(logging.Formatter(GRAY+'%(asctime)s - %(levelname)s' +RESET+ ' - %(message)s'))

    if not args.matrix: 
        if args.unfilteredBC:
            logging.warning("Count matrix simulation skipped. Transcripts will be generated randomly.")
        else:
            logging.error("\033[91mError: Please provide the path to the count matrix file using the '--matrix' or cell barcode counts using '--unfilteredBC' argument.\033[0m")
            sys.exit(1)
            
    if args.full_length:
        length_dist = None
        
    else :
        try:
            shape, loc, scale = args.length_dist.split(",")
            length_dist = [float(x) for x in [shape, loc, scale]]
        except:
            length_dist = [0.37, 0.0, 824.94]


    transcriptome = Transcriptome(args.transcriptome)
    generator = TemplateGenerator(transcriptome = transcriptome, 
                                  adapter = args.adapter, 
                                  len_dT = args.len_dT, 
                                  TSO = args.TSO,
                                  full_length = args.full_length,
                                  length_dist = length_dist,
                                  outFasta = args.outFasta, 
                                  threads = args.threads,
                                  amp = args.amp)
    if args.matrix:
        matrix = pd.read_csv(args.matrix, header=0, index_col=0) 

    transcripts_index = None
    if args.features != "transcript_id":
        if args.gtf:
            logging.info("Parsing GTF file...")
            transcripts_index = parse_gtf(args.gtf, args.features)
            logging.info("Parsing GTF file completed successfully.")
        else:
            logging.error("\033[91mError: Please provide the path to the GTF file using the '--gtf' argument.\033[0m")
            sys.exit(1)
    
            
    logging.info("Creating template sequences...")
    if args.unfilteredBC:
        unfiltered_bc = pd.read_csv(args.unfilteredBC)
        unfiltered_bc.columns =['BC', 'counts']
    
    if args.threads == 1:
        if args.matrix:
            logging.info("Filtered Matrix...")
            for cell in tqdm(matrix.columns):
                umi_counts = matrix[cell]
                generator.make_template(cell, umi_counts, transcripts_index, args.features)
            
        if args.unfilteredBC:
            logging.info("Unfiltered BC counts ...")
            unfiltered_bc = unfiltered_bc[~unfiltered_bc['BC'].isin(matrix.columns)]
            for _, row in tqdm(unfiltered_bc.iterrows(), total=len(unfiltered_bc)):
                generator.make_random_template(row['BC'], row['counts'])
                
    else:
        if not os.path.exists("tmp_cells"):
            os.makedirs("tmp_cells")
        
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            if args.matrix:
                logging.info("Filtered Matrix...")
                for cell in tqdm(matrix.columns):
                    umi_counts = matrix[cell]
                    future = executor.submit(generator.make_template, cell, umi_counts, transcripts_index, args.features)
                    futures.append(future)
    
                for future in futures:
                    future.result()
                
            futures = []
            if args.unfilteredBC:
                logging.info("Unfiltered BC counts ...")
                unfiltered_bc = unfiltered_bc[~unfiltered_bc['BC'].isin(matrix.columns)]
                for idx, (_, row) in tqdm(enumerate(unfiltered_bc.iterrows()), total=len(unfiltered_bc)):
                    future = executor.submit(generator.make_random_template, row['BC'], row['counts'])

                for future in futures:
                    future.result()
                    
        with open(args.outFasta, "w") as outfile:
            for filename in glob.glob("tmp_cells/*.fa"):
                with open(filename) as infile:
                    outfile.write(infile.read())
                    
        logging.info("Removing temporary files...")
        for filename in glob.glob("tmp_cells/*.fa"):
            os.remove(filename)
            
    count_unfiltered_bc = len(unfiltered_bc) if args.unfilteredBC else 0

    if args.log:
        log_df = {
        "Simulated Cell BC": len(matrix.columns), 
        "Simulated Filtered-Out": count_unfiltered_bc,
        "Simulated UMI counts": generator.counter,
        "Unknown transcript counts": generator.unfound}
        log_df = pd.DataFrame(log_df, index=[0])
        log_df.to_csv(args.log, index=False)
    
    logging.info("Completed successfully. Have a great day!")
    logging.info("Stats : "+
                 "\nSimulated Cell BC: "+
                 str(len(matrix.columns))+
                 "\nSimulated Filtered-Out Cell BC: "+
                 str(count_unfiltered_bc)+
                 "\nSimulated UMI counts: "+
                 str(generator.counter)+
                 "\nUnknown transcript counts: "+str(generator.unfound))
    
if __name__ == "__main__":
    args = setup_template_parameters(None).parse_args()
    if args:
        template_maker(args)
