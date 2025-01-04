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
from intron_retention import get_sequence_with_intron, gtf_to_ref_structure, read_IR_model
from intron_retention import get_unspliced_sequence
import pysam

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
    parser.add_argument('-g','--gtf', type=str, help="Path to the annotation GTF file. or annotation file in GFF format containing the intron coordinate information. Required when intron retention mode is enabled.")
    parser.add_argument('-b','--unfilteredBC', type=str, help="Path to the unfiltered barcode counts CSV file.")

    parser.add_argument('--ref_genome', type=str, help="Path to the reference genome file. Required when intron retention mode is enabled.")

    parser.add_argument('-o','--outFasta', type=str, default="template.fa", help="Output file name for the generated template sequences.")
    parser.add_argument('-t','--threads', type=int, default=1, help="Number of threads to use.")
    parser.add_argument('-f','--features', type=str, default="transcript_id", help="Feature rownames in input matrix.")
    parser.add_argument('--on_disk', action='store_true', help="Whether to use Faidx to index the FASTA or load the FASTA in memory.")
    parser.add_argument('-c','--amp', type=int, default=0, help="amplification rate.")
    parser.add_argument('--unspliced_ratio', type=float, default=0.0, help="Unspliced transcript ratio.")

    parser.add_argument('--full_length', action='store_true', help="Simulate a full length untruncated transcripts.")
    parser.add_argument('--length_dist', type=str, default="0.37,0.0,824.94", help="amplification rate.")
    parser.add_argument('--truncation_model', type=str, default="truncation_default_model.csv", help="truncation model.")

    parser.add_argument('--intron_retention', action='store_true', help="Simulate intron retention events.")
    parser.add_argument('--IR_model', type=str, default='SC3pv3_GEX_Human_IR_markov_model', help="Intron retention model.")

    parser.add_argument('--adapter', type=str, default="ATGCGTAGTCAGTCATGATC", help="Adapter sequence.")
    parser.add_argument('--TSO', type=str, default="ATGCGTAGTCAGTCATGATC", help="TSO sequence.")
    parser.add_argument('--len_dT', type=int, default=15, help="Poly-dT sequence.")
    parser.add_argument('--log', type=str, help="Path to the log file CSV.")
    return parser


class Transcriptome:
    def __init__(self, fasta_file, on_disk):
        self.on_disk = on_disk
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
    def __init__(self, transcriptome, adapter, len_dT, TSO, outFasta, threads, 
                 amp, truncation_model, full_length, length_dist, 
                 intron_retention, dict_ref_structure, IR_markov_model, genome_fai, unspliced_ratio):
        self.COMPLEMENT_MAP  = str.maketrans('ATCG', 'TAGC')
        self.transcriptome = transcriptome
        self.full_length = full_length
        self.length_dist = length_dist
        self.adapter = adapter
        self.len_dT = int(len_dT)
        self.TSO = TSO
        self.amp = amp
        self.truncation_model = truncation_model
        self.intron_retention=intron_retention
        self.dict_ref_structure=dict_ref_structure
        self.IR_markov_model=IR_markov_model
        self.genome_fai=genome_fai
        self.unspliced_ratio = unspliced_ratio
        self.outFasta = outFasta
        self.threads = threads
        self.counter = 0
        self.unfound = 0
        self.total_umi_counts = 0

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
            
        if not self.full_length:
            left_probs = list(self.truncation_model.left)
            right_probs = list(self.truncation_model.right)
        
        with open(filename, mode) as fasta:
            cnt = 0
            unfound = 0
            total_umi_counts = 0
            trns_ids = random.choices(list(self.transcriptome.transcripts.keys()), k=umi_counts)
            for _, trns in enumerate(trns_ids):
                if self.unspliced_ratio > 0 and random.random() < self.unspliced_ratio:
                        cDNA = get_unspliced_sequence(trns,
                                                      self.genome_fai,
                                                      self.dict_ref_structure)
                elif self.intron_retention:
                        i_retained, cDNA = get_sequence_with_intron(trns, 
                                                                    self.genome_fai,
                                                                    self.dict_ref_structure,
                                                                    self.IR_markov_model)
                        if not i_retained:
                            cDNA = self.transcriptome.get_sequence(trns)
                else:
                    cDNA = self.transcriptome.get_sequence(trns)
                if cDNA:
                    if not self.full_length:
                        cDNA = truncate_cDNA(cDNA, probs_3p=right_probs, probs_5p=left_probs)
                    umi = self.generate_random_umi()
                    total_umi_counts += 1
                    dT = "T"*int(np.random.normal(self.len_dT, 3))
                    seq = f"{self.adapter}{cell}{umi}{dT}{cDNA}{self.TSO}"
                    
                    if random.randint(0, 1) == 1:
                        seq = self.complement(seq)[::-1]
                    ir_tag = "-IR" if i_retained else ""

                    fasta.write(f">{cell}-{umi}-{trns}-{cnt}{ir_tag}\n"
                                    f"{seq}\n")
                    cnt += 1
                    if self.amp > 0:
                        for j in range(np.random.poisson(self.amp)):
                            fasta.write(f">{cell}-{umi}-{trns}-{cnt}{ir_tag}\n"
                                        f"{seq}\n")
                            cnt += 1

                else : unfound += 1

            self.total_umi_counts += total_umi_counts
            self.counter += cnt
            self.unfound += unfound


    def make_template(self, cell, umi_counts, transcripts_index, features):
        if self.threads == 1:
            filename = self.outFasta
            mode = 'a'
        else:
            filename = f"tmp_cells/{cell}.fa"
            mode = 'w'

        transcript_id = True if features=="transcript_id" else False

        if not self.full_length:
            left_probs = list(self.truncation_model.left)
            right_probs = list(self.truncation_model.right)
        
        with open(filename, mode) as fasta:
            idx = 0
            unfound = 0
            total_umi_counts = 0
            i_retained = False
            cDNA = None
            for trns, count in umi_counts.items():
                count=int(count)
                if count > 0:
                    trns = trns if transcript_id else fetch_transcript_id_by_gene(trns, transcripts_index)
                    if trns:
                        if self.unspliced_ratio > 0 and random.random() < self.unspliced_ratio:
                                cDNA = get_unspliced_sequence(trns,
                                                              self.genome_fai,
                                                              self.dict_ref_structure)
                        elif self.intron_retention:
                            i_retained, cDNA = get_sequence_with_intron(trns, 
                                                                        self.genome_fai,
                                                                        self.dict_ref_structure,
                                                                        self.IR_markov_model)
                            if not i_retained:
                                cDNA = self.transcriptome.get_sequence(trns)
                        else:
                            cDNA = self.transcriptome.get_sequence(trns)

                    if cDNA:
                        if not self.full_length:
                            cDNA = truncate_cDNA(cDNA, probs_3p=right_probs, probs_5p=left_probs)
                        for _ in range(count):
                            umi = self.generate_random_umi()
                            total_umi_counts += 1
                            dT = "T"*int(np.random.normal(self.len_dT, 3))
                            seq = f"{self.adapter}{cell}{umi}{dT}{cDNA}{self.TSO}"
                            if random.randint(0, 1) == 1:
                                seq = self.complement(seq)[::-1]
                            ir_tag = "-IR" if i_retained else ""

                            fasta.write(f">{cell}-{umi}-{trns}-{idx}{ir_tag}\n"
                                                f"{seq}\n")
                            idx += 1
                            if self.amp > 0:
                                for _ in range(np.random.poisson(self.amp)):
                                    fasta.write(f">{cell}-{umi}-{trns}-{idx}{ir_tag}\n"
                                                f"{seq}\n")
                                    idx += 1

                    else : unfound += 1

            self.total_umi_counts += total_umi_counts
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
    logging.info("Amplification rate: %s", args.amp)
    logging.info("Adapter sequence: %s", args.adapter)
    logging.info("TSO sequence: %s", args.TSO)
    logging.info("Truncation model: %s", args.truncation_model)
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
            logging.error("\033[91mError: Please provide the path to the count matrix file using the '--matrix' or cell barcode counts using '--unfilteredBC' option.\033[0m")
            sys.exit(1)
            
    if args.full_length:
        length_dist = None
        truncation_model = None
        
    else :
        try:
            shape, loc, scale = args.length_dist.split(",")
            length_dist = [float(x) for x in [shape, loc, scale]]
        except:
            length_dist = [0.37, 0.0, 824.94]
    
    if not args.full_length:
        truncation_model = pd.read_csv(args.truncation_model)

    dict_ref_structure = None
    IR_markov_model = None
    genome_fai = None

    if args.intron_retention or args.unspliced_ratio > 0 :
        dict_ref_structure = gtf_to_ref_structure(args.gtf)

        if args.ref_genome:
            genome_fai = pysam.Fastafile(args.ref_genome)
        else:
            logging.error("\033[91mError: Please provide the path to the reference genome file using the '--ref_genome' option.\033[0m")
            sys.exit(1)
        
        if args.intron_retention:
            IR_markov_model = read_IR_model(args.IR_model)

        

    transcriptome = Transcriptome(args.transcriptome, args.on_disk)
    generator = TemplateGenerator(transcriptome = transcriptome, 
                                  adapter = args.adapter, 
                                  len_dT = args.len_dT, 
                                  TSO = args.TSO,
                                  full_length = args.full_length,
                                  length_dist = length_dist,
                                  outFasta = args.outFasta, 
                                  threads = args.threads,
                                  amp = args.amp,
                                  truncation_model=truncation_model,
                                  intron_retention=args.intron_retention,
                                  dict_ref_structure=dict_ref_structure,
                                  IR_markov_model=IR_markov_model,
                                  genome_fai=genome_fai,
                                  unspliced_ratio=args.unspliced_ratio)
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
                generator.make_template(cell, 
                                        umi_counts, 
                                        transcripts_index, 
                                        args.features)
            
        if args.unfilteredBC:
            logging.info("Unfiltered BC counts ...")
            unfiltered_bc = unfiltered_bc[~unfiltered_bc['BC'].isin(matrix.columns)]
            for _, row in tqdm(unfiltered_bc.iterrows(), total=len(unfiltered_bc)):
                generator.make_random_template(row['BC'], 
                                               row['counts'])
                
    else:
        if not os.path.exists("tmp_cells"):
            os.makedirs("tmp_cells")
        
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            if args.matrix:
                logging.info("Filtered Matrix...")
                for cell in tqdm(matrix.columns):
                    umi_counts = matrix[cell]
                    future = executor.submit(generator.make_template, 
                                             cell, 
                                             umi_counts, 
                                             transcripts_index, 
                                             args.features)
                    futures.append(future)
    
                for future in futures:
                    future.result()
                
            futures = []
            if args.unfilteredBC:
                logging.info("Unfiltered BC counts ...")
                unfiltered_bc = unfiltered_bc[~unfiltered_bc['BC'].isin(matrix.columns)]
                for idx, (_, row) in tqdm(enumerate(unfiltered_bc.iterrows()), total=len(unfiltered_bc)):
                    future = executor.submit(generator.make_random_template, 
                                             row['BC'], 
                                             row['counts'])

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
        "Simulated UMI counts": generator.total_umi_counts,
        "Unknown transcript counts": generator.unfound,
        "Number of reads": generator.counter}
        log_df = pd.DataFrame(log_df, index=[0])
        log_df.to_csv(args.log, index=False)
    
    logging.info("Completed successfully. Have a great day!")
    logging.info("Stats : "+
                 "\nSimulated Cell BC: "+
                 str(len(matrix.columns))+
                 "\nSimulated Filtered-Out Cell BC: "+
                 str(count_unfiltered_bc)+
                 "\nSimulated UMI: "+
                 str(generator.total_umi_counts)+
                 "\nSimulated reads: "+
                 str(generator.counter)+
                 "\nUnknown transcript counts: "+str(generator.unfound))
    
if __name__ == "__main__":
    args = setup_template_parameters(None).parse_args()
    if args:
        template_maker(args)
