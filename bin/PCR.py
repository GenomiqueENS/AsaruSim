import argparse
import logging
import numpy as np
import random
import os, sys, math
from toolkit import multiprocessing_submit

YELLOW = "\033[93m"
GRAY = "\033[90m"
RESET = "\033[0m"

logging.basicConfig(level=logging.INFO, format=GRAY+ '%(message)s' +RESET)
parser = argparse.ArgumentParser(description="Script for generating template sequences for scRNAseq simulation.")

parser.add_argument('-f','--template', type=str, required=True, help="Path to the template FASTA file.")
parser.add_argument('-o','--out', type=str, default="out.fa", help="Path to the output FASTA file.")
parser.add_argument('-c','--cycles', type=int, default=5, help="number of cycles.")
parser.add_argument('-d','--dup', type=float, default=0.7, help="duplication rate.")
parser.add_argument('-e','--error', type=float, default=0.00003, help="error rate.")
parser.add_argument('-t','--thread', type=int, default=4, help="number of threads to use.")
parser.add_argument('-b','--batch_size', type=int, default=500, help="batch size.")
parser.add_argument('-n','--totalNamber', type=int, default=None, help="total number of sequence to select from the finel pool.")
parser.add_argument('-s','--seed', type=int, default=2024, help="seed value.")

class Molecule ():
    def __init__(self, name, length, seq, root, inherited_mut):
        self.name = name
        self.length = length
        self.seq = seq
        self.root = root
        self.nucleotides = ['A', 'C', 'G', 'T'] 
        self.inherited_mut = inherited_mut

    def mutate(self, error_rate):
        """
        Apply random mutations based on the given error rate, directly altering the sequence.
    
        Args:
            mutations (dict): A dictionary mapping nucleotides to possible mutations.
            sequence (list): The current sequence of the molecule as a list of characters.
            error_rate (float): The mutation rate per base pair per cycle.
        """
        mutation_count = math.ceil(self.length * error_rate)
        if mutation_count > 0:
            mutation_positions = np.random.choice(self.length, size=mutation_count, replace=True)
    
            for pos in mutation_positions:
                if pos in self.inherited_mut:
                    current_base = self.inherited_mut[pos]
                else:
                    current_base = self.seq[int(pos)]
                    
                possible_mutations = [nuc for nuc in self.nucleotides if nuc != current_base]
                
                if possible_mutations:
                    self.inherited_mut[pos] = np.random.choice(possible_mutations)


    def duplicate(self, mol_counter):
        """
        Duplicates this molecule by creating a new instance with the same properties.
        Returns:
            Molecule: A new instance of Molecule with the same data.
        """
        return Molecule(
            name=str(self.root) + "_" + str(mol_counter),
            length=self.length,
            seq=self.seq,
            root=self.root,
            inherited_mut=self.inherited_mut.copy())

        
def sequencing(pool, outfile, totalNumber=None):
    if totalNumber:
        pool = random.sample(pool, totalNumber)
    with open(outfile, 'a') as outFasta:
        for mol in pool:
            root_seq = list(str(mol.seq))
            for pos in mol.inherited_mut:
                root_seq[pos] = str(mol.inherited_mut[pos])
            outFasta.write(f">{mol.name}\n"
                           f"{''.join(root_seq)}\n")


def PCR(template, cycles, dup_rate, error_rate, totalNumber, outfile):
    pool_list=[] 
    mol_counter=0

    for cycle in range(cycles):
        if mol_counter == 0:
            for name, seq in template:
                pool_list.append(Molecule(name=name, 
                                         length=len(seq),
                                         seq=seq,
                                         root=name, 
                                         inherited_mut=dict()))
        print("cycle : ", cycle+1)
        product = []
        
        if dup_rate < 1:
            dup_number = len(pool_list) * dup_rate
            pool_list = random.sample(pool_list, int(dup_number))
        
        for mol in pool_list:
            child_mol = mol.duplicate(mol_counter)
            child_mol.mutate(error_rate)

            product.append(child_mol)
            mol_counter += 1

        pool_list.extend(product)

    sequencing(pool_list, outfile, totalNumber)


def read_template(fasta, thread, batch_size):
    """
    Function to read and extract sequences from a FASTQ file. It supports both plain text and gzipped files.
    Inputs:
    fastq: String, path to the FASTQ file.
    Returns: List of sequences encoded in ASCII.
    """
    if thread==1:
        with open(fasta, 'rt') as f:
            while True:
                header = f.readline()
                if not header: break
                name = header[1:].strip()
                seq = f.readline().strip()
                yield name, seq
    else:
        batch = []
        with open(fasta, 'r') as f:
                while True:
                    header = f.readline()
                    if not header: break
                    name = header[1:].strip()
                    seq = f.readline().strip()
                    batch.append([name, seq])
                    if len(batch) == batch_size:
                        yield batch
                        batch = []
                if len(batch) > 0:
                    yield batch


def main():
    args = parser.parse_args()
    logging.info("_____________________________")
    logging.info("Template file name : %s" , args.template)
    logging.info("PCR cycles         : %s" , args.cycles)
    logging.info("Error rate         : %s" , args.error)
    logging.info("Duplication rate   : %s" , args.dup)
    logging.info("Output file name   : %s" , args.out)
    logging.info("Thread             : %s" , args.thread)
    logging.info("Batch size         : %s" , args.batch_size)
    logging.info("_______________________________")

    random.seed(args.seed)
    np.random.seed(args.seed)

    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        handler.setFormatter(logging.Formatter(YELLOW+'%(asctime)s - %(levelname)s' +RESET+ ' - %(message)s'))


    if os.path.exists(args.out):
        logging.error(f"File '{args.out}' already exists.")
        sys.exit()

    if args.dup >1 or args.dup < 0:
        logging.error(f"Invalid number of --dup. It takes values between 0 and 1.")
        sys.exit()
    
    if args.thread == 1:
        template = read_template(args.template,
                                 thread=args.thread,
                                 batch_size=args.batch_size)

        pcr_pool = PCR(template,
                       cycles=args.cycles, 
                       dup_rate=args.dup,
                       error_rate=args.error,
                       totalNumber=args.totalNamber,
                       outfile=args.out)
    
    else:
        template_chunks = read_template(args.template,
                                        thread=args.thread,
                                        batch_size=args.batch_size)

        results = multiprocessing_submit(PCR,
                                        template_chunks, 
                                        cycles=args.cycles,
                                        dup_rate=args.dup,
                                        error_rate=args.error,
                                        totalNumber=args.totalNamber,
                                        outfile=args.out,
                                        n_process=args.thread, 
                                        pbar_update=args.batch_size)
                            
        pcr_pool = []
        for f in results:
            pass
            
    logging.info("DONE!_______________________________")

if __name__ == "__main__":
    main()  