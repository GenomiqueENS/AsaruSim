# Author: Ali Hamraoui
# Date: 25-04-2024

import pandas as pd
import numpy as np
from scipy import stats
import argparse
import gzip
import sys
import logging

from libs import *

YELLOW = "\033[93m"
GRAY = "\033[90m"
RESET = "\033[0m"

logging.basicConfig(level=logging.INFO, format=GRAY+ '%(message)s' +RESET)

parser = argparse.ArgumentParser(description="Script for estimating real data identity and read lengths params.")
parser.add_argument('-r','--params', type=str, default="identity", help="read parametres to estimate: (identity or lengths)")
parser.add_argument('-p','--paf', type=str, help="a PAF (Pairwise Alignment Format) file")

parser.add_argument('-f','--fastq', type=str, help="trainning FASTQ file. can also be a compressed .gz FASTQ file")
parser.add_argument('-t','--thread', type=int, help="number of threads")

col_names = [
    'query_name', 'query_length', 'query_start', 'query_end',
    'strand', 'target_name', 'target_length', 'target_start',
    'target_end', 'residue_matches', 'alignment_block_length', 'mapping_quality',
    'col13', 'col14', 'col15', 'col16', 'col17', 'col18',
    'col19', 'col20', 'col21', 'col22', 'col23'
]

def calc_identity(paf_file):
    """
    Calculate the sequence identity percentage for alignments in a PAF file.

    Parameters:
    paf_file (str): The file path to the PAF file.

    Returns:
    Series: A pandas Series containing the identity percentages normalized by 100.01 for each alignment.
    """
    paf_df = pd.read_csv(paf_file, sep='\t',on_bad_lines='skip', header=None, names=col_names)

    paf_df['identity_percentage'] = (paf_df['residue_matches'] / paf_df['alignment_block_length']) * 100

    paf_df = paf_df[paf_df['identity_percentage']>79]

    return paf_df['identity_percentage']/100.01


def read_fastq(fastq):
    """
    Function to read and extract sequences from a FASTQ file. It supports both plain text and gzipped files.
    Inputs:
    fastq: String, path to the FASTQ file.
    Returns: List of sequences encoded in ASCII.
    """
    open_fn = gzip.open if fastq.endswith('.gz') else open
    with open_fn(fastq, 'rt') as f:
        records = []
        while True:
            header = f.readline()
            if not header: break
            seq = f.readline().strip().encode('ascii')
            f.readline()
            f.readline()
            records.append(seq)
    return records


def estimate_beta_parames(data):
    """
    Estimate the parameters of a beta distribution given data.

    Parameters:
    data (array-like): The data to fit the beta distribution to.

    Returns:
    tuple: A tuple containing the mean, standard deviation, and mode of the estimated beta distribution.
    """
    alpha, beta, loc, scale = stats.beta.fit(data, floc=0, fscale=1)
    mean = alpha / (alpha + beta)
    variance = (alpha * beta) / ((alpha + beta) ** 2 * (alpha + beta + 1))
    std_dev = np.sqrt(variance)
    if alpha > 1 and beta > 1:
        mode = (alpha - 1) / (alpha + beta - 2)
    else:
        mode = "NA" 
    return mean, std_dev, mode


def main():
    args = parser.parse_args()

    if args.params == "length":
        lengths_real = get_read_lengths(args.fastq, thread = args.thread, batch_size = 500)
        #lower_bound = np.percentile(lengths_real, 0.5)
        #upper_bound = np.percentile(lengths_real, 99.5)

        #filtered_lengths_real = [length for length in lengths_real if lower_bound < length < upper_bound]
        shape, loc, scale = stats.lognorm.fit(lengths_real, floc=0)
        print("{},{},{}".format(round(shape,2), round(loc,2), round(scale,2)))

    elif args.params == "identity":
        identities = calc_identity(args.paf)
        mean, std_dev, mode = estimate_beta_parames(identities)
        print("{},{},{}".format(round(mean,2)*100, round(std_dev,2)*100, round(mode,2)*100))
    
    else:
        logging.error("\033[91mError: Please choose parameres to estimate with -- params argument (identity or length).\033[0m")
        sys.exit(1)


if __name__ == "__main__":
    main()
 