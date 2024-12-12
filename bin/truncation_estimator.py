#import pysam
import numpy as np
import pandas as pd
import argparse
import logging
import sys
import gzip

YELLOW = "\033[93m"
GRAY = "\033[90m"
RESET = "\033[0m"

logging.basicConfig(level=logging.INFO, format=GRAY+ '%(message)s' +RESET)

def setup_truncation_parameters(parent_parser):
    description="truncation model estimator."
    if parent_parser:
        parser = argparse.ArgumentParser(parents=[parent_parser], 
                                        description=description)
    else :
        parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-p','--paf', type=str, help="Path to the alignement PAF file.")
    parser.add_argument('-b','--bam', type=str, help="Path to the alignement BAM/SAM file.")
    return parser


def campute_truncation_from_bam(bam):

    mode =  "r" if bam.endswith(".sam") else "rb"
    bamfile = pysam.AlignmentFile(bam, mode)

    left_truncations = []
    right_truncations = []
    for alignment in bamfile:
        if alignment.is_secondary or alignment.is_supplementary or alignment.is_unmapped:
              continue
    
        transcript_len = bamfile.get_reference_length(bamfile.get_reference_name(alignment.reference_id))
        left_truncations.append(float(alignment.reference_start) / float(transcript_len))
        right_truncations.append(float(transcript_len - alignment.reference_end) / float(transcript_len))
    return left_truncations, right_truncations


def compute_truncation_from_paf(paf_file):
    left_truncations = []
    right_truncations = []
    
    with gzip.open(paf_file, "r") as file:
        logging.info("reading alignement file...")
        for line in file:
            fields = line.decode('utf-8').strip().split('\t')
            
            optional_tags = {tag.split(':')[0]: tag.split(':')[2] for tag in fields[12:] if ':' in tag}
            if optional_tags.get('tp') != 'P':  # Skip if not primary alignment
                continue

            query_length = int(fields[1])
            query_start = int(fields[2]) #0-based
            query_end = int(fields[3])
            
            left_truncations.append(float(query_start) / float(query_length))
            right_truncations.append(float(query_length - query_end) / float(query_length))
    
    return left_truncations, right_truncations


def comp_probs(truncations):
    logging.info("Computing truncation probabilities...")
    bins = [0.0] + [0.005 + 0.01 * i for i in range(0,100)] + [1.0]
    hist = np.histogram(truncations, bins=bins, density=True)
    
    probs = []
    for i, w in enumerate(hist[0]):
        probs.append(w * (bins[i+1] - bins[i]))
    return probs


def truncation_estimator(args):
    if args.paf:
        left_truncations, right_truncations = compute_truncation_from_paf(args.paf)
    elif args.paf:
        left_truncations, right_truncations = compute_truncation_from_paf(args.bam)
    else:
        logging.error("\033[91mError: Please provide the path to the PAF file using the '--paf' argument.\033[0m")
        sys.exit(1)

    left_prob = comp_probs(left_truncations)
    right_prob = comp_probs(right_truncations)

    df = pd.DataFrame({
    "left": left_prob,
    "right": right_prob
    })
    df.to_csv("truncation_model.csv")

if __name__ == "__main__":
    args = setup_truncation_parameters(None).parse_args()
    if args:
        truncation_estimator(args)