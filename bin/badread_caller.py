"""
This script is based on the original work from https://github.com/youyupei/SLSim/
All original rights are acknowledged.
"""
import sys, os
import argparse
from tqdm import tqdm
import logging
from badread.simulate import sequence_fragment, ErrorModel, QScoreModel, Identities
from toolkit import multiprocessing_submit, read_template

YELLOW = "\033[93m"
GRAY = "\033[90m"
RESET = "\033[0m"

logging.basicConfig(level=logging.INFO, format=GRAY+ '%(message)s' +RESET)

def setup_badread_command(parent_parser):
    description="Nanopore sequencing error simulator."
    if parent_parser:
        parser = argparse.ArgumentParser(parents=[parent_parser], 
                                        description=description)
    else :
        parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-t', '--template', type=str,required=True,
                    help='Template fastq file.')
    parser.add_argument('--badread-error-model', type=str, default='nanopore2020',
                        help="Error model from babread.")
    parser.add_argument('--badread-qscore-model', type=str, default='nanopore2020',
                        help="Error model from babread.")
    parser.add_argument('--badread-identity', type=str, default='87.5,5,97.5',
                        help="Identity/accuracy (in percentage) parameter pass to badread: format mean,st,max.")              
    parser.add_argument('-o','--output',  type=str, default='sim_read.fq',
                        help="Filename of simulated reads with errors.")
    parser.add_argument('--batch-size',  type=int, default=500,
                        help= "Batch size")
    parser.add_argument('--thread',  type=int, default = 1, help="Thread")
    return parser


def write_fastq(output, name, seq, quals):
    with open(output, 'a') as outFasta:
        outFasta.write(f"@{name}\n"
                        f"{seq}\n"
                        "+\n"
                        f"{quals}\n")


def error_simulator(template, identities, err_model, q_model):
    reads=[]
    for name, read in template:
        err_seq, quals, _, _= sequence_fragment(read, identities, err_model, q_model)
        reads.append((name, err_seq, quals))
    return reads


def badread_caller(args):
    stderr = sys.stderr

    logging.info("_____________________________")
    logging.info("Input file name : %s" , args.template)
    logging.info("Badread error model: %s", args.badread_error_model)
    logging.info("Badread Qscore model: %s", args.badread_error_model)
    logging.info("Badread identity parameters: %s", args.badread_identity)
    logging.info("Output file name : %s", args.output)
    logging.info("Threads: %d", args.thread)
    logging.info("_______________________________")

    if os.path.exists(args.output):
        sys.exit(f'Oups! Output file {args.output} already exists.')

    error_model = ErrorModel(args.badread_error_model, stderr)
    qscore_model = QScoreModel(args.badread_qscore_model, stderr)

    badread_identity = args.badread_identity.split(',')
    mean, sd, mx = [float(x) for x in badread_identity]
    identities = Identities(mean, sd, mx, stderr).get_identity()

    template_chunks = read_template(args.template,
                                    thread=args.thread,
                                    batch_size=args.batch_size)

    results = multiprocessing_submit(error_simulator,
                                    template_chunks, 
                                    identities=identities, 
                                    err_model=error_model, 
                                    q_model=qscore_model,
                                    n_process=args.thread, 
                                    pbar_update=args.batch_size)                   

    for f in results:
        reads_res = f.result()
        for name, err_seq, quals in reads_res:
            write_fastq(args.output, name, err_seq, quals)

if __name__ == '__main__':
    args = setup_badread_command(None).parse_args()
    badread_caller(args)