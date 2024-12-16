"""
This script is based on the original work from https://github.com/youyupei/SLSim/
All original rights are acknowledged.
"""
import sys, os
import argparse
from tqdm import tqdm
import logging
from badread.simulate import sequence_fragment, ErrorModel, QScoreModel, Identities
from toolkit import read_template, process_template_chunk
from concurrent.futures import ProcessPoolExecutor


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
    parser.add_argument('--batch-size',  type=int, default=5000,
                        help= "Batch size")
    parser.add_argument('--thread',  type=int, default = 1, help="Thread")
    return parser


def write_fastq(output, name, seq, quals):
    with open(output, 'a') as outFasta:
        outFasta.write(f"@{name}\n"
                        f"{seq}\n"
                        "+\n"
                        f"{quals}\n")

def write_fastq_batch(output, batch_reads):
    with open(output, 'a') as outFasta:
        for name, seq, quals in batch_reads:
            outFasta.write(f"@{name}\n{seq}\n+\n{quals}\n")


def ln_error_simulator(template, identities, err_model, q_model, output):
    with open(output, 'w') as outFasta:
        for name, read in tqdm(template, total=100000):
            err_seq, quals, _, _= sequence_fragment(read, identities, err_model, q_model)
            outFasta.write(f"@{name}\n{err_seq}\n+\n{quals}\n")


def badread_caller(args):
    stderr = sys.stderr

    logging.info("_____________________________")
    logging.info("Input file name : %s" , args.template)
    logging.info("Badread error model: %s", args.badread_error_model)
    logging.info("Badread Qscore model: %s", args.badread_error_model)
    logging.info("Badread identity parameters: %s", args.badread_identity)
    logging.info("Output file name : %s", args.output)
    logging.info("Batch size : %s", args.batch_size)
    logging.info("Threads: %d", args.thread)
    logging.info("_______________________________")

    if os.path.exists(args.output):
        sys.exit(f'Oups! Output file {args.output} already exists.')

    error_model = ErrorModel(args.badread_error_model, stderr)
    qscore_model = QScoreModel(args.badread_qscore_model, stderr)

    badread_identity = args.badread_identity.split(',')
    mean, sd, mx = [float(x) for x in badread_identity]
    identities = Identities(mean, sd, mx, stderr).get_identity()


    if args.thread  > 1 :
    
        template_chunks = read_template(args.template,
                                        thread=args.thread,
                                        batch_size=args.batch_size)

        with ProcessPoolExecutor(max_workers=args.thread) as executor:
            futures = [
                executor.submit(process_template_chunk, chunk, identities, error_model, qscore_model)
                for chunk in template_chunks
            ]

            for future in tqdm(futures, desc="Processing"):
                batch_results = future.result()
                write_fastq_batch(args.output, batch_results)

    else:
        template=read_template(args.template,
                                        thread=1,
                                        batch_size=0)
        ln_error_simulator(template, 
                                identities=identities, 
                                err_model=error_model, 
                                q_model=qscore_model,
                                output=args.output)


if __name__ == '__main__':
    args = setup_badread_command(None).parse_args()
    badread_caller(args)