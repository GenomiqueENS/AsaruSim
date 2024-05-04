"""
This script is based on the original work from https://github.com/youyupei/SLSim/
All original rights are acknowledged.
"""
import sys
import argparse
import Bio.SeqIO
import textwrap
from tqdm import tqdm
from badread.simulate import sequence_fragment, ErrorModel, QScoreModel, Identities
import os
import multiprocessing as mp
import random
import numpy as np

import helper
#import config

def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        Descript come later
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Required positional argument
    parser.add_argument('-t', '--template-fasta', type=str,required=True,
                        help='Template fastq file.')
    # Required positional argument
    parser.add_argument('--badread-error-model', type=str, default='nanopore2020',
                        help=textwrap.dedent('''
                        Error model from babread.[reference needed]
                        '''
                        ))
    parser.add_argument('--badread-qscore-model', type=str, default='nanopore2020',
                        help=textwrap.dedent('''
                        Error model from babread.[reference needed]
                        '''
                        ))
    parser.add_argument('--badread-identity', type=str, default='87.5,5,97.5',
                        help=textwrap.dedent('''
                        Identity/accuracy (in percentage) parameter pass to badread: format mean,st,max.
                        '''
                        ))

                        
    parser.add_argument('-o','--output-filename',  type=str, default='sim_read.fq',
                        help=textwrap.dedent('''
                        Filename of simulated reads with errors.
                        '''
                        ))
    parser.add_argument('--batch-size',  type=int, default=500,
                        help=textwrap.dedent('''
                        Batch size
                        '''
                        ))
    parser.add_argument('--thread',  type=int, default = mp.cpu_count()-1,
                        help=textwrap.dedent('''
                        Thread
                        '''
                        ))
    parser.add_argument('--test-error-dist',  action='store_true',
                        help=textwrap.dedent('''
                        Test distribution of badread
                        '''
                        ))
    

    args = parser.parse_args()

    if not args.test_error_dist:
        helper.check_exist([args.template_fasta])
    return args


def sim_read(perfect_read, args, error_model, qscore_model, identities):
    """Simulate error into perfect read using Badread error and qscore model

        read (str): perfect reads
    """
    #random.seed("123")
    #np.random.seed("123")
    output=sys.stderr
    seq, quals, actual_identity, identity_by_qscores = \
                    sequence_fragment(perfect_read, identities.get_identity(), 
                                        error_model, qscore_model)
    return seq, quals



def read_batch_generator(fasta_fns, batch_size):
    """Generator of barches of reads from list of fasta files
    Args:
        fasta_fns (list): fasta filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fasta_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fasta = Bio.SeqIO.parse(handle, "fasta")
                read_batch = helper.batch_iterator(fasta, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            fasta = Bio.SeqIO.parse(fn, "fasta")
            read_batch = helper.batch_iterator(fasta, batch_size=batch_size)
            for batch in read_batch:
                yield batch


def sim_read_batch(read_batch, args, error_model, qscore_model, identities):
    fastq_lines = []
    for read in read_batch:
        seq, qscore = sim_read(str(read.seq), args, error_model, 
                                        qscore_model, identities)
        # write read.id, se quals\
        fastq_lines.extend([f'@{read.id}', seq, '+' , qscore])
    return fastq_lines


def main(output=sys.stderr):
    args = parse_arg()

    # load BadRead models
    error_model = ErrorModel(args.badread_error_model, output)
    qscore_model = QScoreModel(args.badread_qscore_model, output)
    
    # error rate distribution, .get_identity method generate random error rate 
    # from beta distribution
    mean, sd, maxi = [float(x) for x in args.badread_identity.split(',')]
    identities = Identities(mean, sd, maxi, output)
    if args.test_error_dist:
        sys.exit(0)

    # input template file
    if os.path.isdir(args.template_fasta):
        fns = helper.get_files(fastq_dir, ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz',
                                    '*.fasta', '*.fa', '*.fasta.gz', '*.fa.gz'])
    elif os.path.isfile(args.template_fasta):
        fns = [args.template_fasta]
    else:
        helper.err_msg("Invalid input of template fasta file")
        sys.exit(1)

    read_batchs = read_batch_generator(fns, args.batch_size)
    

    # single threading mode (good for debuging)
    if args.thread == 1:
        pbar = tqdm(unit = 'read', desc='Processed')
        for batch in read_batchs:
            fastq_records = sim_read_batch(batch, args, error_model, 
                                                qscore_model, identities)
            
            with open(args.output_filename, 'w') as f:
                f.write('\n'.join(fastq_records) + '\n')
            pbar.update(args.batch_size)
        
        helper.green_msg(f'Simulation finished!')
        return None
    else:
        # multiprocessing mode
        pass
    
    rst_futures = helper.multiprocessing_submit(sim_read_batch,
                    read_batchs, n_process=args.thread, 
                    pbar_update=args.batch_size,
                    args=args,
                    error_model=error_model,
                    qscore_model=qscore_model,
                    identities=identities)
    
    for idx, f in enumerate(rst_futures):
        fastq_records = f.result() #write rst_df out
        if idx == 0:
            with open(args.output_filename, 'w') as f:
                f.write('\n'.join(fastq_records)+'\n')
        else:
            with open(args.output_filename, 'a') as f:
                f.write('\n'.join(fastq_records)+'\n')

    helper.green_msg(f'Simulation finished!')
    return None    


if __name__ == '__main__':
    main()
