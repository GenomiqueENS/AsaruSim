"""
This script is based on the original work from https://github.com/youyupei/SLSim/tree/main.
It has been modified and adapted
to meet specific requirements and functionalities not present in the original version.
All original rights are acknowledged.
"""

import pandas as pd
import numpy as np
import argparse
import textwrap
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint
from datetime import datetime
import os
import sys
import shutil
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
from multiprocessing import RLock

import helper as helper
import config

def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        This script can be used to simulate reads (perfect reads with no errors)
        from 10X single-cell Nanopore RNA sequencing data. 
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required positional argument
    parser.add_argument('-r', '--trans-ref', type=str,required=True,
                        help='Reference transcriptome to simulate reads from.')
    # Required positional argument
    parser.add_argument('-i', '--umi-count-file', type=str, default='../data/SR_bc.csv',
                        help='Filename of the BC UMI count csv file. ')
    parser.add_argument('--amp-rate', type=int, default=None,
                        help=textwrap.dedent('''
                        The number of reads (X) a unique molecule will be 1 if --amp-rate 
                        is not specified. Otherwise, X will be randomly drawn from poisson 
                        distribution with specified rate.
                        '''
                        ))
    parser.add_argument('-o','--output-filename',  type=str, default='template.fa',
                        help=textwrap.dedent('''
                        Filename of simulated reads.
                        '''
                        ))
    parser.add_argument('--mRNA-len',  type=int, default=200,
                        help=textwrap.dedent('''
                        Filename of simulated reads.
                        '''
                        ))
    parser.add_argument('-f', '--fit_distribution',  type=str, default=None,
                        help=textwrap.dedent('''
                        Filename of reference fastq to fit.
                        '''
                        ))
    parser.add_argument('--thread',  type=int, default=1,
                        help=textwrap.dedent('''
                        Filename of simulated reads.
                        '''
                        ))
    
    args = parser.parse_args()

    if args.amp_rate  and args.amp_rate <= 0:
        helper.err_msg(textwrap.dedent(
                f'''
                Invalid number of --amp-rate. It takes positive value only...
                '''))

    # check file 
    helper.check_exist([args.umi_count_file])
    return args

def polyT(length):
    return 'T'*length
    # do I need to randomize the length?

def random_seq(length):
    '''
    Generate random DNA seq with specified length
    '''
    return ''.join(np.random.choice(['A','G','C', 'T'],length))

def artificial_template(BC_list,
                        repeats,
                        tran_len,
                        tran_ref,
                        adaptor,
                        TSO,
                        output_fn,
                        fastq_dist,
                        thread,
                        amp_rate = None,
                        pbar_pos = 1
                        ):
    '''
    Generate a barcode 
        BC_list: list of unique BC
        repeats: list of count each BC is simulated
        adaptor/TSO: adapter/TSO sequence to attach before/after the barcodes, 
                            (use 10X 3' solution V3 adaptor by default)
        output_fn: output filename (append to the end of the file if exist)
        num_read: if specified, each UMI will be repected x time where x~pois(amp_rate)
    output:
        FASTA:
            seqname: barcode sequence
            seq: adapter + BC + 12nt random seq + 15nt polyT + 500 random sequence + TSO
    ''' 
    if repeats is None:
        repeats = [1]*len(BC_list)
    
    lock = tqdm.get_lock()

    # v3 adaptor
    if os.path.exists(output_fn):
        print(f"The file '{output_fn}' already exists. Exiting.")
        sys.exit()
    with open(output_fn, 'a') as f:
        if fastq_dist is None:
            random_frag = get_random_frag(tran_ref, tran_len)
        else:
            distribution_dict = helper.get_distribution(fastq_dist, thread=thread)
            print(sum(distribution_dict.values()))
            distribution_dict = scale_dist(distribution_dict, repeats)
            random_frag = get_random_frag_with_distribution(tran_ref, distribution_dict)
            print(sum(distribution_dict.values()))
        read_index = 0
        print(repeats.sum())
        #pbar = tqdm(total=sum(repeats), desc='Total UMI count')

        pbar = tqdm(total=sum(repeats), desc=f'Thread #{pbar_pos}: Total UMI count', position = pbar_pos, leave=False )
        for BC,count in zip(BC_list, repeats):
            for i in range(count):
                next_frag =  next(random_frag)
                UMI = random_seq(12)
                seq = adaptor+BC+UMI+ polyT(15) + next_frag[1] +TSO
                
                if amp_rate is None:
                    # seq (randomly reversed) 
                    read_index += 1
                    read_id = BC + '_' +  UMI + '_' + str(read_index) + '_' + str(next_frag[0])
                    f.write(">" + read_id + '\n')
                    if np.random.choice([False,True]):
                        f.write(seq + '\n')
                    else:
                        f.write(helper.reverse_complement(seq) + '\n')

                elif amp_rate<=0:
                    print('invalide number of amp_rate, please use amp_rate > 0')
                    sys.exit()
                else:
                    # random repeat each UMI x times (x~poi(amp_rate))
                    for j in range(np.random.poisson(amp_rate)):
                        read_index += 1
                        read_id = BC + '_' +  UMI + '_' +  str(read_index) + '_' + str(next_frag[0])
                        f.write(">" + read_id + '\n')
                        if np.random.choice([False,True]):
                            f.write(seq + '\n')
                        else:
                            f.write(helper.reverse_complement(seq) + '\n')
                lock.acquire()
                pbar.update(1)
                lock.release()


def get_random_frag(ref, frag_len):                    
    # There should be one and only one record, the entire genome:
    records = list(SeqIO.parse(ref, "fasta"))
    num_records = len(records)
    limits = [len(r.seq) for r in records]
    while True:
        i = randint(0, num_records-1)
        record = records[i]
        limit = limits[i]
        if frag_len > limit:
            frag = record.seq
            name = record.name.split(' ')[0]
            if 'N' in frag:
                continue
            yield name, helper.reverse_complement(frag)
        else:
            start = randint(0, limit - frag_len)
            end = start + frag_len
            frag = record.seq[start:end]
            name = record.name.split(' ')[0]

            if 'N' in frag:
                continue
            #yield helper.reverse_complement(frag)
            yield name, helper.reverse_complement(frag)


def get_random_frag_with_distribution(ref, frag_len_dist):                    
    # There should be one and only one record, the entire genome:
    records = list(SeqIO.parse(ref, "fasta"))
    num_records = len(records)
    lengths = [len(r.seq) for r in records]
    while True: 
        i = randint(0, num_records-1)
        record = records[i]
        length = lengths[i]
        round_len = helper.round_to_nearest(length+94)
        if frag_len_dist.get(round_len, 0) > 0:
            frag = record.seq
            name = record.name.split(' ')[0]
            if 'N' in frag:
                continue
            frag_len_dist[round_len]-=1
            yield name, helper.reverse_complement(frag)
        else :
            continue


def scale_dist(distribution_dict, repeats):
    total_dist = sum(distribution_dict.values())
    total_umi = repeats.sum()
    scaling_factor = total_umi / total_dist
    scaled_distribution = {length: int(count * scaling_factor) for length, count in distribution_dict.items()}
    scaled_distribution[max(scaled_distribution, key=scaled_distribution.get)] += total_umi - sum(scaled_distribution.values())
    return scaled_distribution


def chunk_generator(df, num_chunk):
    '''split umi count dataframe into chunks with simular UMI count
    '''
    total_umi = df.counts.sum()
    
    umi_per_chunk = total_umi/num_chunk
    current_count = 0
    current_idx = 0
    for idx, i in enumerate(df.counts):
        current_count += i
        if current_count > umi_per_chunk:
            yield df.loc[current_idx: idx]
            current_idx = idx+1
            current_count = 0
    yield df.loc[current_idx:]

def main():
    args = parse_arg()
    BC_df = pd.read_csv(args.umi_count_file)
    #total_umi = BC_df.counts.sum()
    tqdm.set_lock(RLock())
    if args.thread == 1:
        artificial_template(BC_list=BC_df.BC,
                            repeats=BC_df.counts,
                            tran_len=args.mRNA_len,
                            tran_ref=args.trans_ref,
                            adaptor=config.ADPTER_SEQ,
                            TSO=config.TSO_SEQ,
                            output_fn=args.output_filename,
                            fastq_dist=args.fit_distribution,
                            thread=args.thread,
                            amp_rate=args.amp_rate)
    
    else: # multiprocessing
        # create tmp dir
        tmp_dirname = f".tmp_{datetime.now()}"
        os.mkdir(tmp_dirname)
 
        chunks = chunk_generator(BC_df, args.thread)
        executor = concurrent.futures.ProcessPoolExecutor(args.thread)
        futures = {}
        
        for idx, chunk in enumerate(chunks):
            futures[executor.submit(artificial_template,
                            BC_list=chunk.BC,
                            repeats=chunk.counts,
                            tran_len=args.mRNA_len,
                            tran_ref=args.trans_ref,
                            adaptor=config.ADPTER_SEQ,
                            TSO=config.TSO_SEQ,
                            output_fn=f'{tmp_dirname}/tmp_{idx}',
                            fastq_dist=args.fit_distribution,
                            amp_rate=args.amp_rate,
                            pbar_pos = idx+1)] = None
        concurrent.futures.wait(futures)

        # cat the tmp file into one
        with open(args.output_filename, 'wb') as outf:
            for fn in os.listdir(tmp_dirname):
                with open(f'{tmp_dirname}/{fn}', 'rb') as readfile:
                    shutil.copyfileobj(readfile, outf)
        
        
        # clean tmp dir
        print("Cleaning the temp directories...")
        shutil.rmtree(tmp_dirname)
        print("Finished.")
if __name__=='__main__':
    main()
