import multiprocessing as mp
import gzip

from math import log
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

from QC_tools.fqread.fastq_reader import SequenceStats


def process_fastq_file(filename, pos, reverse):
    stats = SequenceStats(pos)
    sequence_count = 0 
    
    open_fn = gzip.open if filename[0].endswith('.gz') else open
    with open_fn(filename[0], 'rb') as f:
        while sequence_count < 100000:
            header = f.readline()
            if not header: break
            if reverse:
                seq = f.readline().strip()[::-1][:pos]
                f.readline() 
                quality = f.readline().strip()[::-1]
            else:
                seq = f.readline().strip()[:pos]
                f.readline()
                quality = f.readline().strip()
            
            stats.add_sequence(seq, quality)
            sequence_count += 1

    return stats.calculate_stats()

