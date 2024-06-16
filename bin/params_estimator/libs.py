import multiprocessing as mp
from math import log
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed


def get_read_lengths(filename, thread = 15, batch_size = 500):
    read_batchs = read_batch_generator(filename, batch_size)
    
    results = multiprocessing_submit(calculate_lengths_batch,
                                      read_batchs, n_process=thread, 
                                      pbar_update=batch_size)

    length_lst = []
    
    for f in results:
        length = f.result()
        length_lst.extend(length)
        
    return length_lst


def calculate_lengths_batch(read_batch):
    lenghts = []
    for seq, quality in read_batch:
        lenghts.append(len(seq))
    return lenghts
    

def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,pbar = True, pbar_update = 500,  *arg, **kwargs):
    executor = ProcessPoolExecutor(n_process)

    max_queue = n_process * 2
    if pbar:
        pbar = tqdm(unit = 'read', desc='Processed')

    futures = {}
    n_job_in_queue = 0
    while True:
        while n_job_in_queue < max_queue:
            i = next(iterator, None)
            if not i:
                break
            futures[executor.submit(func, i, *arg, **kwargs)] = None
            n_job_in_queue += 1

        job = next(as_completed(futures), None)
 
        if job is None:
            break
        else:
            n_job_in_queue -= 1
            pbar.update(pbar_update)
            yield job
            del futures[job]


def read_batch_generator(fastq_fn, batch_size):
    open_fn = gzip.open if fastq_fn.endswith('.gz') else open
    with open(fastq_fn, 'rt') as f:
        batch = []
        line_num = 0
        while True:
            header = f.readline()
            if not header: break
            seq = f.readline().strip().encode('ascii')
            f.readline() 
            quality = f.readline().strip().encode('ascii')
            batch.append([seq, quality])
            if len(batch) == batch_size:
                yield batch
                batch = []
        if len(batch) > 0:
            yield batch