"""
This script is based on the original work from https://github.com/youyupei/SLSim/tree/main.
It has been modified and adapted
to meet specific requirements and functionalities not present in the original version.
All original rights are acknowledged.
"""
import numpy as np
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
from tqdm import tqdm
import os
import sys

def reverse_complement(seq):
	'''
	Args: <str>
		queried seq
	Returns: <str>
		reverse_complement seq
	'''
	comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
					'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	letters = \
		[comp[base] if base in comp.keys() else base for base in seq]
	return ''.join(letters)[::-1]

def err_msg(msg):
	CRED = '\033[91m'
	CEND = '\033[0m'
	print(CRED + msg + CEND)	

def warning_msg(msg):
	CRED = '\033[93m'
	CEND = '\033[0m'
	print(CRED + msg + CEND)

def green_msg(msg):
    CRED = '\033[92m'
    CEND = '\033[0m'
    print(CRED + msg + CEND)

def sliding_window_sum(array, window) :
    cum = np.cumsum(array)  
    return cum[window:] - cum[:-window]

def sliding_window_mean(array, window) :
    cum = np.cumsum(array)  
    return (cum[window:] - cum[:-window]) / window


class param:
    def __init__(self, **kwargs):
        for key in kwargs.keys():
            self.__dict__[key] = kwargs[key]
    def add(self, attr_name, attr_val, overwrite = True):
        if attr_name in __dict__.keys() and not overwrite:
            pass
        else:    
            self.__dict__[attr_name] = attr_val
    def rm(self, attr_name):
        if attr_name in self.__dict__.keys():
            del self.__dict__[attr_name]
    def __str__(self):
        return str(self.__dict__)
    
    def check(self, attr_list, add_none = True, silent = True):
        """
        Check whether the attributes in a given list are present. 

        Parameters
        ----------
        attr_list : LIST
            list of strings of the attributes to check 
        add_none : BOOL
            whether or not to create the missed attributes with value None
        silent : BOOL
            always return True if silent = True
        Returns
        -------
        True/False if all attributes is present
        """
        try:
            assert isinstance(add_none, bool)
            assert isinstance(silent, bool)
        except (AttributeError, TypeError):
            raise AssertionError("add_none and silent must be bool variable.")
        
        check_res = True
        for attr in attr_list:
            if attr not in self.__dict__.keys():
                check_res = False
                self.__dict__[attr] = None
        return check_res if not silent else True
    
# multiprocessing
def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,pbar = True, pbar_update = 500,  *arg, **kwargs):
    executor = concurrent.futures.ProcessPoolExecutor(n_process)
    
    # A dictionary which will contain the  future object
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

        # will wait until as least one job finished
        job = next(as_completed(futures), None)
        
        # no more job  
        if job is None:
            break
        # otherwise
        else:
            n_job_in_queue -= 1
            # update pregress bar based on batch size
            pbar.update(pbar_update)
            yield job
            del futures[job]


# check file exist
def check_exist(file_list):
    exit_code = 0
    for fn in file_list:
        if not os.path.exists(fn):
            exit_code = 1
            err_msg(f'Error: can not find {fn}')
    if exit_code == 1:
        sys.exit()

# get file with a certian extensions
def get_files(search_dir, extensions, recursive=True):
    files = []
    if recursive:
        for i in extensions:
            files.extend(Path(search_dir).rglob(i))
        return files
    else:
        for i in extensions:
            files.extend(Path(search_dir).glob(i))
        return files
    

# split any iterator in to batches  
def batch_iterator(iterator, batch_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):
        yield batch


def round_to_nearest(value, base=10):
    return base * round(value / base)


def fastq_batch_generator(fastq, batch_size):
    open_fn = gzip.open if fastq.endswith('.gz') else open

    with open_fn(fastq, 'rt') as fq:
        batch = []
        line_num = 0
        for line in fq:
            line_num += 1
            if line_num % 4 == 0:
                batch.append(line.rstrip())
                if len(batch) == batch_size:
                    yield batch
                    batch = []

        if len(batch) > 0:
            yield batch


def length_read_by_batch(read_batch):
    lengths = [len(r) for r in read_batch]
    return lengths


def get_distribution(fastq, thread, step=10, batch_size = 200):
    thread = max(thread,4)
    read_batchs = fastq_batch_generator(fastq, batch_size)
    rst_futures = multiprocessing_submit(length_read_by_batch,
                    read_batchs, n_process=thread, 
                    pbar_update=batch_size)
    
    reads_len_dist = []
    for idx, f in enumerate(rst_futures):
        reads_len_dist.extend(f.result())

    distribution_dict = {}
    for value in reads_len_dist:
        key = (value // step) * step 
        distribution_dict[key] = distribution_dict.get(key, 0) + 1
    return distribution_dict