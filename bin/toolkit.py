# Author: Ali Hamraoui
# Date: 25-04-2024
import multiprocessing as mp
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import random
from scipy import stats
import numpy as np
import HTSeq
from badread.simulate import sequence_fragment

ont_3p = [0.46713188751686213, 0.16673847625320887, 0.03520399186107891, 0.01779148300706967, 0.011167330623367366, 0.007589634154594385, 0.0057037865834126626, 0.005360905206834167, 0.004714618927954307, 0.005824471804774172, 0.003999530793905735, 0.003578824368038042, 0.00447776008229153, 0.003863054982833373, 0.004495806470532503, 0.0033385818245800832, 0.003102850878182368, 0.0023415188742663016, 0.002119322719049316, 0.0030261537281582305, 0.002912235902387086, 0.002373100053688005, 0.003551754785676582, 0.007077567888256763, 0.00419578526602632, 0.002210682559519244, 0.004627770684544622, 0.002146392301410776, 0.0036092776481946843, 0.002611086798615842, 0.0025039363684350624, 0.003856287587243008, 0.0038607991843032513, 0.0028445619464834357, 0.003017130534037744, 0.004093146432905784, 0.002603191503760416, 0.002035858173434814, 0.003538219994495851, 0.003646498323941692, 0.004468736888171044, 0.0032054897113029043, 0.0035348362967006688, 0.0019862306057721376, 0.002995700448001588, 0.0024700993904832372, 0.002491529476519393, 0.002064055655061335, 0.003814555314435757, 0.0023054260977843546, 0.004129239209387731, 0.0024678435919531154, 0.0019704400160612855, 0.001982846907976955, 0.002121578517579438, 0.0027633531993990555, 0.00470559573383382, 0.0034919761246283566, 0.0025039363684350624, 0.004561224627906033, 0.0024982968721097582, 0.0020798462447721867, 0.002897573211941295, 0.0023437746727964235, 0.0026392842802423627, 0.0032111292076282085, 0.006500083464545616, 0.006404212027015444, 0.0049221523927255005, 0.0023291119823506323, 0.002362948960302458, 0.003120897266423342, 0.0026697375603990056, 0.0038326017026767308, 0.003440092758435558, 0.008549476429161158, 0.003399488384893368, 0.0025400291449170095, 0.002512959562555549, 0.001722302177747901, 0.002155415495531263, 0.00236633265809764, 0.0031242809642185237, 0.0027306441207122912, 0.0022715891198325294, 0.003955542722568362, 0.0019287077432540345, 0.001819301514543133, 0.003299105350302954, 0.005128557958231634, 0.0014031066857356834, 0.0016433492291936425, 0.0014741643394345165, 0.0014403273614826914, 0.0015813147696152963, 0.0037209396754357077, 0.0004489039074942139, 0.000371078858205016, 0.00019851027065070765, 4.286017207231188e-05, 0.0]
ont_5p = [0.5807485641842356, 0.13738151418220534, 0.05887183003911555, 0.031937595589462714, 0.023225701666132793, 0.014341239155248566, 0.010606764688632128, 0.008649859463751573, 0.00788063216498008, 0.008692719635823885, 0.006399700429955199, 0.004554457232315667, 0.0038506480909177043, 0.004070588447604567, 0.0035111504121343915, 0.0030069794406521965, 0.002865992032519592, 0.0024012975353145264, 0.0024464135059169595, 0.004743944308845889, 0.002703574538350831, 0.0030520954112546305, 0.0021384970065553505, 0.0024520530022422637, 0.0018644174851455667, 0.0016072564527116957, 0.0018802080748564186, 0.002547924439772435, 0.0015249198063622543, 0.0016659072144948593, 0.002093381035952917, 0.0028253876589774012, 0.0018181736152780722, 0.0014234088725067786, 0.002099020532278221, 0.0011854221275789418, 0.0010275162304704242, 0.0011504572503620557, 0.0021892524734830883, 0.0013636302114585542, 0.0013884439952898927, 0.0010320278275306675, 0.0010884227907837095, 0.0010884227907837095, 0.001266630874663322, 0.0011075970782897437, 0.0010252604319403024, 0.0010184930363499375, 0.0008526918443859941, 0.0009000636135185493, 0.0008696103333619067, 0.0008831451245426367, 0.0008989357142534885, 0.0008921683186631235, 0.001091806488578892, 0.0009260052966149487, 0.00098127236060293, 0.0009429237855908614, 0.0011042133804945614, 0.0009654817708920781, 0.0008323896576148991, 0.0007365182200847278, 0.0006756116597714426, 0.0012790377665789913, 0.0007466693134702754, 0.000650797875940104, 0.0006271119913738264, 0.0006756116597714426, 0.0006429025810846782, 0.0010230046334101807, 0.0009891676554583558, 0.0009316447929402529, 0.0006000424090123664, 0.000609065603132853, 0.0005662054310605412, 0.0005797402222412712, 0.0005425195464942635, 0.0005154499641328033, 0.0005256010575183509, 0.0006282398906388872, 0.0008380291539402033, 0.0006598210700605907, 0.0005650775317954804, 0.0004049158361568412, 0.0004274738214580579, 0.0003879973471809286, 0.00049740357589183, 0.0004658223964701265, 0.0004116832317472062, 0.00032370708907246073, 0.0002594168309639929, 0.00023685884566277616, 0.0002436262412531412, 0.00020414976697601188, 0.0001917428750603426, 0.00012970841548199646, 0.00012632471768681395, 7.105765369883285e-05, 4.737176913255523e-05, 6.767395590365033e-06, 0.0]

def adjust_lambda_cut(base_lambda, seq_len):
    return base_lambda / seq_len


def cut_sequence_proba(sequence, base_lambda=1):
        lambda_param = adjust_lambda_cut(base_lambda, len(sequence))
        prob = selection_probas(range(len(sequence)), lambda_param)
        cut_length = np.random.choice(range(len(sequence)), p=prob[::-1])
        return sequence[-cut_length:]


def adjust_lambda(base_lambda, lengths):
    std_dev = np.std(lengths)
    return base_lambda / std_dev


def selection_probas(lengths, lambda_param):
    probabilities = np.exp(-lambda_param * np.array(lengths)**2)
    return probabilities / probabilities.sum()


def logN_probas(transcript_lengths, shape = 0.39146, scale = 830):
    probabilities = [stats.lognorm.pdf(length, shape, scale=scale) for length in transcript_lengths]
    return probabilities/sum(probabilities)


def fetch_transcript_id_by_gene(gene, transcripts_index, base_lambda=0.001):
    transcript_ids = transcripts_index.get(gene, [])
    
    if transcript_ids:
        trans_len_dict = {trns: length for trns, length in transcript_ids}
        
        lengths = list(trans_len_dict.values())
        transcripts = list(trans_len_dict.keys())

        probabilities = logN_probas(lengths)

        if np.any(np.isnan(probabilities)):
            return random.choice(transcripts)
        
        selected_index = np.random.choice(range(len(lengths)), p=probabilities)
        
        return transcripts[selected_index]
        
    else:
        return None


def fetch_transcript_id_by_gene_old(gene, transcripts_index):
    transcript_ids = transcripts_index.get(gene, [])
    if transcript_ids:
        return transcript_ids[0]
    else:
        return None
   

def truncate_cDNA(cDNA, probs_3p=ont_3p, probs_5p=ont_5p):
    cDNA_len = len(cDNA)
    bp_list=list(range(101))

    del3 = int(np.random.choice(bp_list, p=probs_3p) * cDNA_len / 100.0)
    del5 = int(np.random.choice(bp_list, p=probs_5p) * cDNA_len / 100.0)

    new_len = max(min(50, cDNA_len), cDNA_len - del3 - del5)
    start_pos = max(0, min(del5, cDNA_len - new_len))

    new_cDNA = cDNA[start_pos:start_pos+new_len]
    return new_cDNA


def cut_sequence(sequence, length_distribution):
    """
    Function to truncate sequences to a random length based on a given log-normal distribution.
    Inputs:
    sequences: List of sequences (byte strings).
    length_distribution: Tuple (shape, loc, scale) defining the parameters of the log-normal distribution.
    Returns: List of truncated sequences.
    """
    shape, loc, scale = length_distribution
    cut_length = int(stats.lognorm.rvs(shape, loc, scale))-100
    cut_length = min(cut_length, len(sequence))
    return sequence[-cut_length:]


def parse_gtf(gtf_file, index_by, protein_coding=False):

    if index_by == 'gene_id' or index_by == 'gene_name':
        pass
    else:
        logging.error("\033[91mError: Please specify the feature type in the provided matrix (transcript_id, gene_id or gene_name) using '--features' argument.\033[0m")
        sys.exit(1)
        
    transcripts = {}

    gff_features = HTSeq.GFF_Reader(gtf_file, end_included=True)

    for feature in gff_features:
        if index_by in feature.attr:
                index_key = feature.attr[index_by]
        else:
            continue

        if feature.type == "transcript":
            transcript_id = feature.attr.get("transcript_id")
            if index_by in feature.attr:
                index_key = feature.attr[index_by]
            else:
                continue
            
            length = 0
            
            if index_key not in transcripts:
                    transcripts[index_key] = {}
                
            if not protein_coding:
                transcripts[index_key][transcript_id]=length
                
            elif protein_coding and feature.attr['transcript_type'] == 'protein_coding':
                transcripts[index_key][transcript_id]=length
                
        elif feature.type == "exon":
            transcript_id = feature.attr.get("transcript_id")
            if transcript_id in transcripts[index_key]:
                transcripts[index_key][transcript_id] += feature.iv.length
        
    transcripts = {key: list(sub_dict.items()) for key, sub_dict in transcripts.items()}
    
    return transcripts


def parse_gtf_small(gtf_file, index_by, protein_coding=False):

    if index_by == 'gene_id' or index_by == 'gene_name':
        pass
    else:
        logging.error("\033[91mError: Please specify the feature type in the provided matrix (transcript_id, gene_id or gene_name) using '--features' argument.\033[0m")
        sys.exit(1)
        
    transcripts = {}
        
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            columns = line.strip().split("\t")
            if len(columns) > 8:
                Attribute = columns[8].replace('"', '')
                feature = columns[2]
                if feature != 'transcript': 
                    continue
                    
                attributes = dict(item.strip().split(' ',1) for item in Attribute.strip().split(';') if item.strip())
                transcript_id = attributes['transcript_id']
                transcript_type = attributes.get('transcript_biotype', attributes.get('transcript_type', ''))

                index_key = attributes[index_by]
                
                start = int(columns[3]) 
                end = int(columns[4])  
                current_length = end - start 

                if index_key in transcripts:
                    if len(transcripts[index_key]) > 1:
                        _, length = transcripts[index_key]
                        
                        if current_length < length:
                            if not protein_coding:
                                transcripts[index_key] = [transcript_id, length]
                                
                            elif protein_coding and transcript_type == 'protein_coding':
                                transcripts[index_key] = [transcript_id, length]

                else:
                    transcripts[index_key] = [transcript_id, current_length]
   
    transcripts = {key:[transcript_id[0]]  for key, transcript_id in transcripts.items()}
    
    return transcripts


def parse_gtf_old(gtf_file, index_by, protein_coding=True):
    
    transcripts = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            columns = line.strip().split("\t")
            if len(columns) > 8:
                Attribute = columns[8].replace('"', '')
                feature = columns[2]
                if feature != 'transcript': 
                    continue
                    
                attributes = dict(item.strip().split(' ',1) for item in Attribute.strip().split(';') if item.strip())
                transcript_id = attributes['transcript_id']
                gene_id = attributes['gene_id']
                gene_name = attributes['gene_name']
                transcript_type = attributes.get('transcript_biotype', attributes.get('transcript_type', ''))
                
                index_key = None
                if index_by == 'gene_id':
                    index_key = gene_id
                elif index_by == 'gene_name':
                    index_key = gene_name
                
                if index_key:
                    if index_key not in transcripts:
                        transcripts[index_key] = []
                    
                    if not protein_coding:
                        transcripts[index_key].append(transcript_id)
                    elif protein_coding and transcript_type == 'protein_coding':
                        transcripts[index_key].append(transcript_id)

    return transcripts


def custom_key_function(header):
    """Extrait l'identifiant sans la version."""
    return header.split('.')[0]


def generate_random_umi(length=12):
    """Génère une séquence UMI aléatoire de longueur spécifiée."""
    bases = 'ATCG'
    return ''.join(random.choice(bases) for _ in range(length))


def cDNA_pool(trns, transcriptome):
    if trns in transcriptome:
        sequence = transcriptome[trns]
        return sequence[:]


def read_template(fasta, thread=1, batch_size=0):
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


def process_template_chunk(chunk, identities, err_model, q_model):
    results = []
    for name, read in chunk:
        err_seq, quals, _, _ = sequence_fragment(read, identities, err_model, q_model)
        results.append((name, err_seq, quals))
    return results


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