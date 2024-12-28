import HTSeq
import random
import copy
from six.moves import xrange


def gtf_to_ref_structure(gff_file):
    dict_ref_structure = {}
    gff_features = HTSeq.GFF_Reader(gff_file, end_included=True)
    for feature in gff_features:
        if feature.type == "exon" or feature.type == "intron":
            if "transcript_id" in feature.attr:
                feature_id = feature.attr['transcript_id']
            elif "Parent" in feature.attr:
                info = feature.name.split(":")
                if len(info) == 1:
                    feature_id = info[0]
                else:
                    if info[0] == "transcript":
                        feature_id = info[1]
                    else:
                        continue
            else:
                continue

            feature_id = feature_id.split(".")[0]
            if feature_id not in dict_ref_structure:
                dict_ref_structure[feature_id] = []

            # remove "chr" from chromosome names to be consistent
            if "chr" in feature.iv.chrom:
                feature.iv.chrom = feature.iv.chrom.strip("chr")

            dict_ref_structure[feature_id].append((feature.type, feature.iv.chrom, feature.iv.start,
                                                   feature.iv.end, feature.iv.length, feature.iv.strand))
    return dict_ref_structure


def read_IR_model(model_file):
    IR_markov_model = {}
    with open(model_file, "r") as IR_markov:
        IR_markov.readline()
        for line in IR_markov:
            info = line.strip().split()
            k = info[0]
            IR_markov_model[k] = {}
            IR_markov_model[k][(0, float(info[1]))] = "no_IR"
            IR_markov_model[k][(float(info[1]), float(info[1]) + float(info[2]))] = "IR"
    return IR_markov_model


def update_structure(ref_trx_structure, IR_markov_model, seed=123):
    random.seed(seed)
    count = 0
    for item in ref_trx_structure:
        if item[0] == "intron":
            count += 1

    list_states = []
    flag_ir = False
    prev_state = "start"
    for i in range(0, count):
        p = random.random()
        for key in IR_markov_model[prev_state]:
            if key[0] <= p < key[1]:
                flag = IR_markov_model[prev_state][key]
                if flag == "IR":
                    flag_ir = True
                list_states.append(flag)
                prev_state = flag
                break

    if flag_ir:
        ref_trx_structure_temp = copy.deepcopy(ref_trx_structure)
        j = -1
        for i in xrange(0, len(ref_trx_structure_temp)):
            if ref_trx_structure_temp[i][0] == "intron":
                j += 1
                if list_states[j] == "IR":
                    ref_trx_structure_temp[i] = ("retained_intron",) + ref_trx_structure_temp[i][1:]
    else:
        ref_trx_structure_temp = ref_trx_structure

    return flag_ir, ref_trx_structure_temp


def extract_read_pos(length, ref_len, ref_trx_structure, polya, buffer=10):

    len_before = 0
    for item in ref_trx_structure:
        if item[0] == "exon":
            len_before += item[4]
        elif item[0] == "retained_intron":
            break

    #TODO change the random into something truer
    start_pos = random.randint(0, min(ref_len - length, len_before))  # make sure the retained_intron is included

    list_intervals = []
    ir_list = []
    for item in ref_trx_structure:
        if length == 0:
            break
        chrom = item[1]
        if item[0] in ["exon", "retained_intron"]:
            if start_pos < item[4]:
                start = start_pos + item[2]
                if start + length <= item[3]:
                    end = start + length
                else:
                    end = item[3]
                length -= end - start
                start_pos = 0
                iv = HTSeq.GenomicInterval(chrom, start, end, item[5])
                list_intervals.append(iv)
                if item[0] == "retained_intron":
                    ir_list.append((start, end))
            else:
                start_pos -= item[4]

    if polya and end + buffer >= ref_trx_structure[-1][3]:
        retain_polya = True
    else:
        retain_polya = False

    return list_intervals, retain_polya, ir_list


def calculate_reference_length(ref_trx_structure):
    total_length = 0
    
    for item in ref_trx_structure:
        if item[0] in ["exon", "retained_intron"]:
            total_length += item[4]

    return total_length


def custom_key_function(header):
    return header.split('.')[0]


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq


def fitch_trx_by_tructure(structure, genome_fai):
    new_cDNA = ""

    unspliced_trx_len = calculate_reference_length(structure)
    list_iv, _, _ = extract_read_pos(unspliced_trx_len, unspliced_trx_len, structure, False)

    for interval in list_iv:
        chrom = interval.chrom
        if chrom not in genome_fai.references:
            chrom = "chr" + chrom
            if chrom not in genome_fai.references:
                return None
        start = interval.start
        end = interval.end
        new_cDNA += genome_fai.fetch(chrom, start, end) 

        if interval.strand == '-':
            new_cDNA = reverse_complement(new_cDNA)
        
    return new_cDNA


def get_unspliced_sequence(trx, genome_fai, dict_ref_structure):

    if trx in dict_ref_structure:
        ref_structure = dict_ref_structure[trx]
    else:
        return None

    unspliced_structure = [
    ('retained_intron', *entry[1:]) if entry[0] == 'intron' else entry
    for entry in ref_structure]

    unspliced_structure = sorted(unspliced_structure, key=lambda x: x[2])

    unspliced_trx = fitch_trx_by_tructure(unspliced_structure, genome_fai)

    return unspliced_trx


def get_sequence_with_intron(trx, genome_fai, dict_ref_structure, IR_markov_model):

    new_cDNA = ""

    if trx in dict_ref_structure:
        ref_structure = dict_ref_structure[trx]
    else:
        return False, None

    retained, ref_trx_structure_new = update_structure(ref_structure, IR_markov_model)

    if retained:
        ref_trx_structure_new = sorted(ref_trx_structure_new, key=lambda x: x[2])

        new_cDNA = fitch_trx_by_tructure(ref_trx_structure_new, genome_fai)

        # ref_trx_len = calculate_reference_length(ref_trx_structure_new)

        # list_iv, _, _ = extract_read_pos(ref_trx_len, ref_trx_len, ref_trx_structure_new, False)

        # for interval in list_iv:
        #     chrom = interval.chrom
        #     if chrom not in genome_fai.references:
        #         chrom = "chr" + chrom
        #         if chrom not in genome_fai.references:
        #             break
        #     start = interval.start
        #     end = interval.end
        #     new_cDNA += genome_fai.fetch(chrom, start, end) 

        #     if interval.strand == '-':
        #         new_cDNA = reverse_complement(new_cDNA)
    else: new_cDNA = None

    return retained, new_cDNA


if __name__ == "__main__":
    get_sequence_with_intron()