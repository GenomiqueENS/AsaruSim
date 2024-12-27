from pyfaidx import Fasta
import pysam
import unittest
import sys
import os
import pandas as pd

bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(bin_path)
bin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'bin'))
sys.path.append(bin_path)

from bin.intron_retention import get_sequence_with_intron, gtf_to_ref_structure, read_IR_model
from bin.intron_retention import custom_key_function , update_structure, calculate_reference_length
from bin.intron_retention import extract_read_pos

from bin.template_maker import Transcriptome, TemplateGenerator

test_files = {'transcriptome' : 'data/ENST00000383052.fa',
              'genome' : 'data/chromosome_Y.fa.gz',
              'gtf':  'data/ENST00000383052.gtf',
              'model' : 'data/test_IR_markov_model'}

ref_structure = {"ENST00000383052" : [('exon', 'Y', 2935280, 2935446, 166, '+'),
                                            ('exon', 'Y', 2953908, 2953997, 89, '+'),
                                            ('exon', 'Y', 2961073, 2961646, 573, '+'),
                                            ('exon', 'Y', 2975094, 2975244, 150, '+'),
                                            ('exon', 'Y', 2975510, 2975654, 144, '+'),
                                            ('exon', 'Y', 2976669, 2976822, 153, '+'),
                                            ('exon', 'Y', 2977939, 2978080, 141, '+'),
                                            ('exon', 'Y', 2978809, 2982506, 3697, '+'),
                                            ('intron', 'Y', 2935446, 2953908, 18462, '+'),
                                            ('intron', 'Y', 2953997, 2961073, 7076, '+'),
                                            ('intron', 'Y', 2961646, 2975094, 13448, '+'),
                                            ('intron', 'Y', 2975244, 2975510, 266, '+'),
                                            ('intron', 'Y', 2975654, 2976669, 1015, '+'),
                                            ('intron', 'Y', 2976822, 2977939, 1117, '+'),
                                            ('intron', 'Y', 2978080, 2978809, 729, '+')]}

updated_ref_structure=dict()
updated_ref_structure["ENST00000383052"] = [(('retained_intron' if t[0] == 'intron' else t[0]),
                                                  ) + t[1:] for t in ref_structure["ENST00000383052"]]

test_filtered_transcripts = Fasta(test_files['transcriptome'], key_function=custom_key_function)
test_genome_fai = pysam.Fastafile(test_files['genome'])

IR_markov_model = read_IR_model(test_files['model'])

barcodes = ["ACTGCGCTAGCTTTCGACG", "GTACTAGCTATAACGC"]
transcripts = ["ENST00000383052"]
data = [[2, 1]]

matrix = pd.DataFrame(data, index=transcripts, columns=barcodes)

class TestStringMethods(unittest.TestCase):
    def test_ref_structer(self):
        test_ref_structure = gtf_to_ref_structure(test_files['gtf'])
        self.assertEqual(test_ref_structure, ref_structure)

    def test_update_structure(self):
        test_i_r, test_updated_ref_structure = update_structure(ref_structure["ENST00000383052"], IR_markov_model)
        self.assertEqual(test_updated_ref_structure, updated_ref_structure["ENST00000383052"])

    def test_extract_read_pos(self):
        up_ref_str = updated_ref_structure["ENST00000383052"]
        middle_ref = calculate_reference_length(up_ref_str)
        list_iv, _, _ = extract_read_pos(middle_ref, middle_ref, up_ref_str, False)

        new_read = ""
        for interval in list_iv:
            chrom = interval.chrom
            if chrom not in test_genome_fai.references:
                chrom = "chr" + chrom
                if chrom not in test_genome_fai.references:
                    break
            start = interval.start
            end = interval.end
            new_read += test_genome_fai.fetch(chrom, start, end)

        self.assertEqual(len(new_read), 47226)
    
    def test_get_sequence_with_intron(self):
        _, new_cDNA = get_sequence_with_intron("ENST00000383052",
                                            test_genome_fai,
                                            ref_structure, 
                                            IR_markov_model)

        self.assertEqual(len(new_cDNA), 47226)

    def test_TemplateGenerator(self):
        outFasta = "data/out_test.fa"
        len_dT = 10
        transcriptome = Transcriptome(test_files['transcriptome'], False)
        generator = TemplateGenerator(transcriptome = transcriptome, 
                                        adapter = "ATGC", 
                                        len_dT = len_dT, 
                                        TSO = "ATGC",
                                        full_length = True,
                                        length_dist = None,
                                        outFasta = outFasta, 
                                        threads = 1,
                                        amp = 1,
                                        truncation_model=None,
                                        intron_retention=False,
                                        dict_ref_structure=None,
                                        IR_markov_model=None,
                                        genome_fai=None)
        
        for cell in matrix.columns:
                umi_counts = matrix[cell]
                generator.make_template(cell, 
                                        umi_counts, 
                                        None, 
                                        "transcript_id")
        
        with open (outFasta) as outf:
            for line in outf:
                if not line.startswith(">"):
                    len_trx = len(line.strip())
                    break
        os.remove(outFasta) 
        self.assertTrue(abs(len_trx - 5156) <= len_dT, "The length is not within the tolerance")
    
    def test_TemplateGenerator_with_intron_retention(self):
        outFasta_IR = "data/out_test_IR.fa"
        len_dT = 10
        transcriptome = Transcriptome(test_files['transcriptome'], False)
        generator = TemplateGenerator(transcriptome = transcriptome, 
                                        adapter = "ATGC", 
                                        len_dT = len_dT, 
                                        TSO = "ATGC",
                                        full_length = True,
                                        length_dist = None,
                                        outFasta = outFasta_IR, 
                                        threads = 1,
                                        amp = 1,
                                        truncation_model=None,
                                        intron_retention=True,
                                        dict_ref_structure=ref_structure,
                                        IR_markov_model=IR_markov_model,
                                        genome_fai=test_genome_fai)
        
        for cell in matrix.columns:
                umi_counts = matrix[cell]
                generator.make_template(cell, 
                                        umi_counts, 
                                        None, 
                                        "transcript_id")
        
        with open (outFasta_IR) as outf_IR:
            for line in outf_IR:
                if not line.startswith(">"):
                    len_trx = len(line.strip())
                    break
        os.remove(outFasta_IR) 
        self.assertTrue(abs(len_trx - 47269) <= len_dT, "The length is not within the tolerance")

if __name__ == '__main__':
    unittest.main()

