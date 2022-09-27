import argparse
import sys
import pandas as pd
import numpy as np
from datetime import date
import re
from Bio import SeqIO


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument("--sample_id",  help= "")
    parser.add_argument('--fasta_file', help = '')
    options = parser.parse_args(args)
    return options


def calculate_percent_coverage(sample_id, fasta_file):
    # first check that there is only one sequence in the fasta files
    num_records = 0
    with open(fasta_file, 'r') as fasta_handle:
        inside_text = fasta_handle.read()
    for item in re.finditer('>', inside_text):
        num_records = num_records + 1

    # calculate the percent coverage
    if num_records == 1:
        record = SeqIO.read(fasta_file, 'fasta')
        aligned_bases = len(record.seq)
        if aligned_bases == 0:
            Ns = 29903
            num_non_ambigous_bases = 0
            coverage = 0
        else:
            uncorrected_Ns = record.seq.count('N')
            Ns = (29903 - aligned_bases) + uncorrected_Ns
            num_non_ambigous_bases = aligned_bases - uncorrected_Ns
            coverage = round((1-(Ns/29903)) * 100, 2)
    else:
        aligned_bases = np.NaN
        Ns = np.NaN
        num_non_ambigous_bases = np.NaN
        coverage = np.NaN

    # create pd df with calc_percent_cvg
    df = pd.DataFrame()
    df['accession_id'] = [sample_id]
    df['number_aligned_bases'] = [aligned_bases]
    df['number_N_bases'] = [Ns]
    df['number_non_ambigous_bases'] = [num_non_ambigous_bases]
    df['percent_non_ambigous_bases'] = [coverage]
    df['number_seqs_in_fasta'] = [num_records]

    outfile = '%s_consensus_cvg_stats.csv' % sample_id
    df.to_csv(outfile, index = False)
    
    
if __name__ == '__main__':
    
    options = getOptions()
    calculate_percent_coverage(sample_id = options.sample_id, fasta_file = options.fasta_file)