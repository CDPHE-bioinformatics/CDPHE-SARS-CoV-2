#! /usr/bin/env python


import re
import pandas as pd 
import numpy as np
import sys
import json

import argparse

# note before you can use this script from the command line you must make the script executable
# locate this script in ~/scripts and be sure that ~/scripts is ammended to your $PATH variable
# to use nexclade_json_parser.py <path to json file>
# will create two output files:
# nextclade_variant_summary.csv
# nextclade_results.csv

# updated 2022-06-29
# under data['results'][i]['insertions']
# the new nextclade update got ride of the 'lenght' key 
# only keys under insertions are 'pos' and 'ins'
# i use len(data['results'][i]['insertions']['ins'] to find the length


#### FUNCTIONS #####

def hsn_from_id(sample_name):
    # Horizon hsn has format 2XXXXXXXXX
    if re.search(r'2[0-9]{9}', sample_name) is not None:
        return re.search(r'2[0-9]{9}', sample_name).group()
    else:
        return ''

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("--nextclade_json",  help="nextclade json file")
    parser.add_argument('--project_name')
    parser.add_argument('--workflow_version')
    options = parser.parse_args(args)
    return options


def get_sample_name_from_fasta_header(fasta_header):
    if re.search('CO-CDPHE', fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name 
    else:
        return fasta_header

def extract_variant_list(json_path, project_name, workflow_version):

    # create pd data frame to fill
    df = pd.DataFrame()
    sample_name_list = []
    mutation_list = []
    gene_list = []
    refAA_list = []
    altAA_list = []
    codon_pos_list = []
    nuc_start_list = []
    nuc_end_list = []


    with open(json_path) as f:
        data = json.load(f)

    for i in range(len(data['results'])):
        sample_name = data['results'][i]['seqName']
        print(sample_name)

    #     print(data[i]['seqName'])
        if 'aaDeletions' in data['results'][i].keys():
            aa_deletions = data['results'][i]['aaDeletions']
            for item in aa_deletions:
                
                print(sample_name)
                gene=item['gene']
                refAA= item['refAA']
                altAA= 'del'
                pos=item['codon'] + 1
                nuc_start = item['codonNucRange']['begin']
                nuc_end = item['codonNucRange']['end']

                mutation = '%s_%s%d%s' % (gene, refAA, pos, altAA)
                


                sample_name_list.append(data['results'][i]['seqName'])
                mutation_list.append(mutation)
                gene_list.append(gene)
                refAA_list.append(refAA)
                altAA_list.append(altAA)
                codon_pos_list.append(pos)
                nuc_start_list.append(nuc_start)
                nuc_end_list.append(nuc_end)
                
        if 'insertions' in data['results'][i].keys():
            insertions = data['results'][i]['insertions']
            for item in insertions:
                gene= ''
                refAA= 'ins'
                altAA= item['ins']
                pos=item['pos'] + 1
                nuc_start = item['pos'] + 1
                
                # to find length now that update removed length key
                insert_seq = item['ins']
                length = len(insert_seq)
                #####
                
                nuc_end = item['pos'] + 1 + length

                mutation = '%s_%s%d%s' % (gene, refAA, pos, altAA)
                
                sample_name_list.append(data['results'][i]['seqName'])
                mutation_list.append(mutation)
                gene_list.append(gene)
                refAA_list.append(refAA)
                altAA_list.append(altAA)
                codon_pos_list.append(pos)
                nuc_start_list.append(nuc_start)
                nuc_end_list.append(nuc_end)

        if 'aaSubstitutions' in data['results'][i].keys():
            aa_subs = data['results'][i]['aaSubstitutions']
            for item in aa_subs:
                gene = item['gene']
                refAA = item['refAA']
                if item['queryAA'] == '*':
                    altAA = 'stop'
                else:
                    altAA = item['queryAA']     
                
                pos = item['codon'] + 1
                nuc_start = item['codonNucRange']['begin']
                nuc_end = item['codonNucRange']['end']

                mutation = '%s_%s%d%s' % (gene, refAA, pos, altAA)
                sample_name_list.append(data['results'][i]['seqName'])
                mutation_list.append(mutation)
                gene_list.append(gene)
                refAA_list.append(refAA)
                altAA_list.append(altAA)
                codon_pos_list.append(pos)
                nuc_start_list.append(nuc_start)
                nuc_end_list.append(nuc_end)

    print(sample_name_list)
    print('')
    df['fasta_header'] = sample_name_list
    print(df)
    df['sample_name'] = df.apply(lambda x:get_sample_name_from_fasta_header(x.fasta_header), axis = 1)
    print('')
    print(df)
    df['hsn'] = df.apply(lambda x:hsn_from_id(x.sample_name), axis = 1)
    df['variant_name'] = mutation_list
    df['gene'] = gene_list
    df['codon_position'] = codon_pos_list
    df['refAA'] = refAA_list
    df['altAA'] = altAA_list
    df['start_nuc_pos'] = nuc_start_list
    df['end_nuc_pos'] = nuc_end_list

    # save df    
    path = f'{project_name}_nextclade_variant_summary_{workflow_version}.csv' 
    df.to_csv(path, index=False)
    
    
def get_nextclade(json_path, project_name, workflow_version):

    # create pd data frame to fill
    sample_name_list = []
    clade_list = []
    totalSubstitutions_list = []
    totalDeletions_list = []
    totalInsertions_list = []
    totalAASubstitutions_list = []
    totalAADeletions_list = []
    df = pd.DataFrame()

    with open(json_path) as f:
        data = json.load(f)

    for i in range(len(data['results'])):
        if 'clade' in data['results'][i].keys():
            sample_name_list.append(data['results'][i]['seqName'])
            clade_list.append(data['results'][i]['clade'])
            totalSubstitutions_list.append(data['results'][i]['totalSubstitutions'])
            totalDeletions_list.append(data['results'][i]['totalDeletions'])
            totalInsertions_list.append(data['results'][i]['totalInsertions'])
            totalAASubstitutions_list.append(data['results'][i]['totalAminoacidSubstitutions'])
            totalAADeletions_list.append(data['results'][i]['totalAminoacidDeletions'])
            

    df['fasta_header'] = sample_name_list
    df['sample_name'] = df.apply(lambda x:get_sample_name_from_fasta_header(x.fasta_header), axis = 1)
    df['hsn'] = df.apply(lambda x:hsn_from_id(x.sample_name), axis = 1)
    df['nextclade'] = clade_list
    df['total_nucleotide_mutations'] = totalSubstitutions_list
    df['total_nucleotide_deletions'] = totalDeletions_list
    df['total_nucleotide_insertions'] = totalInsertions_list
    df['total_AA_substitutions'] = totalAASubstitutions_list
    df['total_AA_deletions'] = totalAADeletions_list
    
    # save df to file
    path = f'{project_name}_nextclade_results_{workflow_version}.csv' 
    df.to_csv(path, index = False)

    
if __name__ == '__main__':
    
    options = getOptions()
    nextclade_json = options.nextclade_json
    project_name = options.project_name
    workflow_version = options.workflow_version

    get_nextclade(json_path = nextclade_json, project_name = project_name, workflow_version = workflow_version)
    extract_variant_list(json_path = nextclade_json, project_name = project_name, workflow_version  = workflow_version)
        