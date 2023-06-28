#! /usr/bin/env python

version = '1.0.0'

# updating the script to handle hsn values and dropping the duplicates


# import python modules
import os
import glob
import re
import shutil
import pandas as pd 
from datetime import date
import numpy as np

import sys
import argparse



#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--freyja_aggregate_path")
    parser.add_agrument("--workbook_path")
    parser.add_agrument("--")
    parser.add_argument('--seq_run')
    options = parser.parse_args(args)
    return options


def create_list_from_write_lines_input(write_lines_input):
    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list

def prep_aggregrate_df(aggregrate_file_path, seq_run):
    # will need to rewrite this script once we switch to hsn for id
    # inner functions
    def create_rep_sample_id(file_name):
        rep_sample_id = file_name.split('_')[0]
        return rep_sample_id

    def create_sample_id(file_name):
        sample_id = file_name.split('-')[0]
        return sample_id

    def create_site_id(file_name):
        site_id = re.findall('([A-Z]*)', file_name)[0]
        return site_id

    def get_collection_date(file_name):
        Date = re.findall('\D*([0-9]*)', file_name)[0]
        collection_date = '%s-%s-%s' % (Date[0:4], Date[4:6], Date[6:8])
        return collection_date

    df = pd.read_csv(aggregrate_file_path, sep = '\t')
    df = df.rename(columns = {'Unnamed: 0' : 'file_name'})

    # for lineages and abundancees column with Na, replace with None
    df.lineages = df.lineages.fillna('None')
    df.abundances = df.abundances.fillna('None')

     # add some metadata using the functions defined above
    df['rep_sample_id'] = df.apply(lambda x:create_rep_sample_id(x.file_name), axis = 1)
    df['sample_id'] = df.apply(lambda x:create_sample_id(x.file_name), axis = 1)
    df['site_id'] = df.apply(lambda x:create_site_id(x.file_name), axis = 1)
    df['collection_date'] = df.apply(lambda x:get_collection_date(x.file_name), axis = 1)
    df['seq_run'] = seq_run

    # reset column order
    col_order = ['rep_sample_id', 'sample_id', 'site_id', 'collection_date',
            'seq_run', 'coverage', 'lineages', 'abundances']
    df = df[col_order]

    return df

# this function will become obsolete once we swithc to hsn numbers
# metadata1_dict: key = rep_sample_id (i.e sample duplicate)
def create_metadata1_dictionary(df):
    metadata1_dict = {}

    for row in range(df.shape[0]):
        rep_sample_id = df.rep_sample_id[row]
        sample_id = df.sample_id[row]
        site_id = df.site_id[row]
        collection_date = df.collection_date[row]
        seq_run = df.seq_run[row]
        coverage = df.coverage[row]
        metadata1_dict[rep_sample_id] = {'sample_id' : sample_id,
                                        'site_id' : site_id,
                                        'collection_date' : collection_date,
                                        'seq_run' : seq_run,
                                        'coverage' : coverage}
    return metadata1_dict

#metadata2_dict: key = sample_id (ie. sample)
def create_metadata2_dictionary(df):
    metadata2_dict = {}
        
    for row in range(df.shape[0]):
        sample_id = df.sample_id[row]
        if sample_id not in metadata2_dict.keys():
            site_id = df.site_id[row]
            seq_run = df.seq_run[row]
            collection_date = df.collection_date[row]
            metadata2_dict[sample_id] = {'site_id' : site_id,
                                            'collection_date' : collection_date,
                                            'seq_run' : seq_run}

    return metadata2_dict

def create_nested_dictionary_of_lineages_and_abundances(df):
    nested_dict = {}
    for row in range(df.shape[0]):
        rep_sample_id = df.rep_sample_id[row]
        lineages = df.lineages[row].split(' ')
        abundances = df.abundances[row].split(' ')
        lineage_abudnace_dict = dict(zip(lineages, abundances))
        nested_dict[rep_sample_id] = lineage_abudnace_dict

    return nested_dict

def create_long_df(nested_dict, metadata1_dict):
    # initiate dataframe
    columns = ['rep_sample_id', 'sample_id', 'site_id', 'collection_date', 'coverage', 'seq_run', 'lineage', 'abundance']
    long_df = pd.DataFrame(columns = columns)

    # fill in long_df
    index_counter = 0
    for rep_sample_id in nested_dict.keys(): 
        lineage_dict = nested_dict[rep_sample_id]
        for lineage in lineage_dict.keys():
            abundance = lineage_dict[lineage]
            long_df.at[index_counter, 'rep_sample_id'] = rep_sample_id
            long_df.at[index_counter, 'lineage'] = lineage
            long_df.at[index_counter, 'abundance'] = abundance
            index_counter = index_counter + 1


    # # fill in metadata
    for row in range(long_df.shape[0]):
        rep_sample_id = long_df.rep_sample_id[row]
        inner_dict = metadata1_dict[rep_sample_id]
        long_df.at[row, 'sample_id'] = inner_dict['sample_id']
        long_df.at[row, 'site_id'] = inner_dict['site_id']
        long_df.at[row, 'collection_date'] = inner_dict['collection_date']
        long_df.at[row, 'seq_run'] = inner_dict['seq_run']
        long_df.at[row, 'coverage'] = inner_dict['coverage']

    # format na values
    long_df.abundance = long_df.abundance.replace({'LowCov': np.nan, 'None': np.nan})
    long_df.coverage = long_df.coverage.replace({'LowCov': np.nan, 'None': np.nan})
    long_df.coverage = long_df.coverage.astype(float)
    long_df.abundance = long_df.abundance.astype(float)

    return long_df

def average_by_sample_id(long_df, metadata2_dict):
    average_df = long_df.groupby(['sample_id', 'lineage']).mean().reset_index(drop = False)

    #add metadata
    for row in range(average_df.shape[0]):
        sample_id = average_df.sample_id[row]
        inner_dict = metadata2_dict[sample_id]
        average_df.at[row, 'site_id'] = inner_dict['site_id']
        average_df.at[row, 'collection_date'] = inner_dict['collection_date']
        average_df.at[row, 'seq_run'] = inner_dict['seq_run']

    col_order= ['sample_id', 'site_id', 'collection_date', 'coverage', 'seq_run', 'lineage', 'abundance']
    average_df = average_df[col_order]
    return average_df

#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    aggregrate_file_path = options.freyja_aggregate_path
    seq_run = options.seq_run
    workbook_path = options.workbook_path


    # read in workbook
    # generate list of samples
    sample_type_dict = {}
    wb = pd.read_csv(workbook_path, dtyep = {'sample_name' : object, "hsn" : object})
    sample_names = wb.sample_name.unique().tolist()
    for row in range(wb.shape[0]):
        sample_name = wb.sample_name[row]
        sample_type = wb.sample_type[row]
        sample_type_dict[sample_name] = sample_type

    # read in freyja aggregate file
    frejya_df = pd.read_csv(aggregrate_file_path, sep = '\t')
    frejya_df = frejya_df.rename(columns = {'Unnamed: 0' : 'file_name'})
    for row in range(frejya_df):
        abundances = frejya_df.abundances[row]
        if re.search('inf', abundances):
            frejya_df.at[row, 'abundances'] = 'LowCov'
            frejya_df.at[row, 'lineages'] = 'LowCov'
            frejya_df.at[row, 'summarized'] = 'LowCov'
            frejya_df.at[row, 'resid'] = 'LowCov'
            frejya_df.at[row, 'coverage'] = 'LowCov'

    # get sample name from file_name
    for row in range(frejya_df.shape[0]):
        file_name = frejya_df.file_name[row]
        if re.search('_variants.tsv', file_name):
            sample_name = file_name.split('_variants.tsv')[0]
        else:
            sample_name = file_name
        frejya_df.at[row, 'sample_name'] = sample_name


    # fill in empty cells
    frejya_df.lineages = frejya_df.lineages.fillna('None')
    frejya_df.abundances = frejya_df.abundances.fillna('None')

    # create coverage dictionary
    # create nested dictionary of lineages and abudnaces
    coverage_dict = {}
    lineaage_abundance_dict = {}
    replicate_dict = {}
    base_sample_name_dict = {}
    for row in range(frejya_df.shape[0]):
        sample_name = frejya_df.sample_name[row]
        base_sample_name = '-'.join(sample_name.split('-')[:-1])
        base_sample_name_dict[sample_name] = base_sample_name
        replicate = sample_name.split('-')[-1]
        replicate_dict[sample_name] = replicate
        coverage = frejya_df.coverage[row]
        coverage_dict[sample_name] = coverage

        abundances = frejya_df.abundances[row].split(' ')
        lineages = frejya_df.lineages[row].split(' ')

        inner_dict = dict(zip(lineages, abundances))
        lineaage_abundance_dict[sample_name] = inner_dict

    # create long df
    columns = ['sample_name', 'base_sample_name', 'replicate', 'sample_type', 'frejya_coverage', 'lineage', 'abundance']
    long_df = pd.DataFrame(columns = columns)

    # fill in long_df
    index_counter = 0
    for sample_name in sample_names: 
        # add metadata
        sample_type = sample_type_dict[sample_name]
        base_sample_name = base_sample_name[sample_name]
        replicate = replicate_dict[sample_name]
        frejya_coverage = coverage_dict[sample_name]
        lineage_dict = lineaage_abundance_dict[sample_name]

        for lineage in lineage_dict.keys():
            abundance = lineage_dict[lineage]
            long_df.at[index_counter, 'sample_name'] = sample_name
            long_df.at[index_counter, 'base_sample_name'] = base_sample_name
            long_df.at[index_counter, 'replicate'] = replicate
            long_df.at[index_counter, 'sample_type'] = sample_type
            long_df.at[index_counter, 'frejya_coverage'] = frejya_coverage
            long_df.at[index_counter, 'lineage'] = lineage
            long_df.at[index_counter, 'abundance'] = abundance
            index_counter = index_counter + 1

    # now we need to average across replicates
    # rules for averaging
    #1- if both samples >60% coverage average samples
    #2 - drop samples < 60% from average; only use samples >60% in average calcuation
    #3 - if all samples are < 60% then get "insufficient coverage" result
    

    average_df = long_df.groupby(['base_sample_name', 'lineage']).mean().reset_index(drop = False)
    cols = ['sample_name', 'base_sample_name', 'replicate', 'sample_type', 'frejya_coverage', 'lineage', 'abundance']












    df = prep_aggregrate_df(aggregrate_file_path = aggregrate_file_path, seq_run= seq_run)






    # pull metadata into dictionary
    metadata1_dict = create_metadata1_dictionary(df=df)
    metadata2_dict = create_metadata2_dictionary(df = df)
    
    # create nested_dict
    nested_dict = create_nested_dictionary_of_lineages_and_abundances(df = df)

    # create long_dict
    long_df = create_long_df(nested_dict = nested_dict, metadata1_dict= metadata1_dict)

    # save long_df
    outfile1 = 'wwt_variant_abundances_long_format_%s.csv' % seq_run
    long_df.to_csv(outfile1, index = False)

    # average by sample_id
    average_df = average_by_sample_id(long_df = long_df, metadata2_dict=metadata2_dict)

    #save average_df
    outfile2 = 'wwt_variant_abundances_long_format_mean_%s.csv' % seq_run
    average_df.to_csv(outfile2, index = False)
