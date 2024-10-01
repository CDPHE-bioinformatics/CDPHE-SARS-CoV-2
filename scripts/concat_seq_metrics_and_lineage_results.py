#! /usr/bin/env python

## aggregate lineage code modified from:
## https://github.com/CDPHE-bioinformatics/cloud-run-aggregate-lineages/blob/main/src/main.py


import argparse
import sys
import pandas as pd
from datetime import date
import re
import requests
import json

###### URLs ########

CDC_GROUPINGS_URL = 'https://data.cdc.gov/resource/jr58-6ysp.json?$query=SELECT%0A%20%20%60variant%60%2C%0A%20%20count(%60variant%60)%20AS%20%60count_variant%60%2C%0A%20%20min(%60creation_date%60)%20AS%20%60min_creation_date%60%2C%0A%20%20max(%60creation_date%60)%20AS%20%60max_creation_date%60%0AGROUP%20BY%20%60variant%60'
PANGO_ALIAS_KEY_URL = 'https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json'

###### data type dictionaries ######
cov_out_data_types = {'#rname': object,
'startpos' : 'Int64',
'endpos' : 'Int64',
'numreads' : 'Int64',
'covbases' : 'Int64',
'coverage' : 'float64',
'meandepth' : 'float64',
'meanbaseseq' : 'float64',
'meanmapq' : 'float64'}


percent_cov_data_type = {'sample_name': object,
 'aligned_bases': 'Int64',
 'N_bases': 'Int64',
 'non_ambiguous_bases': 'Int64',
 'percent_coverage': 'float64'}


pangolin_data_types = {'taxon': object,
 'lineage': object,
 'conflict': 'float64',
 'ambiguity_score': 'float64',
 'scorpio_call': 'float64',
 'scorpio_support': 'float64',
 'scorpio_conflict': 'float64',
 'scorpio_notes': 'object',
 'version': object,
 'pangolin_version': object,
 'scorpio_version': object,
 'constellation_version': object,
 'is_designated': 'bool',
 'qc_status': object,
 'qc_notes': object,
 'note': object,
 'expanded_lineage': object}

nextclade_clades_data_types = {'fasta_header': object,
 'sample_name': object,
 'hsn' : object,
 'nextclade': object,
 'total_nucleotide_mutations': 'Int64',
 'total_nucleotide_deletions': 'Int64',
 'total_nucleotide_insertions': 'Int64',
 'total_AA_substitutions': 'Int64',
 'total_AA_deletions': 'Int64'}

nextclade_variants_data_types = {'fasta_header': object,
 'sample_name': object,
 'hsn' : object,
 'variant_name': object,
 'gene': object,
 'codon_position': 'Int64',
 'refAA': object,
 'altAA': object,
 'start_nuc_pos': 'Int64',
 'end_nuc_pos': 'Int64'}

terra_data_table_data_types = {'index_position': 'Int64',
 'hsn': object,
 'sample_name': object,
 'sample_well': object,
 'sample_type': object,
 'index_well': object,
 'index_1': object,
 'index_2': object,
 'index_1_id': object,
 'index_2_id': object,
 'index_kit': object,
 'index_set': 'Int64',
 'plate_name': object,
 'project_name': object,
 'run_name': object,
 'library_prep': object,
 'project_type': object,
 'organism': object,
 'run_date': object,
 'read_length': 'Int64',
 'read_type': object,
 'primer_set': object,
 'platform': object,
 'instrument_id': object,
 'tag': object,
 'note': object,
 'verification_set_name': object,
 'fastq_dir': object,
 'workbook_path': object,
 'terra_data_table_path': object,
 'out_dir': object,
 'download_date': object}


#### Lamda Functions ####
def get_sample_name(fasta_header):
    '''
    drops the CO-CDPHE from the fasta header leaving only the sample 
    name as the fasta header.

    Returns
    -------
    string

    '''
    sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
    return sample_name 

def create_fasta_header(sample_name):
    '''
    Adds CO-CDPHE to the sample name to create a fasta header

    Returns
    ------
    string
    '''
    return 'CO-CDPHE-%s' % sample_name

def get_assembly_pass(percent_coverage):
        
    '''
    this function determines if an assembly was generated for a sample.
    Returns True if the percent coverage > 0 and 
    returns False if the percnet coverage = 0.

    Returns
    -------
    Boolean
    '''
    if percent_coverage == 0:
        return False
    if percent_coverage > 0:
        return True      

def get_cdc_grouping(expanded_lineage, cdc_groupings_df):

    '''
    given an expanded lineage, this function determines the cdc aggregate
    lineage

    Parameters
    ----------
    expanded_lineage: string (from pangolin results df)

    cdc_groupings : json data (from get_json_data function)

    pango_alias_key: json data (from get_json_data function)

    Returns
    -------
    string
    '''


    potential_matches = []
    for cdc_expanded_lineage in cdc_groupings_df['expanded_lineage']:
        if (expanded_lineage == cdc_expanded_lineage 
            or expanded_lineage.startswith(cdc_expanded_lineage + '.')
        ):
            potential_matches.append(cdc_expanded_lineage)

    # get most specific CDC grouping 
    # example: if XBB.1.5.2, aggregate to XBB.1.5, not XBB
    try:
        lineage_match = max(potential_matches, key=len)
        aggregated_lineage = cdc_groupings_df.loc[cdc_groupings_df['expanded_lineage'] == lineage_match, 'lineage'].item()
    except ValueError:
        if expanded_lineage == 'Unassigned':
            aggregated_lineage = 'Unassigned'
        else:
            aggregated_lineage = 'Other'
    
    return aggregated_lineage


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument('--sample_name_array')
    parser.add_argument("--cov_out_files",  help= "txt file with list of bam file paths")
    parser.add_argument('--percent_cvg_files', help = 'txt file with list of percent cvg file paths')
    parser.add_argument('--pangolin_lineage_csv', help = 'csv output from pangolin')
    parser.add_argument('--nextclade_clades_csv', help = 'csv output from nextclade parser')
    parser.add_argument('--nextclade_variants_csv')
    parser.add_argument('--nextclade_version')
    parser.add_argument('--project_name')
    parser.add_argument('--terra_data_table_path')
    parser.add_argument('--workflow_version')
    parser.add_argument('--analysis_date')


    options = parser.parse_args(args)
    return options


def create_list_from_write_lines_input(write_lines_input):

    '''
    when the wdl function write_lines() is used, the wdl will write 
    each line in the declared column as a string in a text file. This
    function reads in the text file and converts each line (which is
    typically a file path) to an item in a list. This list can then 
    be used to read in each file in subsequent functions. 
    '''

    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list


def concat_cov_out(cov_out_file_list):

    '''
    cov_out (or coverage_out) files is the {sample_name}_coverage.txt
    file produced by samtools coverage function. From this file we
    pull the number of mapped reads, the mean depth, the average 
    base quality score, and the average mapping score. We loop through
    the cov_out_file_list, read in the file, pull out the metrics and 
    then append everything in a dataframe.

    Parameters
    ----------
    cov_out_file_list: list of file paths to 
        {sample_name}_ccoverage.txt for each sample

    Returns
    -------
    dataframe
    '''

     # initiate dataframe for concatenation
    df = pd.DataFrame()
    sample_name_list = []
    samtools_mapped_reads_list = []
    samtools_depth_list = []
    samtools_baseq_list = []
    samtools_mapq_list = []

    #loop through bam file stats files and pull data
    for file in cov_out_file_list:
        d = pd.read_csv(file, sep = '\t', dtype = cov_out_data_types)
        if re.search('barcode', file):
            # for nanopore runs
            sample_name = re.findall('/([0-9a-zA-Z_\-\.]+)_barcode', file)[0]
        else:
            # for illumina runs
            sample_name = re.findall('/([0-9a-zA-Z_\-\.]+)_coverage.txt', file)[0]

        # pull data from samtools output
        num_reads = d.numreads[0]
        depth = d.meandepth[0]
        baseq = d.meanbaseq[0]
        mapq = d.meanmapq[0]

        sample_name_list.append(sample_name)
        samtools_mapped_reads_list.append(num_reads)
        samtools_depth_list.append(depth)
        samtools_baseq_list.append(baseq)
        samtools_mapq_list.append(mapq)

    df['sample_name'] = sample_name_list
    df['mapped_reads'] = samtools_mapped_reads_list
    df['mean_depth'] = samtools_depth_list
    df['mean_base_quality'] = samtools_baseq_list
    df['mean_map_quality'] = samtools_mapq_list

    return df

def concat_percent_cvg(percent_cvg_file_list):

    '''
    percent_cvg (or percent_coverage) files is the 
    {sample_name}_consensus_cvg_stats.csv file produced by 
    the custom python function calc_percent_coverage.py We loop through
    the percent_cvg_file_list, read in the file, and append the df to
    a new df.

    Parameters
    ----------
    percent_cvg_file_list: list of file paths to 
        {sample_name}_consensus_cvg_stats.csv for each sample

    Returns
    -------
    dataframe
    '''

    df_list = []
    for file in percent_cvg_file_list:
        d = pd.read_csv(file, dtype = percent_cov_data_type)
        df_list.append(d)

    df = pd.concat(df_list)

    return df

def get_df_spike_mutations(variants_csv):

    '''
    This function creates a string of S gene mutations seperated by ";"
    for each sample and then concantenates everything into a dataframe. 

    Parameters
    ----------
    variants_csv: variants csv file produced from nextclade_json_parser.py
        Each row contains a single AA change.

    Returns
    -------
    dataframe

    '''
    
    # read in variants file
    variants = pd.read_csv(variants_csv, dtype = nextclade_variants_data_types)
    variants = variants.drop(columns = 'fasta_header')

    #### filter AA changes to those in the S gene
    crit = variants.gene == 'S'

    spike_variants_df = variants[crit]
    spike_variants_df = spike_variants_df.reset_index(drop = True)

    # generate a df of the sample and their spike variants
    sample_name_list_variants = spike_variants_df.sample_name.unique().tolist()


    df = pd.DataFrame()
    sample_name_list = []
    variant_name_list = []

    seperator = '; '

    for sample_name in sample_name_list_variants:
        sample_name_list.append(sample_name)

        crit = spike_variants_df.sample_name == sample_name
        f = spike_variants_df[crit]
        f = f.reset_index()

        mutations = []
        for row in range(f.shape[0]):
            mutations.append(f.variant_name[row])
        mutations_string = seperator.join(mutations)
        variant_name_list.append(mutations_string)

    df['sample_name'] = sample_name_list
    df['spike_mutations'] = variant_name_list

    return df

def get_json_data(url):
    '''
    pulls json data from a url

    Returns
    --------
    json object

    '''
    response = requests.get(url)
    if not response.ok:
        raise Exception(f'Error downloading JSON data from {url}.')
    data = json.loads(response.text)
    return data

def expand_lineage(aliased_lineage, alias_key):

    '''
    This function takes the short pangoling lineage and using
    the alias key generates the expanded lineage name.
    Used to generate the expanded lineage from the cdc aggregated lineages.

    Returns
    -------
    string: expanded lineage
    '''
    if aliased_lineage == 'Unassigned':
        return 'Unassigned'
    if aliased_lineage == 'Other':
        return 'Other'
    if pd.isna(aliased_lineage):
        return ''
    
    # expand letter part of alias
    alias_split = aliased_lineage.split('.')
    alias_letters = alias_split[0]
    alias_letters_expanded = alias_key[alias_letters]
    # ignore "X" lineages or lineages that don't have an alias
    # examples: XBB (recombinant) and A.1 (doesn't have an alias)
    if alias_letters.startswith('X') or alias_letters_expanded == '':
        return aliased_lineage
    
    # stitch together expanded name
    expanded_lineage = alias_letters_expanded + '.' + '.'.join(alias_split[1:])
    return expanded_lineage


def generate_cdc_groupings_df (cdc_groupings, pango_alias_key):
    
    '''
    this function converts the cdc groupings json data into a dataframe,
    and expands the cdc aggregate linege into it's expanded form 
    using the pango_alias_key json data

    Parameters
    ----------
    cdc_groupings: json data (produced from get_json_data function)

    pangol_alias_key: json data (produced from get_json_data function)
    
    Returns
    -------
    dataframe

    '''
    cdc_groupings_df = pd.DataFrame(cdc_groupings)
    cdc_groupings_df = cdc_groupings_df.rename(columns={'variant': 'lineage'})

    # only include lineages in the CDC datset that exist in the PANGO dataset
    valid_aliases = pango_alias_key.keys()
    cdc_groupings_df = cdc_groupings_df[cdc_groupings_df['lineage'].isin(valid_aliases)]

    cdc_groupings_df['expanded_lineage'] = cdc_groupings_df['lineage'] \
        .apply(expand_lineage, alias_key=pango_alias_key)

    return cdc_groupings_df

def format_pangolin_df(pangolin_lineage_csv, cdc_groupings_df):
    '''
    This function reads in the pangolin lineage csv file, renames headers
    and sets the sample_name to the index to create the pangolin_df. 
    Then this function creates the aggregated_lineage_df from the pangolin_df. 
    It copies the expanded lineage column and uses the get_cc_groupings lambda
    function to determine the aggregated lineage. The aggregated lineage df 
    has only the aggregated lineage column.

    Parameters
    ---------
    pangolin_lineage_csv: pangolin_lineage.csv file 
        (produed from pangolin task in wdl workflow)

    cdc_groupings_df: dataframe of cdc aggregate data 
        (output of generate_cdc_groupings_df function )

    Returns
    --------
    dataframe 1 (pangolin), dataframe 2 (aggregrated lineage)
    '''

    # read in panlogin results
    pangolin = pd.read_csv(pangolin_lineage_csv, dtype = pangolin_data_types)
    pangolin = pangolin.rename(columns = {'lineage': 'pangolin_lineage',
                                          'conflict' : 'pangoLEARN_conflict',
                                          'taxon' : 'fasta_header',
                                          'ambiguity_score' : 'pangolin_ambiguity_score',
                                          'scorpio_call' : 'pangolin_scorpio_call',
                                          'scorpio_support' : 'pangolin_scorpio_support',
                                          'scorpio_conflict' : 'pangolin_scorpio_conflict',
                                          'scorpio_notes' : 'pangolin_scorpio_notes',
                                          'version' : 'pango_designation_version',
                                          'scorpio_version' : 'pangolin_scorpio_version',
                                          'constellation_version' : 'pangolin_constellation_version',
                                          'is_designated' : 'pangolin_is_designated',
                                          'qc_status' : 'pangolin_qc_status',
                                          'qc_notes' : 'pangolin_qc_notes',
                                          'note' : 'pangolin_note'})

    sample_name = pangolin.apply(lambda x:get_sample_name(x.fasta_header), axis = 1)
    pangolin.insert(value = sample_name, column = 'sample_name', loc = 0)
    drop_columns = ['fasta_header', 'pango_designation_version', 'pangolin_scorpio_version',
                    'pangolin_constellation_version']
    pangolin = pangolin.drop(columns = drop_columns)
    pangolin = pangolin.set_index('sample_name')

    # determine  aggregrated lineage by using the expanded lineage column in pangolin df
    aggregated_lineage_df = pangolin[['expanded_lineage']].copy()
    aggregated_lineage_df['aggregated_lineage'] = aggregated_lineage_df.apply(lambda x:get_cdc_grouping(x.expanded_lineage, cdc_groupings_df), axis = 1 )
    aggregated_lineage_df.drop('expanded_lineage', axis=1, inplace=True)

    return pangolin, aggregated_lineage_df

def concat_results(sample_name_list, terra_data_table_path, project_name, 
                    pangolin_df,aggregated_lineage_df,
                    nextclade_clades_csv, nextclade_version,
                   cov_out_df, percent_cvg_df, spike_variants_df, 
                   workflow_version):

    '''
    Joins together all the dataframes to create a single summary output dataframe.
    This output file is written to a csv file called 
    {project_name}_sequencing_results_{workflow_version}.csv

    Returns
    -------
    dataframe: j
    '''

    # create dataframe and fill with constant strings
    df = pd.DataFrame()
    df['sample_name'] = sample_name_list
    df = df.set_index('sample_name')
    df['analysis_date'] = str(date.today())


    # read in terra_data_table
    terra_data_table = pd.read_csv(terra_data_table_path, sep = '\t', dtype = terra_data_table_data_types)
    entity_col = terra_data_table.columns.tolist()[0]

    # Subsitute the sample_name column with sample names in the entity column
    # (first column) of the Terra data table. This is to fix issues caused by
    # renaming the controls when creating the concatenated input data tables.
    # For example, POS becomes POS-cov_2050_grid.
    terra_data_table['sample_name'] = terra_data_table[entity_col]

    terra_data_table = terra_data_table.drop(columns = entity_col)
    terra_data_table = terra_data_table.set_index('sample_name')

    # read in nextclade csv
    nextclade = pd.read_csv(nextclade_clades_csv, dtype = nextclade_clades_data_types)
    nextclade = nextclade.drop(columns = ['fasta_header', 'hsn'])
    nextclade['nextclade_version'] = nextclade_version
    nextclade = nextclade.set_index('sample_name')

    # set index on sample_names to prepare for joining
    cov_out_df = cov_out_df.set_index('sample_name')
    percent_cvg_df = percent_cvg_df.set_index('sample_name')
    spike_variants_df = spike_variants_df.set_index('sample_name')

    # join
    j = df.join(terra_data_table, how = 'left')
    j = j.join(percent_cvg_df, how = 'left')
    j = j.join(cov_out_df, how = 'left')
    j = j.join(nextclade, how = 'left')
    j = j.join(pangolin_df, how = 'left')
    j = j.join(aggregated_lineage_df, how = 'left')
    j = j.join(spike_variants_df, how = 'left')
    j = j.reset_index()

    # add fasta header
    j['fasta_header'] = j.apply(lambda x:create_fasta_header(x.sample_name), axis=1)

    # add assembled column and fill in failed assembles with 0% coveage
    j.percent_coverage = j.percent_coverage.fillna(value = 0)

    j['assembly_pass'] = j.apply(lambda x:get_assembly_pass(x.percent_coverage), axis = 1)

    # order columns
    columns = j.columns.tolist()
    columns.sort()
    primary_columns = ['hsn', 'sample_name', 'project_name', 'plate_name', 
                       'run_name', 'analysis_date', 'run_date', 'assembly_pass', 
                       'percent_coverage', 'nextclade', 'pangolin_lineage', 
                       'expanded_lineage', 'aggregated_lineage', 'spike_mutations']
    for column in columns:
         if column not in primary_columns:
              primary_columns.append(column)
    
    j = j[primary_columns]


    outfile = f'{project_name}_sequencing_results_{workflow_version}.csv' 
    j.to_csv(outfile, index = False)

    return j




if __name__ == '__main__':

    options = getOptions()

    # read in input parameters
    sample_name_array = options.sample_name_array
    terra_data_table_path = options.terra_data_table_path
    cov_out_files = options.cov_out_files
    percent_cvg_files = options.percent_cvg_files
    project_name = options.project_name
    
    pangolin_lineage_csv = options.pangolin_lineage_csv

    nextclade_clades_csv = options.nextclade_clades_csv
    nextclade_variants_csv = options.nextclade_variants_csv
    nextclade_version = options.nextclade_version

    workflow_version = options.workflow_version
    analysis_date = options.analysis_date

    # create lists from the column table txt file input
    sample_name_list = create_list_from_write_lines_input(write_lines_input=sample_name_array)
    cov_out_file_list = create_list_from_write_lines_input(write_lines_input = cov_out_files)
    percent_cvg_file_list = create_list_from_write_lines_input(write_lines_input=percent_cvg_files)
    
    # concat cov_out files and percent_cvg files
    cov_out_df = concat_cov_out(cov_out_file_list=cov_out_file_list)
    percent_cvg_df = concat_percent_cvg(percent_cvg_file_list=percent_cvg_file_list)

    # get df of spike mutations from nextclade file
    spike_variants_df = get_df_spike_mutations(variants_csv = nextclade_variants_csv)

    # set up pangolin dataframe for aggregated lineage
    pango_alias_key = get_json_data(PANGO_ALIAS_KEY_URL)
    cdc_groupings = get_json_data(CDC_GROUPINGS_URL)
    cdc_groupings_df = generate_cdc_groupings_df(cdc_groupings = cdc_groupings,
                                                 pango_alias_key = pango_alias_key)
    pangolin_df, aggregated_lineage_df = format_pangolin_df(pangolin_lineage_csv = pangolin_lineage_csv,
                                  cdc_groupings_df = cdc_groupings_df)

    # create results file
    results_df = concat_results(sample_name_list = sample_name_list,
                                terra_data_table_path = terra_data_table_path,
                                project_name = project_name,
                                pangolin_df=pangolin_df,
                                aggregated_lineage_df = aggregated_lineage_df,
                                nextclade_clades_csv=nextclade_clades_csv,
                                nextclade_version=nextclade_version,
                                cov_out_df=cov_out_df,
                                percent_cvg_df=percent_cvg_df, 
                                spike_variants_df = spike_variants_df,
                                workflow_version = workflow_version)
    

    

