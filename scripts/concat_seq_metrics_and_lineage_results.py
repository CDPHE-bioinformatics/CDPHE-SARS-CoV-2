#! /usr/bin/env python

import argparse
import sys
import pandas as pd
from datetime import date
import re

# data type dictionaries
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
 'pangolin_version': 'float64',
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




#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument('--sample_name_array')
    parser.add_argument("--cov_out_files",  help= "txt file with list of bam file paths")
    parser.add_argument('--percent_cvg_files', help = 'txt file with list of percent cvg file paths')
    parser.add_argument('--pangolin_lineage_csv', help = 'csv output from pangolin')
    parser.add_argument('--cdc_lineage_groups_json', help = 'json file containing lineage groups for aggregating lineages')
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
    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list


def concat_cov_out(cov_out_file_list):

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

    df_list = []
    for file in percent_cvg_file_list:
        d = pd.read_csv(file, dtype = percent_cov_data_type)
        df_list.append(d)

    df = pd.concat(df_list)

    return df

def get_df_spike_mutations(variants_csv):

    def get_sample_name(fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name 
    
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

def concat_results(sample_name_list, terra_data_table_path, project_name, 
                    pangolin_lineage_csv, cdc_lineage_groups_json,
                    nextclade_clades_csv, nextclade_version,
                   cov_out_df, percent_cvg_df, spike_variants_df, 
                   workflow_version):

    # set some functions for getting data formatted
    def get_sample_name_from_fasta_header(fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name

    def create_fasta_header(sample_name):
        return 'CO-CDPHE-%s' % sample_name
    

    # create dataframe and fill with constant strings
    df = pd.DataFrame()
    df['sample_name'] = sample_name_list
    df = df.set_index('sample_name')
    df['analysis_date'] = str(date.today())


    # read in terra_data_table
    terra_data_table = pd.read_csv(terra_data_table_path, sep = '\t', dtype = terra_data_table_data_types)
    drop_col = terra_data_table.columns.tolist()[0]
    terra_data_table = terra_data_table.drop(columns = drop_col)
    terra_data_table = terra_data_table.set_index('sample_name')

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

    sample_name = pangolin.apply(lambda x:get_sample_name_from_fasta_header(x.fasta_header), axis = 1)
    pangolin.insert(value = sample_name, column = 'sample_name', loc = 0)
    drop_columns = ['fasta_header', 'pango_designation_version', 'pangolin_scorpio_version'
                    'pangolin_constellation_version']
    pangolin = pangolin.drop(columns = drop_columns)
    pangolin = pangolin.set_index('sample_name')

    # read in list of CDC lineage groups
    cdc_lineage_groups_df = pd.read_json(cdc_lineage_groups_json)
    aggregated_lineage_df = aggregate_lineage(pangolin, cdc_lineage_groups_df)


    # read in nextclade csv
    nextclade = pd.read_csv(nextclade_clades_csv, dtype = nextclade_clades_data_types)
    nextclade = nextclade.drop(columns = ['fasta_header', 'hsn', 'nextclade_version'])
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
    j = j.join(pangolin, how = 'left')
    j = j.join(aggregated_lineage_df, how = 'left')
    j = j.join(spike_variants_df, how = 'left')
    j = j.reset_index()

    # add fasta header
    j['fasta_header'] = j.apply(lambda x:create_fasta_header(x.sample_name), axis=1)

    # add assembled column and fill in failed assembles with 0% coveage
    j.percent_coverage = j.percent_coverage.fillna(value = 0)

    def get_assembly_pass(percent_coverage):
        if percent_coverage == 0:
            return False
        if percent_coverage > 0:
            return True      
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


    outfile = f'{project_name}_sequencing_results_v{workflow_version}.csv' 
    j.to_csv(outfile, index = False)

    return j

def aggregate_lineage(pangolin_df, cdc_lineage_groups_df):
    aggregated_lineage_df = pangolin_df[['expanded_lineage']].copy()

    def get_cdc_grouping(expanded_lineage, cdc_lineage_groups_df):
        potential_matches = []
        for cdc_expanded_lineage in cdc_lineage_groups_df['expanded_lineage']:
            if (expanded_lineage == cdc_expanded_lineage 
                or expanded_lineage.startswith(cdc_expanded_lineage + '.')
            ):
                potential_matches.append(cdc_expanded_lineage)

        # get most specific CDC grouping 
        # example: if XBB.1.5.2, aggregate to XBB.1.5, not XBB
        try:
            match = max(potential_matches, key=len)
            aggregated_lineage = cdc_lineage_groups_df.loc[cdc_lineage_groups_df['expanded_lineage'] == match, 'pangolin_lineage'].item()
        except ValueError:
            if expanded_lineage == 'Unassigned':
                aggregated_lineage = 'Unassigned'
            else:
                aggregated_lineage = 'Other'
        
        return aggregated_lineage


    aggregated_lineage_df['aggregated_lineage'] = aggregated_lineage_df.apply(lambda x: get_cdc_grouping(x['expanded_lineage'],
                                                                                                         cdc_lineage_groups_df), axis=1)
    aggregated_lineage_df.drop('expanded_lineage', axis=1, inplace=True)
    
    return aggregated_lineage_df

def make_wgs_horizon_output (results_df, project_name, pangolin_version, 
                             analysis_date, workflow_version):

    results_df['report_to_epi'] = ''
    results_df['Run_Date'] = str(date.today())

    # rename columns 
    results_df = results_df.rename(columns = {'hsn' : 'accession_id', 'pango_designation_version' : 'pangoLEARN_version'})
    
    results_df['pangolin_version'] = pangolin_version
    results_df['workflow_version'] = workflow_version
    results_df['anlaysis_date'] = analysis_date

    col_order = ['accession_id', 'percent_coverage', 'pangolin_lineage', 'pangolin_version',
                 'report_to_epi', 'Run_Date', 'pangoLEARN_version']
    
   
    results_df = results_df[col_order]

    outfile = "%s_wgs_horizon_report.csv" % project_name
    results_df.to_csv(outfile, index = False)


    # for new horizon parser when ready...
    # results_df['pangolin_version'] = pangolin_version
    # results_df['workflow_version'] = workflow_version
    # results_df['anlaysis_date'] = analysis_date

    # col_order2 = ['hsn', 'percent_coverage', 'mean_depth', 'pangolin_lineage',
    #               'aggregated_lineage', 'expanded_lineage', 'pangolin_version',
    #               'project_name', 'platform', 'workflow_version', 'anlaysis_date', 'run_date']

    # results_df = results_df[col_order2]

    # outfile = "%s_wgs_horizon_report.csv" % project_name
    # results_df.to_csv(outfile, index = False)


if __name__ == '__main__':

    options = getOptions()

    sample_name_array = options.sample_name_array
    terra_data_table_path = options.terra_data_table_path
    cov_out_files = options.cov_out_files
    percent_cvg_files = options.percent_cvg_files
    project_name = options.project_name

    pangolin_lineage_csv = options.pangolin_lineage_csv
    cdc_lineage_groups_json = options.cdc_lineage_groups_json

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

    # create results file
    results_df = concat_results(sample_name_list = sample_name_list,
                                terra_data_table_path = terra_data_table_path,
                                project_name = project_name,
                                pangolin_lineage_csv=pangolin_lineage_csv,
                                cdc_lineage_groups_json=cdc_lineage_groups_json,
                                nextclade_clades_csv=nextclade_clades_csv,
                                nextclade_version=nextclade_version,
                                cov_out_df=cov_out_df,
                                percent_cvg_df=percent_cvg_df, 
                                spike_variants_df = spike_variants_df,
                                workflow_version = workflow_version)
    
    # create wgs horizon output
    make_wgs_horizon_output(project_name = project_name,
                            results_df=results_df,
                            workflow_version = workflow_version,
                            analysis_date = analysis_date)
    

