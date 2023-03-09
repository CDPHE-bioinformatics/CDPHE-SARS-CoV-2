#! /usr/bin/env python

# updated 2022-03-04
## added assembler version

#updated 2022-04-04
## adjust columns to account for pangolin v4.0 major update

# updated 2022-07-07
## gets rid of the report_to_epi stuff in the wgs_horizon_csv (the column will just be blank).

# update 2022-12-29
## adds in the expanded pangolin lineage into the results output (no big changes)

# update 2023-03-01
## major update; refactoring to match universal naming conventions of the workbook generator;
## also reads in the workbook path to pull various values

import argparse
import sys
import pandas as pd
from datetime import date
import re

#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument('--sample_name_array')
    parser.add_argument('--workbook_path')
    parser.add_argument("--cov_out_files",  help= "txt file with list of bam file paths")
    parser.add_argument('--percent_cvg_files', help = 'txt file with list of percent cvg file paths')
    parser.add_argument('--assembler_version')
    parser.add_argument('--pangolin_lineage_csv', help = 'csv output from pangolin')
    parser.add_argument('--nextclade_clades_csv', help = 'csv output from nextclade parser')
    parser.add_argument('--nextclade_variants_csv')
    parser.add_argument('--nextclade_version')
    parser.add_argument('--project_name')

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
        d = pd.read_csv(file, sep = '\t')
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
        d = pd.read_csv(file, dtype = {'accession_id' : object})
        d = d.rename(columns = {'accession_id' : 'sample_name'})
        df_list.append(d)

    df = pd.concat(df_list)

    return df

def get_df_spike_mutations(variants_csv):

    def get_sample_name(fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name 
    
    # read in variants file
    variants = pd.read_csv(variants_csv, dtype = {'sample_name' : object})
    variants = variants.rename(columns = {'sample_name' : 'fasta_header'})

    variants['sample_name'] = variants.apply(lambda x:get_sample_name(x.fasta_header), axis = 1)
    variants = variants.drop(columns = 'fasta_header')
    # print(variants)

    #### filter variants for spike protein varaints in rbd and pbcs #####
    crit = variants.gene == 'S'
    critRBD = (variants.codon_position >= 461) & (variants.codon_position <= 509)
    critPBCS = (variants.codon_position >= 677) & (variants.codon_position <= 694)
    crit732 = variants.codon_position == 732
    crit452 = variants.codon_position == 452
    crit253 = variants.codon_position == 253
    crit13 = variants.codon_position == 13
    crit145 = variants.codon_position == 145 # delta plus AY.4.2
    crit222 = variants.codon_position == 222 # delta plus AY.4.2
    critdel = variants.variant_name.str.contains('del') # mainly for 69/70 del

    spike_variants_df = variants[crit & (critRBD | critPBCS | crit732 | critdel | crit452 | crit253 | crit13 | crit145 | crit222)]
    spike_variants_df = spike_variants_df.reset_index(drop = True)
    # print(spike_variants_df)

    # generate a df of the sample and their spike variants
    sample_name_list_variants = spike_variants_df.sample_name.unique().tolist()
    # print(sample_name_list)

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
    print(df)

    return df

def concat_results(sample_name_list, workbook_path, project_name, 
                   assembler_version, pangolin_lineage_csv,
                    nextclade_clades_csv, nextclade_version,
                   cov_out_df, percent_cvg_df, spike_variants_df):

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
    df['assembler_version'] = assembler_version
    # print(df)
    # read in workbook
    workbook = pd.read_csv(workbook_path, sep = '\t')
    workbook = workbook.set_index('sample_name')

    # read in panlogin results
    pangolin = pd.read_csv(pangolin_lineage_csv, dtype = {'taxon' : object})
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
    pangolin = pangolin.drop(columns = 'fasta_header')
    pangolin = pangolin.set_index('sample_name')


    # read in nextclade csv
    nextclade = pd.read_csv(nextclade_clades_csv, dtype = {'accession_id' : object})
    nextclade = nextclade.rename(columns = {'sample_name' : 'fasta_header'})
    sample_name = nextclade.apply(lambda x:get_sample_name_from_fasta_header(x.fasta_header), axis = 1)
    nextclade.insert(value = sample_name, column = 'sample_name', loc = 0)
    nextclade = nextclade.drop(columns = 'fasta_header')
    nextclade['nextclade_version'] = nextclade_version
    nextclade = nextclade.set_index('sample_name')


    # set index on the samtools_df and percent_cvg_df and variants_df to prepare for joining
    cov_out_df = cov_out_df.set_index('sample_name')
    percent_cvg_df = percent_cvg_df.set_index('sample_name')
    spike_variants_df = spike_variants_df.set_index('sample_name')
    print(spike_variants_df)

    # join
    j = df.join(workbook, how = 'left')
    j = j.join(percent_cvg_df, how = 'left')
    j = j.join(cov_out_df, how = 'left')
    j = j.join(nextclade, how = 'left')
    j = j.join(pangolin, how = 'left')
    j = j.join(spike_variants_df, how = 'left')
    j = j.reset_index()
    # print(j)
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
                 'expanded_lineage', 'spike_mutations']
    for column in columns:
         if column not in primary_columns:
              primary_columns.append(column)
    
    j = j[primary_columns]


#     # add in 'failed assembly" in missing columns
#     j.spike_mutations = j.spike_mutations.fillna(value = '')
#     j.nextclade = j.nextclade.fillna(value = '')
#     j.nextclade_version = j.nextclade_version.fillna(value = next_version)
    # j.pangolin_lineage = j.pangolin_lineage.fillna(value = 'not assembled')
    # j.expanded_lineage = j.expanded_lineage.fillna(value = 'not assembled')
    # j.percent_non_ambigous_bases = j.percent_non_ambigous_bases.fillna(value = 0)
#     j.mean_depth = j.mean_depth.fillna(value = 0)
#     j.number_aligned_bases = j.number_aligned_bases.fillna(value = 0)
#     j.number_seqs_in_fasta = j.number_seqs_in_fasta.fillna(value = 0)
#     j.num_reads = j.num_reads.fillna(value = 0)
#     j.mean_base_quality = j.mean_base_quality.fillna(value = 0)
#     j.mean_map_quality = j.mean_map_quality.fillna(value = 0)
#     j.number_N_bases = j.number_N_bases.fillna(value = 29903)
#     j.pangolin_version = j.pangolin_version.fillna(value = pangolin_version )
# #     j.pangolin_designation_version = j.pangolin_designation_version.fillna(value = pango_learn_version )
#     j.assembler_version = j.assembler_version.fillna(value = assembler_version_txt)

    outfile = '%s_sequencing_results.csv' % project_name
    j.to_csv(outfile, index = False)

    return j

# def get_project_name(workbook_path):
#     df = pd.read_csv(workbook_path, sep = '\t')
#     project_name = df.project_name[0]

#     return project_name


def make_wgs_horizon_output (results_df, project_name):

    results_df['report_to_epi'] = ''
    results_df['Run_Date'] = str(date.today())

    # rename columns 
    results_df = results_df.rename(columns = {'hsn' : 'accession_id', 'pango_designation_version' : 'pangoLEARN_version'})

    # # add pangolin to the as prefix to version
    # # pull out panglin version 
    
    # for row in range(results_df.shape[0]):
    #     version = results_df.pangolin_version[row]
    #     new_version_txt = f'pangolin {version}'
    #     results_df.at[row, 'pangolin_version'] = new_version_txt

    col_order = ['accession_id', 'percent_coverage', 'pangolin_lineage', 'pangolin_version',
                 'report_to_epi', 'Run_Date', 'pangoLEARN_version']
    
    results_df = results_df[col_order]

    outfile = "%s_wgs_horizon_report.csv" % project_name
    results_df.to_csv(outfile, index = False)


if __name__ == '__main__':

    options = getOptions()

    sample_name_array = options.sample_name_array
    workbook_path = options.workbook_path
    cov_out_files = options.cov_out_files
    percent_cvg_files = options.percent_cvg_files
    assembler_version = options.assembler_version
    project_name = options.project_name

    pangolin_lineage_csv = options.pangolin_lineage_csv
    # pangolin_version = options.pangolin_version

    nextclade_clades_csv = options.nextclade_clades_csv
    nextclade_variants_csv = options.nextclade_variants_csv
    nextclade_version = options.nextclade_version

    # create lists from the column table txt file input
    sample_name_list = create_list_from_write_lines_input(write_lines_input=sample_name_array)
    cov_out_file_list = create_list_from_write_lines_input(write_lines_input = cov_out_files)
    percent_cvg_file_list = create_list_from_write_lines_input(write_lines_input=percent_cvg_files)
    
    # concat cov_out files and percent_cvg files
    cov_out_df = concat_cov_out(cov_out_file_list=cov_out_file_list)
    percent_cvg_df = concat_percent_cvg(percent_cvg_file_list=percent_cvg_file_list)

    # get df of relavant spike mutations (inlcuding 69/70 del) from nextclade file
    spike_variants_df = get_df_spike_mutations(variants_csv = nextclade_variants_csv)

    # create results file
    results_df = concat_results(sample_name_list = sample_name_list,
                                workbook_path = workbook_path,
                                project_name = project_name,
                                assembler_version = assembler_version,
                                pangolin_lineage_csv=pangolin_lineage_csv,
                                nextclade_clades_csv=nextclade_clades_csv,
                                nextclade_version=nextclade_version,
                                cov_out_df=cov_out_df,
                                percent_cvg_df=percent_cvg_df, 
                                spike_variants_df = spike_variants_df)
    
    # create wgs horizon output
    make_wgs_horizon_output(project_name = project_name,
                            results_df=results_df)
    

    print('DONE!')


#def get_df_spike_mutations(variants_csv):

#     def get_accession_id(fasta_header):
#         accession_id = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
#         return accession_id

#     variants = pd.read_csv(variants_csv, dtype = {'accession_id' : object})
#     if variants.shape[0] == 0 :
#         print('there are no accession_ids in nextclade variants file; likely all your sequences were bad')
#     variants = variants.rename(columns = {'accession_id' : 'fasta_header'})

#     accession_id = variants.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
#     variants.insert(value = accession_id, loc = 0, column = 'accession_id')
#     variants = variants.drop(columns = 'fasta_header')

#     #### filter variants for spike protein varaints in rbd and pbcs #####

#     crit = variants.gene == 'S'
#     critRBD = (variants.codon_position >= 461) & (variants.codon_position <= 509)
#     critPBCS = (variants.codon_position >= 677) & (variants.codon_position <= 694)
#     crit732 = variants.codon_position == 732
#     crit452 = variants.codon_position == 452
#     crit253 = variants.codon_position == 253
#     crit13 = variants.codon_position == 13
#     crit145 = variants.codon_position == 145 # delta plus AY.4.2
#     crit222 = variants.codon_position == 222 # delta plus AY.4.2
#     critdel = variants.variant_name.str.contains('del')

#     variants_general = variants[crit & (critRBD | critPBCS | crit732 | critdel | crit452 | crit253 | crit13 | crit145 | crit222)]

#     ##### special filter for omicron variants ######
#     omicron_spike_mutations_list = ['S_G339D', 'S_S371L', 'S_S373LP', 'S_S375F', 'S_N440K', 'S_G446S', 'S_E484K', 'S_Q493K',
#                                    'S_Q486S', 'S_Q498R', 'S_Y505H', 'S_T547K', 'S_N764K', 'S_N856K', 'S_Q954H', 'S_N969K', 'S_L981F',
#                                    '_ins22205GAGCCAGAA', 'S_G142del', 'S_V143del', 'S_Y144del', 'S_N211del']

#     crit_omicron = variants.variant_name.isin(omicron_spike_mutations_list)
#     crit_omicron_insertion_1 = variants.variant_name.str.contains('ins')
#     crit_omicron_insertion_2 = variants.codon_position == 22205

#     omicron_variants = variants[crit_omicron | (crit_omicron_insertion_1 & crit_omicron_insertion_2) ]


#     #### special filter for delta + (AY.4.2) mutations #####
#     delta_plus_spike_mutations_list = ['S_Y145H', 'S_A222V']
#     crit_delta_plus = variants.variant_name.isin(delta_plus_spike_mutations_list)
#     delta_plus_variants = variants[crit_delta_plus]


#     ####### get list of spike mutations in the poly cleavage site and the rbd
#     ## use variants_general
#     accession_ids = variants_general.accession_id.unique().tolist()

#     df = pd.DataFrame()
#     accession_id_list = []
#     variant_name_list = []

#     seperator = '; '

#     for accession_id in accession_ids:
#         accession_id_list.append(accession_id)

#         crit = variants_general.accession_id == accession_id
#         f = variants_general[crit]
#         f = f.reset_index()

#         mutations = []
#         for row in range(f.shape[0]):
#             mutations.append(f.variant_name[row])
#         mutations_string = seperator.join(mutations)
#         variant_name_list.append(mutations_string)

#     df['accession_id'] = accession_id_list
#     df['spike_mutations'] = variant_name_list


#     #### make special column for omicron varints ######
#     ## use omicron_variants

#     if omicron_variants.shape[0] != 0:
#         accession_ids = omicron_variants.accession_id.unique().tolist()

#         df_omicron = pd.DataFrame()
#         accession_id_list = []
#         omicron_variant_list = []

#         seperator = '; '

#         for accession_id in accession_ids:
#             accession_id_list.append(accession_id)

#             crit = omicron_variants.accession_id == accession_id
#             f = omicron_variants[crit]
#             f = f.reset_index()

#             mutations = []
#             for row in range(f.shape[0]):
#                 mutations.append(f.variant_name[row])
#             mutations_string = seperator.join(mutations)
#             omicron_variant_list.append(mutations_string)

#         df_omicron['accession_id'] = accession_id_list
#         df_omicron['omicron_spike_mutations'] = omicron_variant_list

#         df_omicron = df_omicron.set_index('accession_id')
#         df = df.set_index('accession_id')
#         df = df.join(df_omicron, how = 'left')
#         df = df.reset_index()

#     else:
#         df['omicron_spike_mutations'] = ''

#     #### make special column for delta plus varints ######
#     ### use delta_plus_variants
#     if delta_plus_variants.shape[0] != 0:
#         accession_ids = delta_plus_variants.accession_id.unique().tolist()

#         df_delta_plus = pd.DataFrame()
#         accession_id_list = []
#         delta_plus_variants_list = []

#         seperator = '; '

#         for accession_id in accession_ids:
#             accession_id_list.append(accession_id)

#             crit = delta_plus_variants.accession_id == accession_id
#             f = delta_plus_variants[crit]
#             f = f.reset_index()

#             mutations = []
#             for row in range(f.shape[0]):
#                 mutations.append(f.variant_name[row])
#             mutations_string = seperator.join(mutations)
#             delta_plus_variants_list.append(mutations_string)

#         df_delta_plus['accession_id'] = accession_id_list
#         df_delta_plus['delta_plus_spike_mutations'] = delta_plus_variants_list

#         df_delta_plus = df_delta_plus.set_index('accession_id')
#         df = df.set_index('accession_id')
#         df = df.join(df_delta_plus, how = 'left')
#         df = df.reset_index()

#     else:
#         df['delta_plus_spike_mutations'] = ''


#     return df