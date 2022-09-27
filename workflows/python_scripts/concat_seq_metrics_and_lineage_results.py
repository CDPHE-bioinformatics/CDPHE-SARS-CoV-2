#! /usr/bin/env python

# updated 2022-03-04
## added assembler version

#updated 2022-04-04
## adjust columns to account for pangolin v4.0 major update

# updated 2022-07-07
## gets rid of the report_to_epi stuff in the wgs_horizon_csv (the column will just be blank).

import argparse
import sys
import pandas as pd
import numpy as np
from datetime import date
import re
# import os


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")

    parser.add_argument('--sample_list')
    parser.add_argument('--plate_name_file_list')
    parser.add_argument('--plate_sample_well_file_list')
    parser.add_argument('--primer_set_file_list')
    parser.add_argument('--tech_platform_file_list')
    parser.add_argument('--read_type_file_list')
    parser.add_argument("--bam_file_list",  help= "txt file with list of bam file paths")
    parser.add_argument('--percent_cvg_file_list', help = 'txt file with list of percent cvg file paths')

    parser.add_argument('--pangolin_lineage_csv', help = 'csv output from pangolin')
    parser.add_argument('--pangolin_version')

    parser.add_argument('--assembler_version_table_list')

    parser.add_argument('--nextclade_clades_csv', help = 'csv output from nextclade parser')
    parser.add_argument('--nextclade_variants_csv')
    parser.add_argument('--nextclade_version')

    parser.add_argument('--seq_run_file_list')
    options = parser.parse_args(args)
    return options


def concat_samtools(bam_file_list):

    # get the input paths from text file into a python list
    with open(bam_file_list, 'r') as f:
        file_paths = []
        for line in f:
            file_paths.append(line.strip())

     # initiate dataframe for concatenation
    df = pd.DataFrame()
    accession_id_list = []
    samtools_numreads_list = []
    samtools_depth_list = []
    samtools_baseq_list = []
    samtools_mapq_list = []

    #loop through bam file stats files and pull data
    for file_path in file_paths:
        d = pd.read_csv(file_path, sep = '\t')
        if re.search('barcode', file_path):
            # for nanopore runs
            sequence_name = re.findall('/([0-9a-zA-Z_\-\.]+)_barcode', file_path)[0]
        else:
            # for illumina runs
            sequence_name = re.findall('/([0-9a-zA-Z_\-\.]+)_coverage.txt', file_path)[0]

        # pull data from samtools output
        num_reads = d.numreads[0]
        depth = d.meandepth[0]
        baseq = d.meanbaseq[0]
        mapq = d.meanmapq[0]

        accession_id_list.append(sequence_name)
        samtools_numreads_list.append(num_reads)
        samtools_depth_list.append(depth)
        samtools_baseq_list.append(baseq)
        samtools_mapq_list.append(mapq)

    df['accession_id'] = accession_id_list
    df['num_reads'] = samtools_numreads_list
    df['mean_depth'] = samtools_depth_list
    df['mean_base_quality'] = samtools_baseq_list
    df['mean_map_quality'] = samtools_mapq_list

    return df

def concat_percent_cvg(percent_cvg_file_list):

    # get the input paths from text file into a python list
    with open(percent_cvg_file_list, 'r') as f:
        file_paths = []
        for line in f:
            file_paths.append(line.strip())

    df_list = []
    for file in file_paths:
        d = pd.read_csv(file, dtype = {'accession_id' : object})
        df_list.append(d)

    df = pd.concat(df_list)

    return df


def get_df_spike_mutations(variants_csv):

    def get_accession_id(fasta_header):
        accession_id = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return accession_id

    variants = pd.read_csv(variants_csv, dtype = {'accession_id' : object})
    if variants.shape[0] == 0 :
        print('there are no accession_ids in nextclade variants file; likely all your sequences were bad')
    variants = variants.rename(columns = {'accession_id' : 'fasta_header'})

    accession_id = variants.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
    variants.insert(value = accession_id, loc = 0, column = 'accession_id')
    variants = variants.drop(columns = 'fasta_header')

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
    critdel = variants.variant_name.str.contains('del')

    variants_general = variants[crit & (critRBD | critPBCS | crit732 | critdel | crit452 | crit253 | crit13 | crit145 | crit222)]

    ##### special filter for omicron variants ######
    omicron_spike_mutations_list = ['S_G339D', 'S_S371L', 'S_S373LP', 'S_S375F', 'S_N440K', 'S_G446S', 'S_E484K', 'S_Q493K',
                                   'S_Q486S', 'S_Q498R', 'S_Y505H', 'S_T547K', 'S_N764K', 'S_N856K', 'S_Q954H', 'S_N969K', 'S_L981F',
                                   '_ins22205GAGCCAGAA', 'S_G142del', 'S_V143del', 'S_Y144del', 'S_N211del']

    crit_omicron = variants.variant_name.isin(omicron_spike_mutations_list)
    crit_omicron_insertion_1 = variants.variant_name.str.contains('ins')
    crit_omicron_insertion_2 = variants.codon_position == 22205

    omicron_variants = variants[crit_omicron | (crit_omicron_insertion_1 & crit_omicron_insertion_2) ]


    #### special filter for delta + (AY.4.2) mutations #####
    delta_plus_spike_mutations_list = ['S_Y145H', 'S_A222V']
    crit_delta_plus = variants.variant_name.isin(delta_plus_spike_mutations_list)
    delta_plus_variants = variants[crit_delta_plus]


    ####### get list of spike mutations in the poly cleavage site and the rbd
    ## use variants_general
    accession_ids = variants_general.accession_id.unique().tolist()

    df = pd.DataFrame()
    accession_id_list = []
    variant_name_list = []

    seperator = '; '

    for accession_id in accession_ids:
        accession_id_list.append(accession_id)

        crit = variants_general.accession_id == accession_id
        f = variants_general[crit]
        f = f.reset_index()

        mutations = []
        for row in range(f.shape[0]):
            mutations.append(f.variant_name[row])
        mutations_string = seperator.join(mutations)
        variant_name_list.append(mutations_string)

    df['accession_id'] = accession_id_list
    df['spike_mutations'] = variant_name_list


    #### make special column for omicron varints ######
    ## use omicron_variants

    if omicron_variants.shape[0] != 0:
        accession_ids = omicron_variants.accession_id.unique().tolist()

        df_omicron = pd.DataFrame()
        accession_id_list = []
        omicron_variant_list = []

        seperator = '; '

        for accession_id in accession_ids:
            accession_id_list.append(accession_id)

            crit = omicron_variants.accession_id == accession_id
            f = omicron_variants[crit]
            f = f.reset_index()

            mutations = []
            for row in range(f.shape[0]):
                mutations.append(f.variant_name[row])
            mutations_string = seperator.join(mutations)
            omicron_variant_list.append(mutations_string)

        df_omicron['accession_id'] = accession_id_list
        df_omicron['omicron_spike_mutations'] = omicron_variant_list

        df_omicron = df_omicron.set_index('accession_id')
        df = df.set_index('accession_id')
        df = df.join(df_omicron, how = 'left')
        df = df.reset_index()

    else:
        df['omicron_spike_mutations'] = ''

    #### make special column for delta plus varints ######
    ### use delta_plus_variants
    if delta_plus_variants.shape[0] != 0:
        accession_ids = delta_plus_variants.accession_id.unique().tolist()

        df_delta_plus = pd.DataFrame()
        accession_id_list = []
        delta_plus_variants_list = []

        seperator = '; '

        for accession_id in accession_ids:
            accession_id_list.append(accession_id)

            crit = delta_plus_variants.accession_id == accession_id
            f = delta_plus_variants[crit]
            f = f.reset_index()

            mutations = []
            for row in range(f.shape[0]):
                mutations.append(f.variant_name[row])
            mutations_string = seperator.join(mutations)
            delta_plus_variants_list.append(mutations_string)

        df_delta_plus['accession_id'] = accession_id_list
        df_delta_plus['delta_plus_spike_mutations'] = delta_plus_variants_list

        df_delta_plus = df_delta_plus.set_index('accession_id')
        df = df.set_index('accession_id')
        df = df.join(df_delta_plus, how = 'left')
        df = df.reset_index()

    else:
        df['delta_plus_spike_mutations'] = ''


    return df


def concat_results(sample_list, plate_name_file_list, plate_sample_well_file_list, primer_set_file_list,
                   tech_platform_file_list, read_type_file_list,
                   samtools_df, percent_cvg_df, spike_mut_df, nextclade_clades_csv,
                   pangolin_lineage_csv, next_version, pangolin_version, seq_run_file_list, assembler_version_table_list):


    # get the list of samples and create a df using the list of samples
    all_sample_accession_ids = []
    with open(sample_list, 'r') as f:
        for line in f:
            all_sample_accession_ids.append(line.strip())


    # print(all_sample_accession_ids)
    df = pd.DataFrame(all_sample_accession_ids)
    df = df.rename(columns = {0:'accession_id'})
    df = df.set_index('accession_id')
    # print(df)



    # create lists from file lists....and then create a column for each in results df
    plate_name_list = []
    with open(plate_name_file_list, 'r') as f:
        for line in f:
            plate_name_list.append(line.strip())

    plate_sample_well_list = []
    with open(plate_sample_well_file_list, 'r') as f:
        for line in f:
            plate_sample_well_list.append(line.strip())

    primer_set_list = []
    with open(primer_set_file_list, 'r') as f:
        for line in f:
            primer_set_list.append(line.strip())

    tech_platform_list = []
    with open(tech_platform_file_list, 'r') as f:
        for line in f:
            tech_platform_list.append(line.strip())

    read_type_list = []
    with open(read_type_file_list, 'r') as f:
        for line in f:
            read_type_list.append(line.strip())

    seq_run_list = []
    with open(seq_run_file_list, 'r') as f:
        for line in f:
            seq_run_list.append(line.strip())

    assembler_version_list = []
    with open(assembler_version_table_list, 'r') as f:
        for line in f:
            assembler_version_list.append(line.strip())
    for item in assembler_version_list:
        if item != '':
            assembler_version_txt = item
            break # exit loop once find an assembler version that isn't blank

    df['plate_name'] = plate_name_list
    df['plate_sample_well'] = plate_sample_well_list
    df['primer_set'] = primer_set_list
    df['tech_platform'] = tech_platform_list
    df['read_type'] = read_type_list
    df['seq_run'] = seq_run_list
    df['analysis_date'] = str(date.today())
    df['assembler_version'] = assembler_version_txt



    def get_accession_id(fasta_header):
        accession_id = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return accession_id

    def create_fasta_header(accession_id):
        return 'CO-CDPHE-%s' % accession_id

    # read in nextclade clade results
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

    accession_id = pangolin.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
    pangolin.insert(value = accession_id, column = 'accession_id', loc = 0)
    pangolin = pangolin.drop(columns = 'fasta_header')
    pangolin = pangolin.drop(columns = 'pangolin_version')

    pangolin['pangolin_version'] = pangolin_version
    pangolin = pangolin.set_index('accession_id')

    # pull out the pangoLEARN_version....
#     pango_learn_version = pangolin.pangolin_designation_version[0]

    # read in nextclade csv
    nextclade = pd.read_csv(nextclade_clades_csv, dtype = {'accession_id' : object})
    nextclade = nextclade.rename(columns = {'accession_id' : 'fasta_header'})


    accession_id = nextclade.apply(lambda x:get_accession_id(x.fasta_header), axis = 1)
    nextclade.insert(value = accession_id, column = 'accession_id', loc = 0)
    nextclade = nextclade.drop(columns = 'fasta_header')

    nextclade['nextclade_version'] = next_version
    nextclade = nextclade.set_index('accession_id')


    # set index on the samtools_df and percent_cvg_df and variants_df
    samtools_df = samtools_df.set_index('accession_id')
    percent_cvg_df = percent_cvg_df.set_index('accession_id')
    spike_mut_df = spike_mut_df.set_index('accession_id')


    # join
    j = df.join(percent_cvg_df, how = 'left')
    j = j.join(samtools_df, how = 'left')
    j = j.join(nextclade, how = 'left')
    j = j.join(pangolin, how = 'left')
    j = j.join(spike_mut_df, how = 'left')
    j = j.reset_index()

    # add fasta header
    fasta_header = j.apply(lambda x:create_fasta_header(x.accession_id), axis=1)
    j.insert(value = fasta_header, column = 'fasta_header', loc = 0)


    col_order = [ 'accession_id',
                 'plate_name',
                 'plate_sample_well',
                 'primer_set',
                 'percent_non_ambigous_bases',
                 'nextclade',
                 'pangolin_lineage',
                 'assembler_version',

                 'omicron_spike_mutations',
                 'delta_plus_spike_mutations',

                 'spike_mutations',
                 'total_nucleotide_mutations',
                 'total_AA_substitutions',
                 'total_AA_deletions',
                 'mean_depth',

                 'number_aligned_bases',
                 'number_non_ambigous_bases',

                 'number_seqs_in_fasta',

                 'total_nucleotide_deletions',
                 'total_nucleotide_insertions',

                 'num_reads',
                 'mean_base_quality',
                 'mean_map_quality',
                 'number_N_bases',

                 'nextclade_version',
                 'pangolin_version',
                 'pangoLEARN_conflict',
                 'pangolin_ambiguity_score',
                 'pangolin_scorpio_call',
                 'pangolin_scorpio_support',
                 'pangolin_scorpio_conflict',
                 'pangolin_scorpio_notes',
                 'pango_designation_version',
                 'pangolin_scorpio_version',
                 'pangolin_constellation_version',
                 'pangolin_is_designated',
                 'pangolin_qc_status',
                 'pangolin_qc_notes',
                 'pangolin_note',

                 'seq_run',
                 'tech_platform',
                 'read_type',
                'fasta_header',
                'analysis_date']

    j = j[col_order]

    # add in 'failed assembly" in missing columns
    j.spike_mutations = j.spike_mutations.fillna(value = '')
    j.nextclade = j.nextclade.fillna(value = '')
    j.nextclade_version = j.nextclade_version.fillna(value = next_version)
    j.pangolin_lineage = j.pangolin_lineage.fillna(value = 'Unassigned')
    j.percent_non_ambigous_bases = j.percent_non_ambigous_bases.fillna(value = 0)
    j.mean_depth = j.mean_depth.fillna(value = 0)
    j.number_aligned_bases = j.number_aligned_bases.fillna(value = 0)
    j.number_seqs_in_fasta = j.number_seqs_in_fasta.fillna(value = 0)
    j.num_reads = j.num_reads.fillna(value = 0)
    j.mean_base_quality = j.mean_base_quality.fillna(value = 0)
    j.mean_map_quality = j.mean_map_quality.fillna(value = 0)
    j.number_N_bases = j.number_N_bases.fillna(value = 29903)
    j.pangolin_version = j.pangolin_version.fillna(value = pangolin_version )
#     j.pangolin_designation_version = j.pangolin_designation_version.fillna(value = pango_learn_version )
    j.assembler_version = j.assembler_version.fillna(value = assembler_version_txt)

    outfile = '%s_sequencing_results.csv' % seq_run_list[0]
    j.to_csv(outfile, index = False)

    return {"df" : j, "seq_run" : seq_run_list[0]}


def make_assembly_metrics_csv(results_df, seq_run):

    outfile = '%s_sequence_assembly_metrics.csv' % seq_run
    results_df.to_csv(outfile, index = False)


def make_wgs_horizon_output (results_df, seq_run):

    d = results_df.rename(columns = {'percent_non_ambigous_bases' : 'percent_coverage'})

    d['report_to_epi'] = ''
    d['Run_Date'] = str(date.today())

    col_order = ['accession_id', 'percent_coverage', 'pangolin_lineage', 'pangolin_version',
                 'report_to_epi', 'Run_Date', 'pango_designation_version']
    d = d[col_order]
    d = d.rename(columns = {'pango_designation_version' : 'pangoLEARN_version'})


    outfile = "%s_wgs_horizon_report.csv" % seq_run
    d.to_csv(outfile, index = False)


if __name__ == '__main__':

    options = getOptions()

    sam_df = concat_samtools(bam_file_list = options.bam_file_list)

    percent_cvg_df = concat_percent_cvg(percent_cvg_file_list = options.percent_cvg_file_list)

    spike_mut_df = get_df_spike_mutations(variants_csv = options.nextclade_variants_csv)

    results_df_dict = concat_results(sample_list = options.sample_list,
                                plate_name_file_list = options.plate_name_file_list,
                                plate_sample_well_file_list = options.plate_sample_well_file_list,
                                primer_set_file_list = options.primer_set_file_list,
                                tech_platform_file_list = options.tech_platform_file_list,
                                read_type_file_list = options.read_type_file_list,
                                samtools_df = sam_df,
                               percent_cvg_df = percent_cvg_df,
                               spike_mut_df = spike_mut_df,
                               nextclade_clades_csv = options.nextclade_clades_csv,
                               pangolin_lineage_csv = options.pangolin_lineage_csv,
                               next_version = options.nextclade_version,
                               pangolin_version = options.pangolin_version,
                               seq_run_file_list = options.seq_run_file_list,
                                    assembler_version_table_list = options.assembler_version_table_list)

    results_df = results_df_dict['df']
    seq_run = results_df_dict['seq_run']

    make_assembly_metrics_csv(results_df = results_df, seq_run = seq_run)

    make_wgs_horizon_output(results_df = results_df, seq_run = seq_run)
