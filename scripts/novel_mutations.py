# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:14:03 2023

@author: aesmith
"""
import pandas as pd
import numpy as np
import argparse
import sys
import math


def parse_arguments(args = sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('--project_names', nargs = '+', help = 
                        'space-separated list of project names')
    parser.add_argument('--combined_mutations_files', nargs = '+', help = 
                        'space-separated list of combined_mutations_files')
    parser.add_argument('--historical_full', help = 'workspace reference for full\
                        historical mutations data')
    parser.add_argument('--historical_unique', help = 'workspace reference for unique\
                        historical mutations data')
    parser.add_argument('--metadata', help = 'workspace reference for \
                        wastewater metadata')
    parser.add_argument('--gff', help = 'workspace reference for formatted gff file')
    parser.add_argument('--today', help = 'today\'s date in %YYYY-%mm-%dd format')
    parser.add_argument('--sites_to_drop', nargs = '*', help = 'space-separated \
                        list of any wastewater facility site ids to not include')

    parsed_args = parser.parse_args(args)
    
    return parsed_args


def read_reference_files(historical_full, historical_unique, metadata, gff):
   
    df_historical_full = pd.read_csv(historical_full, sep = '\t', 
                         parse_dates = ['collection_date'], 
                         dtype = {'sample_name' : 'category', 
                                  'sample_name_base' : 'category'})
    df_historical_unique = pd.read_csv(historical_unique, sep = '\t', 
                           parse_dates = ['date_first_detected', 'date_last_detected'])
    df_metadata = pd.read_csv(metadata, sep = '\t', parse_dates = ['collection_date'])
    df_gff = pd.read_csv(gff, sep = '\t')
    
    # Change some columns to categorical to decrease memory usage
    hist_category_cols = ['ref_nucl', 'alt_nucl', 'ref_aa', 'alt_aa', 
                          'gff_feature', 'site_id']
    df_historical_full[hist_category_cols] = df_historical_full[hist_category_cols].astype('category')
    
    hist_category_cols = hist_category_cols[:-1]
    hist_category_cols.extend(['mutation_type', 'id_return', 'parent_id', 'parent', 
                               'parent_start', 'parent_end', 'indel_length'])
    df_historical_unique[hist_category_cols] = df_historical_unique[hist_category_cols].astype('category')
        
    return df_historical_full, df_historical_unique, df_metadata, df_gff


def check_collection_date_exists(df, project_name):
    """ Check that all samples have collection dates. List comprehension with 
    set to check for samples duplicated across runs. Output file and exit if 
    any samples are missing dates. """
    
    df_dates = df[~df.collection_date.isna()]
    samples_with_dates = df_dates['sample_name'].tolist()
    samples_no_dates = df[df.collection_date.isna()]['sample_name'].unique().tolist()
    s = set(samples_with_dates)
    no_dates = [x for x in samples_no_dates if x not in s]
    
    # Exit if sample doesn't have collection date in metadata sheet
    if len(no_dates) > 0:
        print('Error! The following sample(s) are not in the wwt_metadata sheet.')
        error_list = []
        for i in no_dates:
            l = df[df.sample_name == i].iloc[0][['sample_name', 'project_name']].tolist()
            error_list.append(l)
            
        no_dates.to_csv(f'{project_name}_missing_dates.tsv', sep = '\t', index = False)
        sys.exit(f'At least one sample from {project_name} didn\'t have a \
collection date associated with it: {error_list}')
       
    return df_dates


def get_mutation_type(alt_nuc):
    """ Lambda function to get mutation type (deletion, insertion, or snp) """

    if '-' in alt_nuc:
        mut_type = 'deletion'
    elif '+' in alt_nuc:
        mut_type = 'insertion'
    else:
        mut_type = 'snp'
    return mut_type


def get_parent(position, feature):
    """ Lambda function to get parent based off gff feature """
    
    # Return other info for ORF1ab since outbreak.info separates them into 
    # ORF1a and ORF1b instead of ensembl's ORF1ab and ORF1a
    if feature == 'cds-YP_009725295.1' or feature == 'cds-YP_009724389.1': #ORF1a
        if position <= 13468:
            d = ['ORF1a polyprotein', 'gene-GU280_gp01', 'ORF1a', 266, 13468]
        else:
            d = ['ORF1b polyprotein', 'gene-GU280_gp01', 'ORF1b', 13468, 21555]
        parent = pd.Series(data = d)

    else:
        parent_info = df_gff[df_gff.id == feature][['id_return', 'parent_id', 'parent']].squeeze()
        parent_data = df_gff[df_gff.id == parent_info.parent_id][['start', 'end']].squeeze()
        parent = pd.concat([parent_info, parent_data], ignore_index = True)

    return parent


def get_gff_id(position):
    """ Get gff feature for positions in CDS """
    cds = df_gff[df_gff['type'] == 'CDS']
    
    # Find all genes it lies in
    df = cds[(cds.start <= position) & (cds.end >= position)]
    # Pick the one with the shortest length if there is overlap
    feature = df[df.length == df.length.min()].id
    if feature.shape[0] == 0:
        feature = None
    else:
        feature = feature.values[0]
    return feature


def assign_gene_coordinates(row):
    """
    Convert genomic coordinate to gene coordinate. Formula is as follows:
        floor((position - parent start) / 3) + 1
        ex: Mutation S:A67V genomic coordinates: 21,761-21,763
            Spike position: 21,563-25,384
            ((21761-21563) / 3) + 1 = 67
            
    The floor must be taken as the nucleotide position given could be at any of
    the three available positions for that amino acid change.

    Parameters
    ----------
    row : dataframe row
        Should include:
            position - genomic position of mutation
            ref_aa - reference amino acid
            alt_aa - alternate amino acid
            parent - parent gene name
            parent_start - genomic position start of parent gene

    Returns
    -------
    gene_coordinate : string
        Calculated gene coordinate

    """
    # Indels start at position after reference position given by IGV
    if row.mutation_type == 'insertion' or row.mutation_type == 'deletion':
        position = row.position + 1
    else:
        position = row.position

    nuc_pos_on_gene = position - row.parent_start + 1
    
    gene_coordinate_number = math.floor((nuc_pos_on_gene - 1)/ 3) + 1
    gene_coordinate_str = str(gene_coordinate_number)
    
    match row.mutation_type:
        
        case 'snp':
            # Some ref/alt aa are nan
            if type(row.alt_aa[0]) is float:
                gene_coordinate =  f'{row.parent}:{gene_coordinate_str}'
            else: 
                gene_coordinate = f'{row.parent}:{row.ref_aa}{gene_coordinate_str}{row.alt_aa}'
                
        case 'insertion':
            gene_coordinate = f'{row.parent}:ins{gene_coordinate_str}'
            
        case 'deletion':
            codon_position = nuc_pos_on_gene % 3
            del_nuc_length = len(row.alt_nucl[1:])

            # If the deletion starts/ends in the middle of an aa, add 1 to gene 
            # coordinate number
            if codon_position != 1 or del_nuc_length < 3:
                gene_coordinate_number += 1
                gene_coordinate_str = str(gene_coordinate_number)

            # If longer than one amino acid
            if del_nuc_length > 3:
                del_aa_length = math.ceil(del_nuc_length / 3)
                del_end = gene_coordinate_number + del_aa_length - 1
                gene_coordinate = f'{row.parent}:del{gene_coordinate_str}-{del_end}'
                
            else:
                gene_coordinate = f'{row.parent}:del{gene_coordinate_str}'
         
    return gene_coordinate


def assign_nuc_coordinates(row):
    """ Lambda function to calculcate nucleotide coordinate for mutation """
    match row.mutation_type:
        case 'insertion':
            position = row.position + 1
            coordinate = f'ins{str(position)}{row.alt_nucl[1:]}'
        case 'deletion':
            position = row.position + 1
            if len(row.alt_nucl) > 2:
                del_length = len(row.alt_nucl) - 1
                del_end = position + del_length - 1
                coordinate = f'del{str(position)}-{str(del_end)}'
            else:
                coordinate = f'del{str(position)}'
        case 'snp':
            coordinate = f'{row.ref_nucl}{str(row.position)}{row.alt_nucl}'
            
    return coordinate


def unique_counts_and_dates(df):
    """
    Calculate counts, earliest date seen, most recent date seen, and length
    of time seen for each unique mutation

    Parameters
    ----------
    df : dataframe
        full mutations.

    Returns
    -------
    df_dates : dataframe
        unique mutations df with added date features.

    """
    
    count_list = ['position', 'ref_nucl', 'alt_nucl']
    mutation_counts = df[count_list].value_counts().reset_index()
    mutation_counts.columns = [*mutation_counts.columns[:-1], 'times_detected']
    date_list = ['position', 'ref_nucl', 'alt_nucl', 'collection_date']
    df_latest_dates = df[date_list].sort_values(by = ['collection_date'], 
                         ascending = False).drop_duplicates(subset = 
                          count_list).reset_index(drop = True).rename(columns = 
                                      {'collection_date' : 'date_last_detected'})
                                                                                  
    df_earliest_dates = df[date_list].sort_values(by = ['collection_date'], 
                         ascending = True).drop_duplicates(subset = 
                          count_list).reset_index(drop = True).rename(columns = 
                                      {'collection_date' : 'date_first_detected'})      
                  
    df_dates = df_earliest_dates.merge(df_latest_dates)
                                                                                  
    df_dates = df_dates.merge(mutation_counts, on = count_list)
    df_dates['length_of_time_seen'] = (df_dates.date_last_detected - 
                                       df_dates.date_first_detected).dt.days + 1
    return df_dates


def parse_project_mutations(project_dict, df_metadata):
    """
    Parse passed combined mutations files. Drop any sites indicated, freyja 
    failed samples, duplicates from replicates, and merge with metadata to add 
    collection date and site ids. Drop verification and control samples. Output 
    small file of unique mutations for each project. Return df of all 'pass' 
    samples across projects with collection dates and site ids.

    Parameters
    ----------
    project_dict : dictionary
        Keys: project names. Values: paths of combined_mutations files.
        
    df_metadata : dataframe
        Metadata file.

    Returns
    -------
    dfs : dataframe
        All samples across projects

    """
    df_list = []
    mut_cols = ['sample_name', 'position', 'ref_nucl', 'alt_nucl', 'ref_aa', 
                'alt_aa', 'pass', 'gff_feature']
    mut_category_cols = mut_cols[:1] + mut_cols[2:6] + mut_cols[-1:]
    
    for project in project_dict:
        path = project_dict[project]
        df = pd.read_csv(path, delim_whitespace = True, usecols = mut_cols,
                         dtype = dict.fromkeys(mut_category_cols, 'category'))
        
        # Drop freyja failures
        df = df[df['pass'] == True]
        
        # Merge with metadata for collection date and site id and check that all 
        # samples have this info. Drop controls and verification samples
        df = df.merge(df_metadata, on = 'sample_name', how = 'left')
        df = df[df.sample_type == 'sample']

        # Drop specified utility site ids if any
        if len(sites_to_drop) > 0:
            df = df[~df['site_id'].isin(sites_to_drop)]
              
        df = check_collection_date_exists(df, project)
        
        df = df.drop(columns = ['pass', 'sample_type'])
        
        # Drop duplicates from replicates (based on site id and collection date)
        duplicate_cols = ['position', 'ref_nucl', 'alt_nucl', 'collection_date', 'site_id']
        df = df.drop_duplicates(subset = duplicate_cols)
        
        # Output file of unique mutations for each project
        df_unique = unique_counts_and_dates(df)
        df_unique.to_csv(f'{project}_unique_mutations.tsv', sep = '\t', index = False)
        
        df_list.append(df)
        
    if len(df_list) > 1:
        dfs = pd.concat(df_list, ignore_index = True)
    else:
        dfs = df_list[0]
        
    return dfs


def calc_recurrent_mutations(row):
    """ Lambda function to determine if mutations detected in this batch haven't
    been seen in longer than 6 months """
    time_between = (row.date_first_detected_new - row.date_last_detected_historical).days
    if time_between > 182:
        recurrent = True
    else:
        recurrent = False
        
    return [recurrent, time_between]


def recurrent_mutations(df):
    """
    For mutations in this batch, check for recurrent mutations not seen 
    historically in more than 6 months. Output file of any found.

    Parameters
    ----------
    df : dataframe
        Unique mutations for this batch, merged with historical unique dates/counts

    Returns
    -------
    None.

    """
    # Check for recurrent mutations not seen in >6mo
    recurrent_cols = ['position', 'ref_nucl', 'alt_nucl', 'ref_aa', 'alt_aa', 
                      'gff_feature', 'mutation_type', 'gene_coordinate', 'nuc_coordinate', 
                      'date_first_detected_new', 'date_last_detected_new', 
                      'times_detected_new', 'date_first_detected_historical', 
                      'date_last_detected_historical', 'times_detected_historical']
    df_recurrent_mutations = df[recurrent_cols].copy()
    df_recurrent_mutations[['recurrent', 'days_between']] = df_recurrent_mutations[['date_first_detected_new', 
                       'date_last_detected_historical']].apply(lambda x: calc_recurrent_mutations(x), 
                           axis = 1, result_type = 'expand')
    df_recurrent_mutations = df_recurrent_mutations[df_recurrent_mutations.recurrent == True]
    df_recurrent_mutations.drop(columns = ['recurrent'], inplace = True)
    df_recurrent_mutations.to_csv(f'recurrent_mutations_{today}.tsv', sep = '\t', index = False)
    
    return 


def add_features_new_mutations(df, df_new_full, cols):
    """
    Add features to mutations not seen before. Fill in ref_aa, alt_aa, and
    gff_feature. Fill in blank gff features if in CDS. Calculate type, 
    indel_length, parent gene info, gene coordinate, and nucleotide coordinate.

    Parameters
    ----------
    df_temp : dataframe
        new unique mutations .
    df_new_full : dataframe
        new full mutations.

    Returns
    -------
    df : dataframe
        new unique mutations with all features filled in.

    """

    df.columns = df.columns.str.replace('_new', '')
    subset_cols = ['position', 'ref_nucl', 'alt_nucl']
    gff_merge_cols = ['position', 'ref_nucl', 'alt_nucl', 'alt_aa', 'ref_aa', 'gff_feature']
    df = df.merge(df_new_full[gff_merge_cols].drop_duplicates(subset = subset_cols), how = 'left')
    
    # Columns to calculate
    # type
    df['mutation_type'] = df.apply(lambda x: get_mutation_type(x['alt_nucl']), axis = 1)
    
    # indel length
    df['indel_length'] = df.apply(lambda x: len(x['alt_nucl']) - 1 
                                  if x['mutation_type'] != 'snp' else np.nan, axis = 1)
    
    # Fill in blank gff_features if applicable
    df.gff_feature.fillna(df.apply(lambda x: get_gff_id(x.position), axis = 1), 
                          inplace = True)
    
    # Get parent gene info for those with gff_features
    parent_cols = ['id_return', 'parent_id', 'parent', 'parent_start', 'parent_end']
    df[parent_cols] = df.apply(lambda x: get_parent(x.position, x.gff_feature) 
                               if(np.all(pd.notnull(x['gff_feature']))) else np.nan, axis = 1)
    
    # Get gene and nuc coordinates
    coordinate_cols = ['position', 'ref_nucl', 'alt_nucl', 'ref_aa', 'alt_aa', \
                       'parent', 'parent_start', 'parent_end', 'mutation_type']
    df['gene_coordinate'] = df.apply(lambda x: 
                                     assign_gene_coordinates(x[coordinate_cols]) 
                                     if(np.all(pd.notnull(x['gff_feature']))) else np.nan, axis = 1)
    
    df['nuc_coordinate'] = df.apply(lambda x: assign_nuc_coordinates(x), axis = 1)
    
    df = df[cols]

    df.to_csv(f'novel_mutations_{today}.tsv', sep = '\t', index = False)
    
    return df


def update_features_prev_seen_mutations(df, cols):
    """
    Update dates/counts for previously seen mutations. Check for recurrent 
    mutations and output file if any found.

    Parameters
    ----------
    df : dataframe
        Unique mutations in this batch that have been previously seen, merged
        with historical data for all dates/counts
    cols : list
        Final column order list

    Returns
    -------
    df : dataframe
        Initial dataframe with updated dates/features

    """
    df = df.loc[:, ~df.columns.str.contains('^length_of_time_seen')]
    
    df['date_first_detected'] = df.loc[:, df.columns.str.contains('^date_first_detected_')].min(axis = 1)
    df['date_last_detected'] = df.loc[:, df.columns.str.contains('^date_last_detected_')].max(axis = 1)
    df['times_detected'] = df.times_detected_new + df.times_detected_historical
    
    # Check for recurrent mutations
    recurrent_mutations(df)
    
    # Remove unnecessary columns
    df = df.loc[:, ~df.columns.str.contains('^date_(fir|la)st_detected_')]
    df = df.loc[:, ~df.columns.str.contains('^times_detected_')]
    
    df['length_of_time_seen'] = (df.date_last_detected - 
                                       df.date_first_detected).dt.days + 1
    df = df.astype({'times_detected' : 'int64'})

    # Put columns back in order
    df = df[cols]
    
    return df


def update_unique_data(df_new_unique, df_new_full, df_historical_unique):
    """
    Update historical unique data with new data

    Parameters
    ----------
    df_new_unique : dataframe
        new unique data before updated/filled features.
    df_new_full : dataframe
        new full data.
    df_historical_unique : dataframe
        historical unique data.

    Returns
    -------
    df_unique : dataframe
        updated historical unique data with all new mutation information.

    """
    # Temp df for appending/recalculating unique mutation info
    duplicate_cols = ['position', 'ref_nucl', 'alt_nucl']
    df_temp = df_new_unique.merge(df_historical_unique, on = duplicate_cols, 
                                  how = 'left', suffixes = ['_new', '_historical'])
    df_temp = df_temp.astype({'alt_nucl' : 'category'})
    
    # Have column order 
    cols = df_historical_unique.columns
    
    # Subset dataframe into novel and previously seen
    has_features = df_temp[~df_temp.nuc_coordinate.isna()]
    filled_in_cols = ['position', 'ref_nucl', 'alt_nucl', 'date_first_detected_new', 
                      'date_last_detected_new', 'times_detected_new', 'length_of_time_seen_new']
    needs_features = df_temp[df_temp.nuc_coordinate.isna()][filled_in_cols]
    
    # Update/fill in features of new mutations. Output files of novel and recurrent mutations
    df_new_updated_features = update_features_prev_seen_mutations(has_features, cols)
    df_new_filled_features = add_features_new_mutations(needs_features, df_new_full, cols)
    
    # Update historical file, first with updated mutations and then with new mutations
    static_cols = ['position', 'ref_nucl', 'alt_nucl', 'ref_aa', 'alt_aa', 
                   'gff_feature', 'mutation_type', 'id_return', 'parent_id', 'parent',
                   'parent_start', 'parent_end', 'gene_coordinate', 
                   'nuc_coordinate', 'indel_length']
    df_unique_updated = df_historical_unique.merge(df_new_updated_features, on = 
                       static_cols, suffixes = ['_old', '_updated'], how = 'left')
    
    updated_cols = ['date_first_detected', 'date_last_detected', 'times_detected', 
                    'length_of_time_seen']
    for c in updated_cols:
        old_col = c + '_old'
        updated_col = c + '_updated'
        df_unique_updated[updated_col].fillna(df_unique_updated[old_col], inplace = True)
        df_unique_updated.drop(columns = [old_col], inplace = True)
        df_unique_updated.rename(columns = {updated_col : c}, inplace = True)
        
    df_unique = pd.concat([df_unique_updated, df_new_filled_features], ignore_index = True)
    
    return df_unique


def main():
    
    # For each project, create file of full passes and file of unique mutations
    project_dict = dict(zip(projects, mutations_files))
    df_new_full = parse_project_mutations(project_dict, df_metadata)
    
    # Append full data to historical
    df_full = pd.concat([df_historical_full, df_new_full], ignore_index = True)
    
    # Check for duplicates again in case samples were re-run
    df_full = df_full.drop_duplicates(subset = ['position', 'ref_nucl', 'alt_nucl', 
                                'collection_date', 'site_id']).reset_index(drop = True)
    
    # For concatenated new mutation data, get dates and counts
    df_new_unique = unique_counts_and_dates(df_new_full)
    
    # Update historical unique data with new information
    df_unique = update_unique_data(df_new_unique, df_new_full, df_historical_unique)
    
    # Write all files to .tsv (updated unique, updated full, each project unique)
    df_unique.to_csv('novel_mutations_historical_unique.tsv', sep = '\t', index = False)
    df_full.to_csv('novel_mutations_historical_full.tsv', sep = '\t', index = False)


if __name__ == '__main__':
    
    # Parse args
    args = parse_arguments()
    projects = args.project_names
    mutations_files = args.combined_mutations_files
    historical_full = args.historical_full
    historical_unique = args.historical_unique
    metadata = args.metadata
    gff = args.gff
    today = args.today
    sites_to_drop = args.sites_to_drop
    
    # Read in reference files
    df_historical_full, df_historical_unique, df_metadata, df_gff = \
        read_reference_files(historical_full, historical_unique, metadata, gff)
        
    main()

