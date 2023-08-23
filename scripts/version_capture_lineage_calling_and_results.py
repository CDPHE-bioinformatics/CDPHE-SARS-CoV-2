import argparse
import sys
import pandas as pd
import numpy as np


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument("--project_name",  help= "")
    parser.add_argument("--pangolin_version",  help= "")
    parser.add_argument("--nextclade_version",  help= "")
    parser.add_argument("--analysis_date",  help= "")
    parser.add_argument('--workflow_version', help = '')
    parser.add_argument('--pangolin_lineage_csv')
    options = parser.parse_args(args)
    return options


if __name__ == '__main__':
    
    options = getOptions()
    project_name = options.project_name
    pangolin_version = options.pangolin_version
    nextclade_version = options.nextclade_version
    analysis_date = options.analysis_date
    workflow_version = options.workflow_version
    pangolin_lineage_csv = options.workflow_version

    # read in pangolin_lineage_csv and pull out versions in file
    pangolin_df = pd.read_csv(pangolin_lineage_csv)
    pango_designation_version = pangolin_df.version.unique().tolist()[0]
    pangolin_scorpio_version = pangolin_df.scorpio_version.unique().tolist()[0]
    pangolin_constellation_version = pangolin_df.constellation_version.unique().tolist()[0]


    df = pd.DataFrame()

    # begin to fill in table
    df.at[0, 'software'] = 'SC2_lineage_calling_and_results'
    df.at[0, 'associated_docker_container'] = ''
    df.at[0, 'version'] = workflow_version

    df.at[1, 'software'] = 'pangolin_version'
    df.at[1, 'associated_docker_container'] = 'staphb/pangolin'
    df.at[1, 'version'] = pangolin_version

    df.at[2, 'software'] = 'pango_designation'
    df.at[2, 'associated_docker_container'] = 'staphb/pangolin'
    df.at[2, 'version'] = pango_designation_version

    df.at[3, 'software'] = 'pangolin_scropio'
    df.at[3, 'associated_docker_container'] = 'staphb/pangolin'
    df.at[3, 'version'] = pangolin_scorpio_version

    df.at[4, 'software'] = 'pangolin_constellation'
    df.at[4, 'associated_docker_container'] = 'staphb/pangolin'
    df.at[4, 'version'] = pangolin_constellation_version

    df.at[5, 'software'] = 'pangolin_version'
    df.at[5, 'associated_docker_container'] = 'staphb/pangolin'
    df.at[5, 'version'] = pangolin_version

    df.at[6, 'software'] = 'nextclade_version'
    df.at[6, 'associated_docker_container'] = 'nextstrain/nextclade'
    df.at[6, 'version'] = nextclade_version

    # add project name and anaysis date
    df['project_name'] = project_name
    df['analysis_date'] = analysis_date

    col_order = ['project_name', 'analysis_date', 'software', 'associated_docker_container', 'version']
    df = df[col_order]
    
    outfile = f'version_capture_lineage_calling_and_results_{project_name}_{workflow_version}.csv'
    df.to_csv(outfile, index = False)
