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
    options = parser.parse_args(args)
    return options


if __name__ == '__main__':
    
    options = getOptions()
    project_name = options.project_name
    pangolin_version = options.pangolin_version
    nextclade_version = options.nextclade_version
    analysis_date = options.analysis_date
    workflow_version = options.workflow_version

    df = pd.DataFrame()

    # begin to fill in table
    df.at[0, 'software'] = 'SC2_lineage_calling_and_results'
    df.at[0, 'associated_docker_container'] = ''
    df.at[0, 'version'] = workflow_version

    df.at[1, 'software'] = 'pangolin_version'
    df.at[1, 'associated_docker_container'] = 'staphb/pangolin'
    df.at[1, 'version'] = pangolin_version

    df.at[2, 'software'] = 'nextclade_version'
    df.at[2, 'associated_docker_container'] = 'nextstrain/nextclade'
    df.at[2, 'version'] = nextclade_version

    # add project name and anaysis date
    df['project_name'] = project_name
    df['analysis_date'] = analysis_date

    col_order = ['project_name', 'analysis_date', 'software', 'associated_docker_container', 'version']
    df = df[col_order]
    
    outfile = f'version_capture_lineage_calling_and_results_{project_name}_v{workflow_version}.csv'
    df.to_csv(outfile, index = False)
