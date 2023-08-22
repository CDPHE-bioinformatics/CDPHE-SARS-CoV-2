import argparse
import sys
import pandas as pd
import numpy as np


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument("--project_name",  help= "")
    parser.add_argument("--samtools_version_staphb",  help= "")
    parser.add_argument("--samtools_version_andersenlabapps",  help= "")
    parser.add_argument("--ivar_version",  help= "")
    parser.add_argument("--freyja_version",  help= "")
    parser.add_argument("--analysis_date",  help= "")
    parser.add_argument('--workflow_version', help = '')
    options = parser.parse_args(args)
    return options


if __name__ == '__main__':
    
    options = getOptions()
    project_name = options.project_name
    samtools_version_staphb = options.samtools_version_staphb
    samtools_version_andersenlabapps = options.samtools_version_andersenlabapps
    ivar_version = options.ivar_version
    freyja_version = options.freyja_version
    analysis_date = options.analysis_date
    workflow_version = options.workflow_version

    df = pd.DataFrame()

    # begin to fill in table
    df.at[0, 'software'] = 'SC2_wastewater_variant_calling'
    df.at[0, 'associated_docker_container'] = ''
    df.at[0, 'version'] = workflow_version

    df.at[1, 'software'] = 'samtools'
    df.at[1, 'associated_docker_container'] = 'staphb/samtools'
    df.at[1, 'version'] = samtools_version_staphb

    df.at[2, 'software'] = 'samtools'
    df.at[2, 'associated_docker_container'] = 'andersenlabapps/ivar'
    df.at[2, 'version'] = samtools_version_andersenlabapps

    df.at[3, 'software'] = 'ivar'
    df.at[3, 'associated_docker_container'] = 'andersenlabapps/ivar'
    df.at[3, 'version'] = ivar_version

    df.at[4, 'software'] = 'freyja'
    df.at[4, 'associated_docker_container'] = 'staphb/freyja'
    df.at[4, 'version'] = freyja_version

    # add project name and anaysis date
    df['project_name'] = project_name
    df['analysis_date'] = analysis_date

    col_order = ['project_name', 'analysis_date', 'software', 'associated_docker_container', 'version']
    df = df[col_order]
    
    outfile = f'version_capture_wastewater_variant_calling_{project_name}_v{workflow_version}.csv'
    df.to_csv(outfile, index = False)
