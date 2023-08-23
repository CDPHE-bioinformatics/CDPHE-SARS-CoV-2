import argparse
import sys
import pandas as pd
import numpy as np


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument("--project_name",  help= "")
    parser.add_argument("--seqyclean_version",  help= "")
    parser.add_argument("--fastqc_version",  help= "")
    parser.add_argument("--bwa_version",  help= "")
    parser.add_argument("--samtools_version_broadinstitute",  help= "")
    parser.add_argument("--ivar_version",  help= "")
    parser.add_argument("--samtools_version_andersenlabapps",  help= "")
    parser.add_argument("--samtools_version_staphb",  help= "")
    parser.add_argument("--analysis_date",  help= "")
    parser.add_argument('--workflow_version', help = '')
    options = parser.parse_args(args)
    return options


if __name__ == '__main__':
    
    options = getOptions()
    project_name = options.project_name
    seqyclean_version = options.seqyclean_version
    fastqc_version = options.fastqc_version
    bwa_version = options.bwa_version
    samtools_version_broadinstitute = options.samtools_version_broadinstitute
    ivar_version = options.ivar_version
    samtools_version_andersenlabapps = options.samtools_version_andersenlabapps
    samtools_version_staphb = options.samtools_version_staphb
    analysis_date = options.analysis_date
    workflow_version = options.workflow_version

    df = pd.DataFrame()

    # begin to fill in table
    df.at[0, 'software'] = 'SC2_illumina_pe_assembly'
    df.at[0, 'associated_docker_container'] = ''
    df.at[0, 'version'] = workflow_version

    df.at[1, 'software'] = 'seqyclean'
    df.at[1, 'associated_docker_container'] = 'staphb/seqyclean'
    df.at[1, 'version'] = seqyclean_version

    df.at[2, 'software'] = 'fastqc'
    df.at[2, 'associated_docker_container'] = 'staphb/fastqc'
    df.at[2, 'version'] = fastqc_version

    df.at[3, 'software'] = 'bwa'
    df.at[3, 'associated_docker_container'] = 'broadinstitute/viral-core'
    df.at[3, 'version'] = bwa_version

    df.at[4, 'software'] = 'samtools'
    df.at[4, 'associated_docker_container'] = 'broadinstitute/viral-core'
    df.at[4, 'version'] = samtools_version_broadinstitute

    df.at[5, 'software'] = 'ivar'
    df.at[5, 'associated_docker_container'] = 'andersenlabapps/ivar'
    df.at[5, 'version'] = ivar_version

    df.at[6, 'software'] = 'samtools'
    df.at[6, 'associated_docker_container'] = 'andersenlabapps/ivar'
    df.at[6, 'version'] = samtools_version_andersenlabapps

    df.at[7, 'software'] = 'samtools'
    df.at[7, 'associated_docker_container'] = 'staphb/samtools'
    df.at[7, 'version'] = samtools_version_staphb

    # add project name and anaysis date
    df['project_name'] = project_name
    df['analysis_date'] = analysis_date

    col_order = ['project_name', 'analysis_date', 'software', 'associated_docker_container', 'version']
    df = df[col_order]
    
    outfile = f'version_capture_illumina_pe_asembly_{project_name}_{workflow_version}.csv'
    df.to_csv(outfile, index = False)
