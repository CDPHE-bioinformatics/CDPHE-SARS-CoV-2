import argparse
import sys
import pandas as pd
import numpy as np


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument("--project_name",  help= "")
    parser.add_argument("--guppy_version",  help= "")
    parser.add_argument("--artic_version",  help= "")
    parser.add_argument("--medaka_version")
    parser.add_argument("--samtools_version",  help= "")
    parser.add_argument("--pyScaf_version",  help= "")
    parser.add_argument("--bcftools_version",  help= "")
    parser.add_argument("--analysis_date",  help= "")
    parser.add_argument("--workflow_version",  help= "")
    options = parser.parse_args(args)
    return options


if __name__ == '__main__':
    
    options = getOptions()
    project_name = options.project_name
    guppy_version = options.guppy_version
    artic_version = options.artic_version
    medaka_version = options.medaka_version
    samtools_version = options.samtools_version
    pyScaf_version = options.pyScaf_version
    bcftools_version = options.bcftools_version
    analysis_date = options.analysis_date
    workflow_version = options.workflow_version

    df = pd.DataFrame()

    # begin to fill in table
    df.at[0, 'software'] = 'SC2_ont_assembly'
    df.at[0, 'associated_docker_container'] = ''
    df.at[0, 'version'] = workflow_version

    df.at[1, 'software'] = 'guppy'
    df.at[1, 'associated_docker_container'] = 'genomicpariscentre/guppy'
    df.at[1, 'version'] = guppy_version

    df.at[2, 'software'] = 'artic'
    df.at[2, 'associated_docker_container'] = 'quay.io/staphb/artic-ncov2019'
    df.at[2, 'version'] = artic_version

    df.at[3, 'software'] = 'medaka'
    df.at[3, 'associated_docker_container'] = 'quay.io/staphb/artic-ncov2019'
    df.at[3, 'version'] = medaka_version

    df.at[4, 'software'] = 'samtools'
    df.at[4, 'associated_docker_container'] = 'staphb/samtools'
    df.at[4, 'version'] = samtools_version

    df.at[5, 'software'] = 'pyScaf'
    df.at[5, 'associated_docker_container'] = 'chrishah/pyscaf-docker'
    df.at[5, 'version'] = pyScaf_version

    df.at[6, 'software'] = 'bcftools'
    df.at[6, 'associated_docker_container'] = 'staphb/bcftools'
    df.at[6, 'version'] = bcftools_version

  
    # add project name and anaysis date
    df['project_name'] = project_name
    df['analysis_date'] = analysis_date

    col_order = ['project_name', 'analysis_date', 'software', 'associated_docker_container', 'version']
    df = df[col_order]
    
    outfile = f'version_capture_ont_asembly_{project_name}_{workflow_version}.csv'
    df.to_csv(outfile, index = False)