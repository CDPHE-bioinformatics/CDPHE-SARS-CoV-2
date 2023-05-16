import argparse
import sys
import pandas as pd
import numpy as np


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    
    parser.add_argument("--project_name",  help= "")
    parser.add_argument('--guppy_version', help = '')
    parser.add_argument('--medaka_version', help = '')
    options = parser.parse_args(args)
    return options


if __name__ == '__main__':
    
    options = getOptions()
    project_name = options.project_name
    guppy_version = options.guppy_version
    medaka_version = options.medaka_version

    df = pd.DataFrame()

    df['project_name'] = project_name
    df['guppy_version'] = guppy_version
    df['medaka_version'] = medaka_version
    
    outfile = f'{project_name}_assembly_software.tsv'
    df.to_csv(outfile, sep = '\t', index = False)