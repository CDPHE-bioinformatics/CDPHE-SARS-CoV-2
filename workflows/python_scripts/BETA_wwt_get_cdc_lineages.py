#! /usr/bin/env python

version = '1.0.0'

# updating the script to handle hsn values and dropping the duplicates


# import python modules
import os
import glob
import re
import shutil
import pandas as pd 
from datetime import date
import numpy as np
import subprocess

import sys
import argparse


#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "--freyja_aggregate_file")
    parser.add_agrument("--project_name")
    parser.add_agrument("--terra_data_table_path")
    parser.add_agrument('--cdc_linage_groups_json')
    options = parser.parse_args(args)
    return options

#### MAIN ####
if __name__ == '__main__':

    options = getOptions()
    freyja_aggregated_file = options.freyja_aggregated_file
    project_name = options.project_name
    terra_data_table_path = options.terra_data_table_path
    cdc_lineage_groups_json = options.cdc_lineage_groups_json
