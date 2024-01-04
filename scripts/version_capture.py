import argparse
import json
import pandas as pd
import sys

def get_options(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('--versions_json')
    parser.add_argument('--workflow_name')
    parser.add_argument('--workflow_version')
    parser.add_argument('--project_name')
    parser.add_argument('--analysis_date')
    options = parser.parse_args(args)
    return options

def create_version_df(options):
    with open(options.versions_json) as f:
        records = json.load(f)['versions']
    
    workflow_record = {
        'software': options.workflow_name,
        'docker': '',  # workflow itself does not have a docker image
        'version': options.workflow_version
    }

    records.insert(0, workflow_record)
    df = pd.DataFrame(records)
    df.insert(0, 'project_name', options.project_name)
    df.insert(1, 'analysis_date', options.analysis_date)
    return df

if __name__ == '__main__':
    options = get_options()
    df = create_version_df(options)
    outfile = f"version_capture_{options.workflow_name}_{options.project_name}_{options.workflow_version}.csv"
    df.to_csv(outfile, index=False)
