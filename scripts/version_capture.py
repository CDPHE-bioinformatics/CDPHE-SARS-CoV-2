import argparse
import json
import pandas as pd
import sys

def get_options(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('--versions_json', help='JSON file with an array of version_info objects (keys should be "software", "docker", and "version")', required=True)
    parser.add_argument('--workflow_name', help='workflow name (e.g. SC2_ont_assembly)', required=True)
    parser.add_argument('--workflow_version', help='workflow version with dashes (e.g. v2-2-0)', required=True)
    parser.add_argument('--project_name', help='project name of batch being analyzed (e.g. cov_2022_grid)', required=True)
    parser.add_argument('--analysis_date', help='date to add to the analysis_date column in the output CSV', required=True)
    parser.add_argument('--sample_name', help='(use for sample-level workflows) sample_name for column and file name')
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
    if options.sample_name:
        df.insert(0, 'sample_name', options.sample_name)
    df = df.rename(columns={'docker': 'associated_docker_container'})  # backwards-compat naming
    return df

if __name__ == '__main__':
    options = get_options()
    df = create_version_df(options)
    if options.sample_name:
        outfile =  f"version_capture_{options.sample_name}_{options.workflow_name}_{options.project_name}_{options.workflow_version}.csv"
    else:
        outfile = f"version_capture_{options.workflow_name}_{options.project_name}_{options.workflow_version}.csv"
    df.to_csv(outfile, index=False)
