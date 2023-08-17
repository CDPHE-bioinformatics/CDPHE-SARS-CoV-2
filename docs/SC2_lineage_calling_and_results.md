# SC2_lineage_calling_and_results Workflow

## Overview
This workflow should be run following assembly with one of the SC2 reference based assembly workflows. The workflow accepts "sample_set" as the root entity type and uses the same data table used for the assembly workflows. All maintained assembly workflows are compatible with this workflow. Breifly, the following steps are preformed during the workflow:
1. consensus sequences are conctenated into a single fasta file. 
2. Lineages are assigned using Pangolin
3. Nextclade is run. 
4. The Nextclade json ouput files are parsed using the ``nextclade_json_parser.py`` script which pulls out clade and nucleotide and amino acid changes and converts the output to a tabular format.
5. Sequenicng assembly metrics (e.g. percent coverage, mean depth), lineage and clade information, and sequence metadata (e.g. plate name, sample well location) are concatenated into a single csv file. This file serves as the summary results file. 
6. A csv file with sequencing assembly metrics and lineage information is generated that can be used to parse sequencing data into our LIMS.
7. Intermediate and summary files are transfered to a user defined google bucket. 
8. Workflow and software versions are recorded.

<br/>

## Inputs and Setup
<br/>

### Terra data table

The terra data table must include the following columns as listed below. Many of the columns are generated and autopopulated with files from the reference based assembly workflows.  

| column header | description | 
|-------------------|-----------------|
| ``entity:{sample}_id``| Column with the list of sample names, where ``{sample}`` can be replaced with a descriptor. For example the header could be: ``entity:covwwt-0203_id``. Each sample name within this column must be unique. Must be the first column. |
| ``cov_out`` | generated from assembly workflow; txt file with mean depth calcuations produced by samtools.  |
| ``out_dir`` | user defined google bucket for where the files will be transfered during the transfer workflows.|
|``percent_cvg_csv`` | generated from assembly workflow; csv file with percent coverage calcualtions |
|``project_name``| The name of the sequencing project. |
|``renamed_consensus``| generated from assembly workflow; consensus sequences as fasta file.|
|``terra_data_table_path``| The GCP bucket path to the terra data table location. |

<br/>

### Terra Workspace Data.

see the SC2_workflows_overview readme for which workspace data elements are required and how to add the correct files to your workspace data. 

<br/>

### Setting up the workflow inputs
Navigate to the workflow launch page within your Terra workspace. On the launch page:
1. Select the correct version of the workflow. Be sure to select the latest release. Other branches or versions are not garanteed to run successfully and are considered develop versions.
2. Choose the data table, being sure to select root entity type as ``sample_set``. Then select the samples from the data table you want to include in set for the anlaysis. 
3. Set up the workflow inputs as follows in the table below (you can also use the available example json input file under the data directory). For the attributes, the ""this.sample{terra_datatable_name}s." syntax refers Terra to pull the variable from the terra datatable as used for sample sets. These variables were either in the original terra datatable as inputs for the assembly workflow (see referece based assembly workflow inputs sections above for more details) or added as outputs during the assemlby workflow (see reference based assembly workflow outputs sections for more details): 

    |Workflow Variable| Type| Attribute (input syntax into workflow) |
    |------------|-----------|---------------------------------------|
    |``cdc_lineage_groups_json``| File | workspace.covid_cdc_lineage_groups_json |
    |``concat_seq_results_py`` | File | workspace.covid_concat_results_py |
    |``cov_out``| Array[File?] this.sample{terra_datatable_name}s.cov_out |
    |``nextclade_json_parser_py``| File | workspace.covid_nextclade_json_parser_py | 
    |``out_dir_array``| Array[String] | this.sample{terra_datatable_name}s.out_dir|
    |``percent_cvg_csv`` | Array[File?] | this.sample{terra_datatable_name}s.percent_cvg_csv |
    |``project_name_array``| Array[String] | this.sample{terra_datatable_name}s.project_name |
    |``renamed_consensus`` | Array[File?] this.sample{terra_datatable_name}s.renamed_consesnus|
    |``sample_name`` | Array[String] | this.sample{terra_datatable_name}s.sample{terra_datatable_name}_id |
    |``terra_data_table_path_array`` | Array[File] | this.sample{terra_datatable_name}s.terra_data_table_path |
    |``version_capture_lineage_calling_and_results_py`` | File | workspace.covid_version_capture_lineage_calling_and_results_py | 


4. Select ``Use Defaults`` for outputs. 

<br/>

## Outputs

The table below lists the following outputs that are generated during the workflow and can be accessed through the terra data table once the workflow finishes. The WDL task name indicates the task in the WDL where the output was generated. The software/program indicates the software/program taht was used in the WDL task to generate the output. The variable name will be the name of the column header generated in the terra data table where the result value or file will be stored/located. The description will state whether the file is transferred to a local GCP storage bucket. If it is not transferred the output can only be accessed via Terra's backend GCP storage buckets and/or the terra data table. 

| WDL task name | software/program | variable name | description |
|---------------|------------------|---------------|-------------|
| concatenate | bash | ``cat_fastas`` | fasta file; transfered to local GCP storage bucket|
| pangolin | pangolin | ``pangolin_version`` | string |
| pangolin | pangolin | ``pangolin_lineage`` | csv file; transfered to local GCP storage bucket | 
| nextclade | nextclade | ``nextclade_version`` | string | 
| nextcalde | nextclade | ``nextclade_json`` | json file; transfered to local GCP storage bucket| 
| nextclade | nextclade | `` nextclade_csv`` | csv file; transfered to local GCP storage bucket |
| parse_nextclade | nextclade_json_parser.py | ``nextclade_variants_csv`` | csv file; transfered to local GCP storage bucket|
|parse_nextclade | nextclade_json_parser.py | ``nexclade_clades_csv`` | csv file; transfered to local GCP storage bucket |
| results_table | concat_seq_results.py | ``sequencing_results_csv`` | csv file; transfered to local GCP storage bucket|
| results_table | concat_seq_results.py | ``wgs_horizon_report_csv`` | csv file; transfered to local GCP storage bucket|
|create_version_capture_file| version_capture_lineage_calling_and_results.py| ``version_capture_lineage_calling_and_results``| csv file; transfered to local GCP storage bucket|
|transfer| gsutil | ``transfer_date_lineage_calling``| string|


<br/>