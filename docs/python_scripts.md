# SARS-CoV-2 Python Scripts

## Table of Contents

1. [Introduction](#introduction)
2. [calc_percent_coverage-py](#calc-percent-coverage-py)
3. [nextclade_json_parser.py](#nextclade-json-parser-py)
4. [concat-seq-metrics-and-lineage-results-py](#concat-seq-metrics-and-lineage-results-py)
5. [novel_mutations.py](novel-mutations-py)
6. [details about working with sample sets](#details-about-working-with-sample-sets)

## Introduction

This repo contains four custom python scripts called in our SC2 wdl workflows. These scripts should be stored in a google bucket and linked to your workspace data within your Terra.bio workspace, so they can be used as inputs within each respective workflow. Listed below are the four scripts and their associated workflows. 

- `calc_percent_coverage.py` is called in `SC2_ont_assembly.wdl`, `SC2_illumina_pe_assembly.wdl`, and `SC2_illumina_se_assembly.wdl`.
- `nextclade_json_parser.py` is called in `SC2_lineage_calling_and_results.wdl`.
- `concat_seq_metrics_and_lineage_results.py` is called in `SC2_lineage_calling_and_results.wdl`.
- `novel_mutations.py` is called in `SC2_novel_mutations.wdl`.


## calc-percent-coverage-py

### Overview

This script is called in the `SC2_illumina_pe_assembly.wdl`, `SC2_illumina_se_assembly.wdl`, and `SC2_nanopore_assembly.wdl` WDL workflows. These workflows act on individual samples (as opposed to sample sets); therefore this script also works on individual samples. The script reads in a consensus genome as a fasta file and calculates the percent coverage of the sample consensus genome. Specifically, the percent coverage is calculated as:

```text
percent_coverage = frac{number_non_ambiguous_bases}{29903}, * 100
```

where the number_non_ambiguous_bases is the number of base pair calls not including Ns in the sample consensus sequence and 29903 is the number of base pairs in the reference genome (NC_045512).

### Inputs

| flag            | description                                                                                                                                |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------|
| `--sample_name` | sample name as recorded in the `entity:sample_ID` column in the terra data table                                                           |
| `--fasta_file`  | renamed consensus sequence fasta file; generated during the ivar_consensus (illumina data) or medaka (ont data) and rename fasta tasks of the assembly workflows |

### Outputs

The script also records the number of aligned bases, the number of ambiguous bases (Ns), and the number of non-ambiguous bases (A,G,C,T). The output is a csv file called `sample_name}_consensus_cvg_stats.csv` and has the following column headers:

| column header                | description                                                                                         |
| ---------------------------- | --------------------------------------------------------------------------------------------------- |
| `sample_name`                | sample name as recorded in the `entity:sample_ID` column in the terra data table                    |
| `number_aligned_bases`       | the total length of the consensus genome (including Ns)                                             |
| `number_N_bases`             | the number of ambiguous (N) bases in the consensus genome                                           |
| `number_non_ambiguous_bases` | the number of non-ambiguous (A, C, G, T) bases in the consensus genome                              |
| `percent_coverage`           | percent coverage of the reference genome represented in the consensus sequence as calculated above  |

There is an example output in the example data directory within this repo.

## nextclade-json-parser-py

### Overview

This script is called in the `SC2_lineage_calling_and_results` WDL workflow. This workflow and script act on sample sets. This script is called in the `parse_nextclade` task within the workflow, which can be seen in the `SC2_lineage_calling_and_results.wdl` workflow diagram in the docs/img directory. Briefly, the workflow concatenates all consensus sequences of the samples in the sample set into a single fasta file (`concatenate` task). The concatenated fasta file is run through nextclade which generates a `nextclade.json` file (`nextclade` task). Within the `nextclade.json` file is data for each sample consensus sequence, including the nextclade clade designation, AA substitutions, deletions, insertions, etc. Generally, this script reads in the `nextclade.json` file, parses the json file to extract the data of interest, formats the data into a table, and saves it as a csv file.

### Inputs

| flag               | description                                                                |
| ------------------ | -------------------------------------------------------------------------- |
| `--nextclade_json` | the nextclade.json file generated in the `nextclade` task of the workflow |
| `--project_name`   | the project name of the sequencing run                                     |

### Outputs

There are two outputs from this script, each accomplished from a separate function within the script. Example outputs can be found in the example data directory within this repo. These functions are:

1. `extract_variant_list()` function : This function generates a summary of the AA substitutions, insertions, and deletions for each sample within the `nextclade.json` file. The output is a csv file called `{project_name}_nextclade_variant_summary.csv` which is one of the files that is transferred to the google bucket as outputs of the workflow. The data is formatted such that each row corresponds to a either an AA substitution, insertion, or deletion, such that each consensus sequence can have more than one row of data. The csv file has the following column headers:

    | column header    | description                                                                                                                           |
    | ---------------- | --------------------------------------------------------------------------------------------------------------------------------------|
    | `sample_name`    | the sample name as listed in the fasta header (therefore there will be a "CO-CDPHE-" prefix added to the sample_name as listed in the entity column in the terra data table)                       |
    | `variant_name`   | the full name of the variant formatted as {gene}\_{refAA}{AApos}{altAA} (e.g. S_L452Q or S_L24del). For insertions: the gene is not listed, the refAA is defined as "ins", the AA position is instead the nucleotide position in the genome, and the altAA is listed as the string of nucleotides (e.g. "\_ins1027T", which would be interpreted as an insertion of a T nucleotide occurred at genome position 1027) |
    | `gene`           | the gene where the AA substitution, deletion or insertion occurs (e.g. N, S, ORF1a, M etc.). For insertions: null.           |
    | `codon_position` | the codon position (or protein position) within the gene where the AA substitution, deletion or insertion occurred. For insertions: nucleotide genome position                              |
    | `refAA`          | the reference AA at the position where the AA substitution, deletion, or insertion occurred. For insertions the refAA is listed as "ins". |
    | `altAA`          | he AA in the consensus sequence at the position where the AA substitution, deletion, or insertion occurred. For insertions: the string of nucleotide base pairs that were inserted       |
    | `start_nuc_pos`  | the starting nucleotide position within the genome where the AA substitution, deletion, or insertion occurred                            |
    | `end_nuc_pos`    | the ending nucleotide position within the genome where the AA substitution, deletion, or insertion occurred (for a single AA substitution the start_ncu_pos and end_nuc_pos will be a difference of 3) |

2. `get_nextclade()` function: This function generates a summary of the nextclade designation, total nucleotide and AA substitutions, total nucleotide and AA deletions, and total nucleotide insertions. The output file is called `{project_name}_nextclade_results.csv` and is used as input for the `concat_seq_metrics_and_lineage_results.py` called in the `results_table` task in the workflow. The output file has the following column headers:

    | column header                 | description                                        |
    | ----------------------------- | ---------------------------------------------------|
    | `sample_name`                 | the sample name as listed in the fasta header (therefore there will be a "CO-CDPHE-" prefix added to the sample_name as listed in the entity column in the terra data table) |
    | `nextclade`                   | nextclade clade designation (e.g. 22C (Omicron))   |
    | `total_nucleotide_mutations`  | number of SNPs in consensus sequence               |
    | `total_nucleotide_deletions`  | number of deletions in the consensus sequence      |
    | `total_nucleotide_insertions` | number of insertions in the consensus sequence     |
    | `total_AA_substitutions`      | number of AA substitutions in the consensus sequence |
    | `total_AA_deletions`          | number of AA deletions in the consensus sequence   |

## concat-seq-metrics-and-lineage-results-py

### Overview

This script is called in the `SC2_lineage_calling_and_results` WDL workflow. This workflow and script act on sample sets. This script is called in the `results_table` task within the workflow which can be seen in the `SC2_lineage_calling_and_results.wdl` workflow diagram in the docs/img directory. Generally, this script pulls together a bunch of metadata and data regarding the consensus sequence and outputs the data to a csv file.  

### Inputs

The script takes the following inputs:

| flag                       | description                                                                                                                     |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------|
| `--sample_name`            | an array of the sample_name variables for the set written to a txt file using the wdl function `write_lines`     |
| `--workbook_path`          | the gcp file path to the workbook. The workbook can include any columns you'd like but must include the following columns at minimum: `hsn`, `sample_name`, `project_name`, `plate_name`, `run_name`. These columns can be left blank if needed. |
| `--cov_out_files`         | the list of the `cov_out` variable (column in the terra data table) for the workflow written to a text file. This variable is a file path to a file with the bam stats generated in the `SC2_ont_assembly.wdl`, `SC2_illumina_se_assembly.wdl` or `SC2_illumina_pe_assembly.wdl` from the bam stats task. |
| `--percent_cvg_files`      | the list of the `percent_cvg_csv` variable (column in the terra data table) for the workflow written to a text file. This variable is a file path to a file with the bam stats generated in the `SC2_ont_assembly.wdl`, `SC2_illumina_se_assembly.wdl` or `SC2_illumina_pe_assembly.wdl` workflows from the `calc_percent_coverage.py` script called during the `calc_percent_cvg` task. |
| `--assembler_version`      | the assembler_version variable (column in the terra data table) for the workflow . This is written to the terra data table during the `SC2_ont_assembly.wdl`, `SC2_illumina_se_assembly.wdl` or `SC2_illumina_pe_assembly.wdl` workflows. |
| `--pangolin_lineage_csv`   | this is the lineage report csv file generated from pangolin during the `pangolin` task                                          |
| `--nextclade_clades_csv`   | this is the `{project_name}_nextclade_results.csv` file generated from the `nextclade_json_parser.py` script during the `parse_nextclade` task. |
| `--nextclade_variants_csv` | this is the `{project_name}_nextclade_variant_summary.csv` file generated from the `nextclade_json_parser.py` script during the `parse_nextclade` task. |
| `--nextclade_version`      | this is the nextclade version which is defined as output during the `nextclade` task.                                           |
| `--project_name`           | project_name from column in terra data table                                                                                    |

### Outputs

There are three outputs from this script. Example outputs can be found in the example data directory within this repo.

1. `{project_name}_sequencing_results.csv`: summary of sequencing metrics for all samples within the sample set. Below is a table of the column headers and their description. Currently all headers from the sequencing workbook that get carried over to the terra data table are included in this output file. The list below includes some, but not all, of the columns in the file.

    | column header name               | description                                                                                                           |
    | -------------------------------- | ----------------------------------------------------------------------------------------------------------------------|
    | `sample_name`                    | sample name                                                                                                           |
    | `hsn`                            | hsn (horizon serial number)                                                                                           |
    | `primer_set`                     | name of primer set used for tiled amplicon sequencing (Artic V3, Artic V4, Artic V4.1, Midnight or COVIDSeqV3)        |
    | `percent_coverage`               | percent coverage; the total proportion of the genome that is covered, not including regions where an N is called for a base call |
    | `nextclade`                      | the nextclade clade assignment                                                                                        |
    | `pangolin_lineage`               | the pangolin lineage assignment                                                                                       |
    | `pangolin_expanded_lineage`      | the expanded pangolin lineage                                                                                         |
    | `assembler_version`              | assembler software version (either bwa or minimap, depending on assembly workflow used)                               |
    | `spike_mutations`                | list of spike mutations in the spike gene sequence that correspond to key spike mutations identified in the sample consensus sequence (this column was created prior to VOCs and includes spike mutations we were watching and has not been updated since) |
    | `total_nucleotide_mutations`     | number of SNPs in the consensus sequence genome                                                                       |
    | `total_AA substitutions`         | number of amino acid substitutions in the consensus sequence genome                                                   |
    | `total_AA_deletions`             | number of deletions in the consensus sequence genome                                                                  |
    | `mean_depth`                     | average number of reads per nucleotide site in the consensus sequence genome                                      |
    | `number_aligned_bases`           | total number of bases aligned to the reference genome (including Ns; so pretty much tells you how much was cut off the ends of the genome) |
    | `number_non_ambiguous_bases`     | total number of non-N bases in the consensus genome sequence                                                          |
    | `number_seqs_in_fasta`           | total number of sequences in the consensus fasta - should always be 1                                                 |
    | `total_nucleotide_deletions`     | number of deletions in the consensus genome sequence                                                                  |
    | `total_nucleotide_insertions`    | number of insertions in the consensus genome sequence                                                                 |
    | `num_reads`                      | total sequencing reads                                                                                                |
    | `mean_base_quality`              | mean quality score across all reads                                                                                   |
    | `mean_map_quality`               | mean mapping quality score for reads mapping to reference genome sequence                                             |
    | `number_N_bases`                 | number of bases called as N in the consensus genome sequence                                                          |
    | `nextclade_version`              | nextclade version                                                                                                     |
    | `pangolin_version`               | pangolin version                                                                                                      |
    | `pangoLEARN_conflict`            | from pangolin lineage report file                                                                                     |
    | `pangolin_ambiguity_score`       | from pangolin lineage report file                                                                                     |
    | `pangolin_scorpio_call`          | from pangolin lineage report file                                                                                     |
    | `pangolin_scorpio_support`       | from pangolin lineage report file                                                                                     |
    | `pangolin_scorpio_conflict`      | from pangolin lineage report file                                                                                     |
    | `pangolin_scorpio_notes`         | from pangolin lineage report file                                                                                     |
    | `pangolin_designation_Version`   | from pangolin lineage report file                                                                                     |
    | `pangolin_scorpio_version`       | from pangolin lineage report file                                                                                     |
    | `pangolin_constellation_version` | from pangolin lineage report file                                                                                     |
    | `pangolin_is_designated`         | from pangolin lineage report file                                                                                     |
    | `pangolin_qc_status`             | from pangolin lineage report file                                                                                     |
    | `pangolin_qc_notes`              | from pangolin lineage report file                                                                                     |
    | `pangolin_note`                  | from pangolin lineage report file                                                                                     |
    | `project_name`                   | sequencing run name                                                                                                   |
    | `tech_platform`                  | sequencing platform (e.g. Illumina MiSeq, Illumina NextSeq, Oxford Nanopore GridION)                                  |
    | `read_type`                      | single or paired end                                                                                                  |
    | `fasta_header`                   | name of the fasta header for GISAID submission (e.g. CO-CDPHE-{accession_id})                                         |
    | `analysis_date`                  | date assembly workflow ran                                                                                            |

2. `{project_name}_wgs_horizon_report.csv`: for internal use, parsing sequencing results into LIMS. Below is a table of the column headers and their description.

    | column header name   | description                                         |
    | -------------------- | --------------------------------------------------- |
    | `accession_id`       | sample name                                         |
    | `percent_coverage`   | percent coverage                                    |
    | `pangolin_lineage`   | pangolin lineage                                    |
    | `pangolin_version`   | pangolin version                                    |
    | `report_to_epi`      | this column is meaningless now but must be kept     |
    | `Run_Date`           | date assembly workflow ran                          |
    | `pangoLEARN_version` | this column is also not used but must be kept       |


## novel-mutations-py

### Overview

This script is called in the `SC2_novel_mutations` WDL workflow. This workflow is slightly different than the previous workflows- it acts on sample sets, but it and the script both create "set" and "sample" level outputs, described in more detail below. This script is called in the `append_new_mutations` task within the workflow which can be seen in the `SC2_novel_mutations.wdl` workflow diagram in the docs/img directory. Generally, this script appends new mutation data to historical files and outputs several files.  

This script follows these steps:
1. For each combined_mutations file, merges with metadata file. Drops samples from sites_to_drop, "fail" (from Freyja) samples, duplicated mutations from replicates, and controls. If any mutations do not have a collection date, it will produce a file with that information and error out.
2. Creates and outputs a unique mutation file for each project. 
3. Concatenates mutations and appends them to historical_full file. Recalculates date-related columns for historical_unique file and appends any new mutations. Checks for recurrent (not seen in > 6 months) and novel mutations and creates files if any found.
4. Output all four files from step 3. **This will overwrite both historical files, so object versioning is highly recommended.**

### Inputs

The script takes the following inputs:

| flag                       | description                                                                                                                     |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------|
| `project_names`            | space-separated list of project names of the sequencing runs                                                                    |
| `combined_mutations_files` | space-separated list of combined mutations file paths generated by Freyja and then transferred to the GCP bucket in `SC2_waste_water_variant_calling.wdl`     |
| `historical_full`          | file created by the previous iteration of this script. More information below                                                   |
| `historical_unique`        | file created by the previous iteration of this script. More information below                                                   |
| `metadata`                 | metadata file for individual samples. More information below                                                                    |
| `gff`                      | a .tsv formatted version of NC_045512-2_reference.gff                                                                           |
| `today`                    | today's date in the format "yyyy-mm-dd"                                                                                         |
| `sites_to_drop`            | space-separated list of wastewater facility site ids to remove from analysis                                                    |

### Private workspace files

Three of the input files for this script are workspace data, but will be internal to your organization.

1. `metadata`: A file that contains the relevant metadata for each sample. Our version of this file is generated by a separate script that allows us to make any manual corrections needed to the metadata. It should include the following columns:

    | column header name         | description                                                                                                                     |
    | -------------------------- | --------------------------------------------------------------------------------------------------------------------------------|
    | `project_name`             | name of sequencing run associated with the sample                                                                               |
    | `sample_name`              | name of the sample                                                                                                              |
    | `sample_type`              | sample type (e.g. "sample", "control"). Normal samples should have the type "sample"                                            |
    | `sample_name_base`         | for samples run in duplicate: remove replicate identifier from sample name (e.g. AA1234-1 and AA1234-2 both become AA1234). If not run in duplicate, just make this a copy of sample_name     |
    | `collection_date`          | sample collection date                                                                                                          |
    | `site_id`                  | wastewater facility id                                                                                                          |

2. `historical_full`: A file that contains all mutations output from Freyja from all wastewater sequencing runs, less the dropped samples described above in step 1 for the script. For the first iteration of this script, you will need to create a dummy table with the following columns:

    | column header name         | description                                                                                                                     |
    | -------------------------- | --------------------------------------------------------------------------------------------------------------------------------|
    | `project_name`             | name of sequencing run associated with the sample                                                                               |
    | `sample_name`              | name of the sample                                                                                                              |
    | `position`                 | nucleotide position of the mutation (Freyja output)                                                                             |
    | `ref_nucl`                 | reference nucleotide for this position (Freyja output)                                                                          |
    | `alt_nucl`                 | alternate nucleotide for this position (Freyja output)                                                                          |
    | `ref_aa`                   | reference amino acid for this position (Freyja output)                                                                          |
    | `alt_aa`                   | reference amino acid for this position (Freyja output)                                                                          |
    | `sample_name_base`         | for samples run in duplicate: see description above                                                                             |
    | `collection_date`          | sample collection date                                                                                                          |
    | `site_id`                  | wastewater facility id                                                                                                          |

3. `historical_unique`: A file of all unique mutations in the historical_full file. It also contains calculated features for each mutation. For the first iteration of this script, you will need to create a dummy table with the following columns:

    | column header name         | description                                                                                                                     |
    | -------------------------- | --------------------------------------------------------------------------------------------------------------------------------|
    | `position`                 | nucleotide position of the mutation (Freyja output)                                                                             |
    | `ref_nucl`                 | reference nucleotide for this position (Freyja output)                                                                          |
    | `alt_nucl`                 | alternate nucleotide for this position (Freyja output)                                                                          |
    | `ref_aa`                   | reference amino acid for this position (Freyja output)                                                                          |
    | `alt_aa`                   | reference amino acid for this position (Freyja output)                                                                          |
    | `date_first_detected`      | earliest collection date for this mutation                                                                                      |
    | `date_last_detected`       | most recent collection date for this mutation                                                                                   |
    | `times_detected`           | number of times seen                                                                                                            |
    | `length_of_time_seen`      | in days: most recent date - earliest date                                                                                       |
    | `gff_feature`              | associated gff feature (Freyja output)                                                                                          |
    | `mutation_type`            | snp, insertion, or deletion                                                                                                     |
    | `id_return`                | name of gene/protein                                                                                                            |
    | `parent_id`                | parent gff id                                                                                                                   |
    | `parent`                   | parent gene name                                                                                                                |
    | `parent_start`             | parent gene start nucleotide position                                                                                           |
    | `parent_end`               | parent gene end nucleotide position                                                                                             |
    | `gene_coordinate`          | SNP: {parent}:{ref_aa}{gene_position}{alt_aa}, deletion: {parent}:del{gene_position}, insertion: {parent}:ins{gene_position}    |
    | `nuc_coordinate`           | SNP: {ref_nucl}{position}{alt_nucl}, deletion: del{position_start}{position_end}, insertion: ins{position}{alt_nucl}            |
    | `indel_length`             | length of indel, null for SNPs                                                                                                  |


### Outputs

There is one optional output for this script that will result in an error for the workflow. There is one project-level output for each project passed. There are two higher-level outputs always generated and two higher-level optional outputs that will be generated automatically if relevant. Output files listed in this order.

1. `{project_name}_missing_dates.tsv`: file of any samples not found in the metadata file that need to be added with a collection date.

2. `{project_name}_unique_mutations.tsv`: file of unique mutations found in individual sequencing run.

    | column header name         | description                                                                                                                     |
    | -------------------------- | --------------------------------------------------------------------------------------------------------------------------------|
    | `position`                 | nucleotide position of the mutation (Freyja output)                                                                             |
    | `ref_nucl`                 | reference nucleotide for this position (Freyja output)                                                                          |
    | `alt_nucl`                 | alternate nucleotide for this position (Freyja output)                                                                          |
    | `date_first_detected`      | earliest collection date for this mutation                                                                                      |
    | `date_last_detected`       | most recent collection date for this mutation                                                                                   |
    | `times_detected`           | number of times seen                                                                                                            |
    | `length_of_time_seen`      | in days: most recent date - earliest date                                                                                       |

3. `novel_mutations_historical_full.tsv`: updated historical_full file with new mutation info. Refer to above for format.

4. `novel_mutations_historical_unique.tsv`: updated historical_unique file with new mutation info. Refer to above for format.

5. `novel_mutations_{today}.tsv`: file of any novel mutations detected in this set of runs. Same format as historical_unique.

6. `recurrent_mutations_{today}.tsv`: file of any mutations detected in this set of runs that haven't been seen in longer than 6 months. Same format as historical_unique, but with different date-related columns, listed below. 
    
    | column header name                | description                                                                                                              |
    | --------------------------------- | -------------------------------------------------------------------------------------------------------------------------|
    | `date_first_detected_new`         | earliest detection within combined_mutations file for this set of sequencing runs                                        |
    | `date_last_detected_new`          | most recent detection within combined_mutations file for this set of sequencing runs                                     |
    | `times_detected_new`              | number of times detected for this set of sequencing runs                                                                 |
    | `date_first_detected_historical`  | earliest detection within historical data                                                                                |
    | `date_last_detected_historical`   | most recent detection within historical data                                                                             |
    | `times_detected_historical`       | number of times seen within historical data                                                                              |
    | `days_between`                    | in days: date_first_detected_new - date_last_detected_historical                                                         |


## Details About Working with Sample Sets

Here I describe a way to create a single summary data table output for sample sets using a wdl workflow in terra. (It's a bit clunky but seems to work). Essentially this method enables one to create python lists from the columns in the terra data table when workflows are run as a sample set. The easiest way to explain this is with an example.

So for example, in the `concat_seq_metrics_and_lineage_results.py`, there is the input flag `--plate_name_file_list`. As input for this flag I use `${write_lines(plate_name)}`. The plate_name corresponds to the plate_name column in the terra data table. Each element in column is a string. The `write_lines()` wdl function will write each element it's own line to a text file. Thus, the input into the python script is really a text file with a list of plate names. I then wrote some code that reads in the text file and generates a python list, with each line being a new element in the list. So it looks something like this:

```python
plate_name_list = []
with open(plate_name_file_list) as f:
  for line in f:
    plate_name_list.append(line.strip())
```

In some cases, instead of string variables being stored in a column within the terra data table, file paths are stored (ie. the data type is a `File` or `Array[File]` in the case of sample sets). For example, in the `concat_seq_metrics_and_lineage_results.py`, there is the input flag `--percent_cvg_file_list`. As input for this flag I use `${write_lines(percent_cvg_csv_non_empty)}`, where the percent_cvg_csv_non_empty variable corresponds to the column percent_cvg_csv. (note the non-empty part just means that the variable may be empty for some samples in the terra data table. To set the variable as input at the beginning of the wdl I use: `Array[File?] percent_cvg_csv`). Similar to above, the script will create a list of file paths from the text file. The script can then loop through the list of file paths, open each file, extract the data from that file, and store it in a list or other data frame to be written out.

## Changelog

- updated 2023-12-12 to add novel mutations script
- updated 2023-03-09 to sync with universal naming switch over
