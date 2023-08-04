# SC2_illumina_pe_assembly Workflow

## Overview
This workflow was developed for the assembly of Illumina 150 bp paired-end read data using the Illumina Nextera XT library prep protocol. The workflow accepts "sample" as the root entity type. The following steps are preformed during the workflow:
1. Fastq files are quality filtered and trimed using Seqyclean 
  - Seqyclean parameters include a minimum read length set to 70 bp and quality trimming set to a minimum Phred quality score of 30.
2. FastQC is run on both the raw (pre filtering and trimming) and cleaned (post filtering and trimming) reads
3. Reads are aligned to the reference genome using bwa and then sorted by coordinates using Samtools.
4. Primer regions of reads are trimed using iVar trim and then the trimmed reads are sorted by coordinates using Samtools.
5. Variants are called relative to the reference genome using iVar variants. 
  - iVar variants parameters include a minimum quality score set to 20, a minimum variant base frequency set to 0.6 and a minimum read depth set to 10.
6. A consensus genome is generated using iVar consensus.
  - iVar consensus parameters include a minimum quality score set to 20, a minimum variant base frequency set to 0.6 and a minimum read depth set to 10.
7. Post assembly statistics, including mean seequencing depth and percent genome coverage, are calculated from the bam file using Samtools flagstat, stats, and coverage and a custom python script (calc_percent_coverage.py). 
8. The fasta header of the consensus sequence is renamed in the GISAID-acceptable format: CO-CDPHE-{sample_name}.
9. Workflow and software versions are recorded.
10. Intermediate files, consensus sequence fastas and post assembly summary stats files are transfered to a local GCP bucket for storage. 


![SC2_illumina_pe_assembly.wdl workflow diagram](./workflow_diagrams/SC2_illumina_pe_assembly.png "SC2_illumina_pe_assembly.wdl workflow diagram")


<br/>

### Inputs
1. Terra data table.

   The terra data table must include the following columns as listed below. Note that optional columns are not neccessary for the assembly workflow but must be present for the SC2_lineage_calling_and results.wdl and Transfer workflows described below under ``Lineage Calling Workflows`` and ``Transfer Workflows``, respecitively.

| column header | description | 
|-------------------|-----------------|
| ``entity:sample_id``| column with the list of sample names. (e.g. ``entity:covwwt-0203_id``) |
| ``fastq_1``| The google bucket path to the R1 fastq file.|
|``fastq_2``| The google bucket path to the R2 fastq file.|
|``out_dir``| User defined google bucket for where the files will be transfered during the transfer workflows. |
|``workbook_path``| (optional; required for lineage calling workflow) | 
|``project_name``| (optional; requried for lineage calling workflow) |

<br/>

2. Terra Workspace Data.

  See [Getting set up](#getting-set-up) above. 

<br/>

3. Setting up the workflow inputs

  For setting up the worklfow inputs, use the ``SC2_illumina_pe_assembly-input.json`` in the ``workflow_inputs`` directory.

  |workflow variable| attribute (input syntax into workflow) |
  |------------|-----------|
  |``adapters_and_contaminants``|workspace.adapters_and_contaminants_fa|
  |``calc_percent_coverage_py``|workspace.covid_calc_percent_coverage_py|
  |``covid_genome``|workspace.covid_genome_fa|
  |``covid_gff``|workspace.covid_genome_gff|
  |``fastq_1``|this.fastq_1|
  |``fastq_2``|this.fastq_2|
  |``primer_bed``|workspace.artic_v4-1_bed|
  |``s_gene_amplicons``|workspace.artic_v4-1_s_gene_amplicons|
  |``sample_name``|this.{entity_name}_id|

<br/>

### Outputs

| WDL task name | software/program | variable name | description |
|---------------|------------------|---------------|-------------|
|seqyclean| seqyclean| ``filtered_reads_1``| file|
|seqyclean| seqyclean| ``filtered_reads_2``| file|
|seqyclean| seqyclean| ``seqyclean_summary``| file|
|fastqc as fastqc_raw| fastqc| ``fastqc_raw1_html``| file|
|fastqc as fastqc_raw| fastqc| ``fastqc_raw1_zip``| file|
|fastqc as fastqc_raw| fastqc| ``fastqc_raw2_html``| file|
|fastqc as fastqc_raw| fastqc| ``fastqc_raw2_zip``| file|
|fastqc as fastqc_cleaned| fastqc| ``fastqc_clean1_html``| file|
|fastqc as fastqc_cleaned| fastqc| ``fastqc_clean1_zip``| file|
|fastqc as fastqc_cleaned| fastqc|``fastqc_clean2_html``| file|
|fastqc as fastqc_cleaned| fastqc| ``fastqc_clean2_zip``| file|
|align_reads|bwa, samtools| ``out_bam``| file|
|align_reads|bwa, samtools| ``out_bamindex``| file|
|align_reads|bwa, samtools|``assembler_version``| string recording the version for bwa, this information is used later for submitting to public repositories.|
|ivar trim |ivar trim, samtools| ``trim_bam`` | file|
|ivar trim |ivar trim, samtools| ``trimsort_bam`` | file|
|ivar trim |ivar trim, samtools| ``trimsort_bamindex`` | file|
|ivar variants| ivar variants| ``variants``| vcf file formatted as a tsv|
|ivar consensus| ivar consnesus| ``consensus``| fasta file of conensus genome, Ns are called in places with less than 10 bp read depth. |
|bam_stats|samtools flagstat, stats, percent_coverage | ``flagstat_out``| file|
|bam_stats|samtools flagstat, stats, percent_coverage | ``stats_out``| file|
|bam_stats|samtools flagstat, stats, percent_coverage |  ``covhist_out``| file|
|bam_stats|samtools flagstat, stats, percent_coverage |  ``cov_out``| file|
|bam_stats|samtools flagstat, stats, percent_coverage | ``cov_s_gene_amplcions_out``| file|
|bam_stats|samtools flagstat, stats, percent_coverage | ``cov_s_gene_out``|file|
|rename_fasta| N/A | ``renamed_consensus``|fasta file; consesnus genome sequence with the fasta header renamed to be CO-CDPHE-{sample_name}|
|calc_percent_cvg|calc_percent_coverage.py| ``percent_cvg_csv``|csv file, see calc_percent_cvg.py script readme for details found in the ./python_scripts directory of this repository.|


<br/>