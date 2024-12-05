# Setup

## Input Workflow from Dockstore
To use the workflow on the Terra platform, first you will need to import the workflow from Dockstore. All workflows can be found under our dockstore organization called CDPHE-bioinformatics.
1. Go to dockstore (https://dockstore.org/).
2. Along the top search bar click on Organizations and search for "CDPHE".
3. Select the workflow.
4. On the right hand side of the workflow description, select "Launch with Terra".
5. Select the Destination workspace and select "Import". 
6. The workflow will now be displayed as a card under your workflows tab in your Terra workspace. 

## Workspace Data
Prior to running any of the workflows, you must set up the Terra workspace data with the correct reference files and custom python scripts. The reference files can be found in this repository in the `data/workspace_data` directory. Python scripts can be found in the `scripts` directory. Workspace variables are named using the following format `{organism}_{description}_{file_type}`, except for the primer bed files which are named as `{description}_{file_type}`. Reference files and python scripts should be copied from this repo into a GCP bucket. The GCP bucket path to the file will serve as the "value" when adding data to the terra workspace data table. 

To add data to the terra workspace data:
1. Navigate to the Data tab in your Terra workspace.
2. In the left hand list of data tables, under "Other Data" select "Workspace Data".
3. Click on the "+" button in the lower right hand corner of the workspace data table. 
4. Fill in the "Key" column with the workspace variable name, the "Value" column with GCP bucket path to the file and the "Description" column with a brief description if desired. 
5. Once complete hit the check mark to the right.

| workspace variable name | workflow|  file name | description | 
|-------------------|-----------------|--------------------|-----------------|
| `adapters_and_contaminants_fa` | `SC2_illumina_pe_assembly` | Adapters_plus_PhiX_174.fasta | adapters sequences and containment sequences removed during fastq cleaning and filtering using SeqyClean. Thanks to Erin Young at Utah Public Health Laboratory for providing this file!  |
| `artic_v3_bed` |  `SC2_illumina_pe_assembly`, `SC2_ont_assembly` |artic_V3_nCoV-2019.primer.bed| primer bed file for the Artic V3 tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| `artic_v4_bed` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` |artic_V4_nCoV-2019.primer.bed | primer bed file for the Artic V4 tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| `artic_v4-1_bed` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` |artic_V4-1_nCoV-2019.primer.bed | primer bed file for the Artic V4.1 tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| `artic_v4-1_s_gene_amplicons` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` |artic_v4_1_s_gene_amplicons.tsv|coordinate positions of S gene amplicons using the artic V4.1 primers|
| `artic_v4-1_s_gene_primer_bed` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` | S_gene_V4-1_nCoV-2021.primer.bed|primer sequences and coordinate positions of primer binding region of Artic v4.1 primers|
| `artic_v5-3-2_bed` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` | artic_v5-3-2_nCoV-2023.primer.bed|primer bed file for the Artic V5.3.2 tiled amplicon primer set. |
| `artic_v5-3-2_s_gene_amplicons` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` |artic_v5-3-2_s_gene_amplicons.tsv| coordinate positions of S gene amplicons using the artic V5.3.2 primers|
| `artic_v5-3-2_s_gene_primer_bed` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` | S_gene_V5-3-2_nCoV-2021.primer.bed| primer sequences and coordinate positions of primer binding region of Artic v5.3.2 primers|
| `midnight_bed` | `SC2_ont_assembly` | Midnight_Primers_SARS-CoV-2.scheme.bed | primer bed file for the Midnight tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| `covid_genome_fa` | `SC2_illumina_pe_assembly`,  `SC2_ont_assembly`, `SC2_wastewater_variant_calling` |MN908947-2_reference.fasta | SARS-CoV-2 whole genome reference sequence in fasta format (we use NCBI genbank ID MN908947.3) |
| `covid_genome_gff` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly`, `SC2_wastewater_variant_calling` | NC_045512-2_reference.gff| whole genome reference sequence annotation file in gff format (we use NCBI genbank ID MN908947.3) |
| `covid_genome_gff_mutations` | `SC2_novel_mutations` | novel_mutations_gff.tsv | tsv formatted version of `covid_genome_gff` for use with `novel_mutations_append_py` |
| `covid_voc_annotations_tsv` | `SC2_wastewater_variant_calling workflow` | SC2_voc_annotations_20220711.tsv | For wastewater only. List of amino acid (AA) substitutions and lineages containing those AA substitutions; for a lineage to be associated with a given AA substitution, 90% of publicly available sequences must contain the AA substitution (the 90% cutoff was determined using outbreak.info) |
| `covid_voc_bed_tsv` | `SC2_wastewater_variant_calling workflow` | SC2_voc_mutations_20220711.tsv |  For wastewater only. List of nucleotide genome positions in relation to the MN908947.3 reference genome of know mutations |
| `covid_calc_per_cov_py` | `SC2_illumina_pe_assembly`, `SC2_ont_assembly` |calc_percent_coverage.py | see detailed description in the readme file found in `./python_scripts/` repo directory|
| `covid_nextclade_json_parser_py` | `SC2_lineage_calling_and_results` | nextclade_json_parser.py | see detailed description in the readme file found in `./python_scripts/` repo directory|
| `covid_concat_results_py` | `SC2_lineage_calling_and_results` | concat_seq_metrics_and_lineages_results.py | see detailed description in the readme file found in `./python_scripts` repo directory |
| `covid_novel_mutations_append_py` | `SC2_novel_mutations` | novel_mutations_append.py | see detailed description in the readme file found in `./python_scripts/` repo directory 
| `covid_version_capture_lineage_calling_py` | `SC2_lineage_calling_and_results` | version_capture_lineage_calling_and_results.py| generates version capture output file for software versions used in the SC2_lineage_calling_and_results workflow|
| `covid_version_capture_wastewater_variant_calling_py` | `SC2_wastewater_variant_calling` | version_capture_wastewater_variant_calling.py| generates version capture output file for software versions used in the SC2_wastewater_variant_calling workflow|
| `novel_mutations_historical_full` | `SC2_novel_mutations` | novel_mutations.py | for wastewater only. See detailed description in the readme file found in `./python_scripts/` repo directory |
| `novel_mutations_historical_unique` | `SC2_novel_mutations` | novel_mutations.py | for wastewater only. See detailed description in the readme file found in `./python_scripts/` repo directory |
| `hostile_human_t2t_hla_bt2` | `SC2_illumina_pe_assembly` | human-t2t-hla.1.bt2, human-t2t-hla.2.bt2, ... | hostile read scrubbing bowtie2 default reference genome index |
| `hostile_human_t2t_hla_fa_gz` | `SC2_ont_assembly` | human-t2t-hla.fa.gz | hostile read scrubbing minimap2 default reference genome |
| `version_capture_py` | `SC2_ont_assembly`, `SC2_illumina_pe_assembly` | version_capture.py | general script for creating version capture CSV file |