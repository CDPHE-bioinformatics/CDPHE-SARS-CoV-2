# CDPHE_SARS-CoV-2 Workflows

## Disclaimer
**Next generation sequencing and bioinformatic and genomic analysis at the Coloardo Department of Public Health and Environment (CDPHE) is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.**

<br/>

## Overview

The following documentation describes the Colorado Department of Public Health and Environments's workflows for the assembly and analysis of whole genome sequencing data of SARS-CoV-2 on GCP's Terra.bio platform. Workflows are written in WDL and can be imported into a Terra.bio workspace through dockstore (see Setup section below; https://dockstore.org/).

Our SARS-CoV-2 whole genome reference-based assembly workflows are highly adaptable and facilitate the assembly and analysis of tiled amplicon based sequencing data of SARS-CoV-2. The workflows can accomodate various amplicon primer schemes including Artic V3, Artic V4, Artic V4.1, Artic V5.3.2 and Midnight, as well as diffent sequencing technology platforms including both Illumina and Oxford Nanopore Technology (ONT). 


<br/>

## Process

<br/>

### Clinical SC2 sequence assembly and lineage calling
Sequence assembly and lineage calling for clincial SC2 samples requries two workflows. We first run either the ``SC2_illumina_pe_assembly`` or ``SC2_ont_assembly`` workflows which perform quality control, trimming, and filtering of raw reads followed by reference-guided whole genome assembly and finally intermediate files and consensus sequences are transferred to a local GCP bucket for storage. Next, we run the ``SC2_lineage_calling_and_results`` which uses Pangolin and Nextclade to peform clade and lineage assignment on the consesnus assemblies and produce a results summary file for the set of sequences analyzed. 

If you already have a multifasta, you can use the ``SC2_multifasta_lineage_calling`` workflow for clade and lineage assignment. 

<br/>

### Wastewater SC2 sequence assembly and variant calling
Sequence assembly and variant calling for wastewater SC2 samples requires three workflows. Similiar to our process for clinical SC2 samples, we first run either the ``SC2_illumina_pe_assembly`` or ``SC2_ont_assembly`` workflows which perform quality control, trimming, and filtering of raw reads followed by reference-guided whole genome assembly and finally intermediate files and consensus sequences are transferred to a local GCP bucket for storage. Next, we run the ``SC2_lineage_calling_and_results`` which uses Pangolin and Nextclade to peform clade and lineage assignment on the consesnus assemblies and produce a results summary file for the set of sequences analyzed. Finally we run the ``SC2_wastewater_variant_calling workflow`` which uses Freyja to recover relative lineage abundances from wastewater sample which is considered a mixed SC2 sample. 

<br/>

![SC2 high level overview workflow diagram](./workflow_diagrams/SC2_overview_workflow_diagram.png "high level overview of SC2 workflow")

<br/>


## Setup

<br/>

### Input Workflow from Dockstore
To use the workflow on the Terra platform, first you will need to import the workflow from Dockstore. All workflows can be found under our dockstore organization called CDPHE-bioinformatics.
1. Go to dockstore (https://dockstore.org/).
2. Along the top search bar click on Organizations and search for "CDPHE".
3. Select the workflow.
4. On the right hand side of the workflow description, select "Launch with Terra"
5. Select the Destination workspace and select "Import". 
6. The workflow will now be displayed as a card under your workflows tab in your Terra workspace. 

<br/>

### Workspace Data
Prior to running any of the workflows, you must set up the Terra workspace data with the correct reference file and custom python scripts. The reference files can be found in this repository in the ``data/workspace_data`` directory. Python scripts can be found in the ``scripts`` directory. Workspace variables are named using the following format ``{organism}_{description}_{file_type}``, except for the primer bed files which are named as ``{desciption}_{file_type}``. Reference files and python scripts should be copied from this repo into a GCP bucket. The GCP bucket path to the file will serve as the "value" when adding data to the terra workspace data table. 

To add data to the terra workspace data:
1. Navigate to the Data tab in your Terra workspace.
2. In the left hand list of data tables, under "Other Data" select "Workspace Data".
3. Click on the "+" button in the lower right hand corner of the workspace data table. 
4. Fill in the "Key" column with the workspace variable name, the "Value" column with GCP bucket path to the file and the "Description" column with a brief description if desired. 
5. Once complete hit the check mark to the right.

Below is a data table detailing the workspace data you will need to set up in order to run the SC2 workflows. 


| workspace variable name | workflow|  file name | description | 
|-------------------|-----------------|--------------------|-----------------|
| ``adapters_and_contaminants_fa`` | ``SC2_illumina_pe_assembly``| Adapters_plus_PhiX_174.fasta | adapaters sequences and contaiment sequences removed during fastq cleaning and filtering using SeqyClean. Thanks to Erin Young at Utah Public Health Laboratory for providing this file!  |
| ``covid_genome_gff`` |``SC2_illumina_pe_assembly``, ``SC2_ont_assembly`` | NC_045512-2_reference.gff| whole genome reference sequence annotation file in gff format (we use NCBI genbank ID MN908947.3) |
| ``covid_genome_fa`` | ``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``|MN908947-2_reference.fasta | SARS-CoV-2 whole genome reference sequence in fasta format (we use NCBI genbank ID MN908947.3) |
| ``artic_v3_bed`` |  ``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``|artic_V3_nCoV-2019.primer.bed| primer bed file for the Artic V3 tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| ``artic_v4_bed`` | ``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``|artic_V4_nCoV-2019.primer.bed | primer bed file for the Artic V4 tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| ``artic_v4-1_bed`` | ``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``|artic_V4-1_nCoV-2019.primer.bed | primer bed file for the Artic V4.1 tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
|``artic_v4-1_s_gene_amplicons``|``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``|artic_v4_1_s_gene_amplicons.tsv|coordinate positions of S gene amplicons using the artic V4.1 primers|
|``artic_v4-1_s_gene_primer_bed``|``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``| S_gene_V4-1_nCoV-2021.primer.bed|primer sequences and coordinate positions of primer binding region of Artic v4.1 primers|
|``artic_v5-3-2_bed``| ``SC2_illumina_pe_assembly``, ``SC2_ont_assembly``| artic_v5-3-2_nCoV-2023.primer.bed|primer bed file for the Artic V5.3.2 tiled amplicon primer set. |
|``artic_v5-3-2_s_gene_amplicons``|``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``|artic_v5-3-2_s_gene_amplicons.tsv| coordinate positions of S gene amplicons using the artic V5.3.2 primers|
|``artic_v5-3-2_s_gene_primer_bed``|``SC2_illumina_pe_assembly``,  ``SC2_ont_assembly``| S_gene_V5-3-2_nCoV-2021.primer.bed| primer sequences and coordinate positions of primer binding region of Artic v5.3.2 primers|
| ``midnight_bed`` | ``SC2_ont_assembly``| Midnight_Primers_SARS-CoV-2.scheme.bed | primer bed file for the Midnight tiled amplicon primer set. Thanks to Theiagen Genomics for providing this file! |
| `` covid_voc_annotations_tsv``|``SC2_wastewater_variant_calling workflow``| SC2_voc_annotations_20220711.tsv | For wastewater only. List of amino acid (AA) subsitiutions and lineages containing those AA substitutions; for a lineage to be associated with a given AA subsitution, 90% of publically available sequences must contain the AA substition (the 90% cutoff was determined using outbreak.info) |
| ``covid_voc_bed_tsv``|``SC2_wastewater_variant_calling workflow``.| SC2_voc_mutations_20220711.tsv |  For wastewater only. List of nucleotide genome positions in relation to the MN908947.3 reference genome of know mutations |
| ``covid_calc_per_cov_py`` |``SC2_illumina_pe_assembly``, ``SC2_ont_assembly`` |calc_percent_coverage.py | see detailed description in the readme file found in ``./python_scripts/`` repo directory|
| ``covid_nextclade_json_parser_py`` | ``SC2_lineage_calling_and_results``| nextclade_json_parser.py | see detailed description in the readme file found in ``./python_scripts/`` repo directory|
| ``covid_concat_results_py`` | ``SC2_lineage_calling_and_results``| concat_seq_metrics_and_lineages_results.py | see detailed description in the readme file found in ``./python_scripts`` repo directory |
|``covid_version_capture_illumina_pe_assembly_py``|``SC2_illumina_pe_assembly``| version_capture_illumina_pe_assembly.py| generates version catpure output file for software versions used in the SC2_illumina_pe_assembly workflow|
|``covid_version_capture_ont_assembly_py``|``SC2_ont_assembly``| version_capture_ont_assembly.py| generates version catpure output file for software versions used in the SC2_ont_assembly workflow|
|``covid_version_capture_lineage_calling_py``|``SC2_lineage_calling_and_results``| version_capture_lineage_calling_and_results.py| generates version catpure output file for software versions used in the SC2_lineage_calling_and_results workflow|
|``covid_version_capture_wastewater_variant_calling_py``|``SC2_wastewater_variant_calling``| version_capture_wastewater_variant_calling.py| generates version catpure output file for software versions used in the SC2_wastewater_variant_calling workflow|
|``covid_version_capture_multifasta_lineage_calling_py``|``SC2_multifasta_lineage_calling``| version_capture_multifasta_lineage_calling.py| generates version catpure output file for software versions used in the SC2_multifasta_lineage_calling workflow|

