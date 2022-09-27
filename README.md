# CDPHE SAR-CoV-2 repository

In this repository you will find all things related to SARS-CoV-2 bioinformatics and analysis. The general directory stucture is as follows:

```
CDPHE_SARS-CoV-2
| README_SC2_overview.md # this document
|____ workflows # terra workflows and asscoiated scripts and workspace data
      | README_CDPHE_SARS-Cov-2_workflow.md
      | dockstore.yml
      |____ wdls # workflows you can find in dockstore and push to your terra workspace
            | SC2_illumina_pe_assembly.wdl
            | SC2_illumina_se_assembly.wdl
            | SC2_ont_assembly.wdl
            | SC2_transfer_illumina_pe_assembly.wdl
            | SC2_transfer_illumina_se_assembly.wdl
            | SC2_transfer_ont_assembly.wdl
            | SC2_lineage_calling_and_results.wdl
            | SC2_multifasta_lineage_calling.wdl
            | SC2_wastewater_variant_calling.wdl
      |____ python_scripts # custom python scripts for analysis used in wdl workflows
            | README_SC2_preprocess_python_script.md
            | calc_percent_coverage.py
            | concat_seq_metrics_and_lineage.py
            | nextclade_json_parser.py
      |____ workspace_data # reference data
            | Adapters_plu_PhiX_174.fasta
            | reference.fasta
            | referece.gff
            | v3_nCoV-2019.primer.bed
            | v4_nCoV-2021.primer.bed
            | v4.1_nCoV-2021.primer.bed
            | Midnight_Primers_SARS-CoV-2.scheme.bed
            | spike_mutations.tsv
            | spike_VOC-VUI_annotations.tsv
      |____ workflow_diagrams
|____ preprocess_python_scripts # python scripts for getting data ready for wdl workflows
      | README_SC2_preprocess_python_script.md
      | SC2_illumina_preprocess.py
      | SC2_ont_preprocess.py
      | environment.yml
      |____ example_data
      |____ workflow_diagrams
|____ SC2_sequence_submission_clincial_samples
      | README_SC2_sequence_submission_clinical_samples.md
|____ SC2_sequence_submission_wwt_samples
      | README_SC2_sequence_submission_wwt_samples.md
|____ SC2_indel_finder
      | README_SC2_indel_finder_workflows.md
      | indel_finder.py
      | environment.yml
```

