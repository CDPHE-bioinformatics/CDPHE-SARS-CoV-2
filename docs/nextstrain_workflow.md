# Nextstrain Workflow

We use the publicly available [Nextstrain workflow](https://dockstore.org/workflows/github.com/broadinstitute/viral-pipelines/sarscov2_nextstrain:master?tab=info) to generate Nextstrain builds, and then transfer the results using this transfer workflow.

This workflow transfers the output file generated from the publicly available [sarscov2_nextstrain workflow](https://dockstore.org/workflows/github.com/broadinstitute/viral-pipelines/sarscov2_nextstrain:master?tab=info) to a user specified google bucket. Below is a summary of the workflow input variables along with the syntax used for the attribute column when setting up the workflow to run on Terra.bio. For the attributes, the "this.sample{terra_data_table_name}s." syntax tells Terra to pull the variable from the sample-level terra data table. The Google Bucket path describes where in the user google bucket the output file is transferred to.  

## Inputs

| workflow variable      | attribute (input syntax into workflow)       | google bucket path                                     |
| ---------------------- | -------------------------------------------- | ------------------------------------------------------ |
| `auspice_input_json`   | this.auspice_input_json                      | `gs://{user_defined_gcp_bucket}/auspice_input_json/`   |
| `combined_assemblies`  | this.combined_assemblies                     | `gs://{user_defined_gcp_bucket}/combined_assemblies/`  |
| `keep_list`            | this.keep_list                               | `gs://{user_defined_gcp_bucket}/keep_list/`            |
| `metadata_merged`      | this.metadata                                | `gs://{user_defined_gcp_bucket}/metadata_merged/`      |
| `ml_tree`              | this.ml_tree                                 | `gs://{user_defined_gcp_bucket}/ml_tree/`              |
| `multiple_alignment`   | this.multiple_alignment                      | `gs://{user_defined_gcp_bucket}/multiple_alignment/`   |
| `node_data_jsons`      | this.node_data_jsons                         | `gs://{user_defined_gcp_bucket}/node_data_jsons/`      |
| `root_sequence_json`   | this.root_sequence_json                      | `gs://{user_defined_gcp_bucket}/root_sequence_json/`   |
| `subsampled_sequences` | this.subsampled_sequences                    | `gs://{user_defined_gcp_bucket}/subsampled_sequences/` |
| `time_tree`            | this.time_tree                               | `gs://{user_defined_gcp_bucket}/time_tree/`            |
| `tip_frequencies_json` | this.tip_frequencies_json                    | `gs://{user_defined_gcp_bucket}/tip_frequencies_json/` |
| `unmasked_snps`        | this.unmasked_snps                           | `gs://{user_defined_gcp_bucket}/unmasked_snps/`        |
| `out_dir`              | `gs://{user_defined_gcp_bucket}/assemblies/` | NA                                                     |
