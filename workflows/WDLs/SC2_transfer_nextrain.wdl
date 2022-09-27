version 1.0

workflow SC2_transfer_nextrain {

    input {
        Array[File] auspice_input_json
        Array[File] combined_assemblies
        Array[File] keep_list
        Array[File] metadata_merged
        Array[File] ml_tree
        Array[File] multiple_alignment
        Array[File] node_data_jsons
        Array[File] root_sequence_json
        Array[File] subsampled_sequences
        Array[File] time_tree
        Array[File] tip_frequencies_json
        Array[File] unmasked_snps
        String out_dir
    }

    call transfer_outputs {
        input:
            auspice_input_json = auspice_input_json,
            combined_assemblies = combined_assemblies,
            keep_list = keep_list,
            metadata_merged = metadata_merged,
            ml_tree = ml_tree,
            multiple_alignment = multiple_alignment,
            node_data_jsons = node_data_jsons,
            root_sequence_json = root_sequence_json,
            subsampled_sequences = subsampled_sequences,
            time_tree = time_tree,
            tip_frequencies_json = tip_frequencies_json,
            unmasked_snps = unmasked_snps,
            out_dir = out_dir
    }

    output {
        String transfer_date = transfer_outputs.transfer_date
    }
}

task transfer_outputs {
    input {
        String out_dir
        Array[File] auspice_input_json
        Array[File] combined_assemblies
        Array[File] keep_list
        Array[File] metadata_merged
        Array[File] ml_tree
        Array[File] multiple_alignment
        Array[File] node_data_jsons
        Array[File] root_sequence_json
        Array[File] subsampled_sequences
        Array[File] time_tree
        Array[File] tip_frequencies_json
        Array[File] unmasked_snps
        String out_dir
    }

    String outdir = sub(out_dir, "/$", "")

    command <<<

        gsutil -m cp ~{sep=' ' auspice_input_json} ~{outdir}/auspice_input_json/
        gsutil -m cp ~{sep=' ' combined_assemblies} ~{outdir}/combined_assemblies/
        gsutil -m cp ~{sep=' ' keep_list} ~{outdir}/keep_list/
        gsutil -m cp ~{sep=' ' metadata_merged} ~{outdir}/metadata_merged/
        gsutil -m cp ~{sep=' ' ml_tree} ~{outdir}/ml_tree/
        gsutil -m cp ~{sep=' ' multiple_alignment} ~{outdir}/multiple_alignment/
        gsutil -m cp ~{sep=' ' node_data_jsons} ~{outdir}/node_data_jsons/
        gsutil -m cp ~{sep=' ' root_sequence_json} ~{outdir}/root_sequence_json/
        gsutil -m cp ~{sep=' ' subsampled_sequences} ~{outdir}/subsampled_sequences/
        gsutil -m cp ~{sep=' ' time_tree} ~{outdir}/time_tree/
        gsutil -m cp ~{sep=' ' tip_frequencies_json} ~{outdir}/tip_frequencies_json/
        gsutil -m cp ~{sep=' ' unmasked_snps} ~{outdir}/unmasked_snps/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}
