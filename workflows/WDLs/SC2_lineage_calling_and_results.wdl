version 1.0

workflow SC2_lineage_calling_and_results {

    input {
        Array[String] sample_id
        Array[File?] assembly_fastas
        Array[File?] cov_out_txt
        Array[File?] percent_cvg_csv
        File nextclade_json_parser_script
        File concat_results_script
        Array[String] out_dir
        Array[String] plate_name
        Array[String] plate_sample_well
        Array[String] primer_set
        Array[String] seq_run
        Array[String] tech_platform
        Array[String] read_type
        Array[String?] assembler_version

    }

    call concatenate {
        input:
            assembly_fastas_non_empty = select_all(assembly_fastas)
    }

    call pangolin {
        input:
            cat_fastas = concatenate.cat_fastas
    }

    call nextclade {
        input:
            multifasta = concatenate.cat_fastas
    }

    call parse_nextclade {
      input:
        seq_run = seq_run,
        nextclade_json_parser_script = nextclade_json_parser_script,
        nextclade_json = nextclade.nextclade_json
    }

    call results_table {
      input:
        sample_id = sample_id,
        plate_name =  plate_name,
        plate_sample_well = plate_sample_well,
        primer_set = primer_set,
        tech_platform = tech_platform,
        read_type = read_type,
        concat_results_script = concat_results_script,
        cov_out_txt_non_empty = select_all(cov_out_txt),
        percent_cvg_csv_non_empty = select_all(percent_cvg_csv),
        pangolin_lineage_csv = pangolin.lineage,
        pangolin_version = pangolin.pangolin_version,
        nextclade_clades_csv = parse_nextclade.nextclade_clades_csv,
        nextclade_variants_csv = parse_nextclade.nextclade_variants_csv,
        nextclade_version = nextclade.nextclade_version,
        seq_run = seq_run,
        assembler_version_non_empty = select_all(assembler_version)

    }

    call transfer {
      input:
          out_dir = out_dir,
          cat_fastas = concatenate.cat_fastas,
          pangolin_lineage = pangolin.lineage,
          nextclade_json = nextclade.nextclade_json,
          nextclade_csv = nextclade.nextclade_csv,
          nextclade_clades_csv = parse_nextclade.nextclade_clades_csv,
          nextclade_variants_csv = parse_nextclade.nextclade_variants_csv,
          sequencing_results_csv = results_table.sequencing_results_csv,
          sequence_assembly_metrics_csv = results_table.sequence_assembly_metrics_csv,
          wgs_horizon_report_csv = results_table.wgs_horizon_report_csv
    }

    output {
        File cat_fastas = concatenate.cat_fastas
        String pangolin_version = pangolin.pangolin_version
        File pangolin_lineage = pangolin.lineage
        String nextclade_version = nextclade.nextclade_version
        File nextclade_json = nextclade.nextclade_json
        File nextclade_csv = nextclade.nextclade_csv
        File nextclade_clades_csv = parse_nextclade.nextclade_clades_csv
        File nextclade_variants_csv = parse_nextclade.nextclade_variants_csv
        File sequencing_results_csv = results_table.sequencing_results_csv
        File sequence_assembly_metrics_csv = results_table.sequence_assembly_metrics_csv
        File wgs_horizon_report_csv = results_table.wgs_horizon_report_csv
    }
}

task concatenate {

    input {
        Array[File] assembly_fastas_non_empty
    }

    command <<<

        cat ~{sep=" " assembly_fastas_non_empty} > concatenate_assemblies.fasta

    >>>

    output {

        File cat_fastas = "concatenate_assemblies.fasta"

    }

    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task pangolin {

    input {

        File cat_fastas
    }

    command {

        pangolin --version > VERSION

        pangolin --skip-scorpio --expanded-lineage --threads 32 \
            --outfile pangolin_lineage_report.csv ${cat_fastas}

    }

    output {

        String pangolin_version = read_string("VERSION")
        File lineage = "pangolin_lineage_report.csv"

    }

    runtime {
        cpu:    32
        memory:    "16 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/pangolin"
    }
}

task nextclade {

    input {
        File multifasta
    }

    command {
        nextclade --version > VERSION
        nextclade dataset get --name='sars-cov-2' --reference='MN908947' --output-dir='data/sars-cov-2'
        nextclade run --input-dataset data/sars-cov-2 --output-json nextclade.json --output-csv nextclade.csv ${multifasta}
    }

    output {
        String nextclade_version = read_string("VERSION")
        File nextclade_json = "nextclade.json"
        File nextclade_csv = "nextclade.csv"
    }

    runtime {
        docker: "nextstrain/nextclade"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task parse_nextclade {

    input {
      File nextclade_json_parser_script
      File nextclade_json
      Array[String] seq_run
    }

    command {
      python ~{nextclade_json_parser_script} \
          --nextclade_json ~{nextclade_json} \
          --seq_run_file_list ${write_lines(seq_run)}

    }

    output {
      File nextclade_clades_csv = '${seq_run[0]}_nextclade_results.csv'
      File nextclade_variants_csv = '${seq_run[0]}_nextclade_variant_summary.csv'
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task results_table {

    input {
      Array[String] sample_id
      Array[String] plate_name
      Array[String] plate_sample_well
      Array[String] primer_set
      Array[String] tech_platform
      Array[String] read_type
      File concat_results_script
      Array[File] cov_out_txt_non_empty
      Array[File] percent_cvg_csv_non_empty
      File pangolin_lineage_csv
      String pangolin_version
      File nextclade_clades_csv
      File nextclade_variants_csv
      String nextclade_version
      Array[String] seq_run
      Array[String] assembler_version_non_empty
    }

    command {
      python ~{concat_results_script} \
          --sample_list ${write_lines(sample_id)} \
          --plate_name_file_list ${write_lines(plate_name)} \
          --plate_sample_well_file_list ${write_lines(plate_sample_well)}\
          --primer_set_file_list ${write_lines(primer_set)}\
          --tech_platform_file_list ${write_lines(tech_platform)}\
          --read_type ${write_lines(read_type)}\
          --bam_file_list ${write_lines(cov_out_txt_non_empty)} \
          --percent_cvg_file_list ${write_lines(percent_cvg_csv_non_empty)} \
          --pangolin_lineage_csv ~{pangolin_lineage_csv} \
          --pangolin_version "~{pangolin_version}" \
          --nextclade_clades_csv ~{nextclade_clades_csv} \
          --nextclade_variants_csv ~{nextclade_variants_csv} \
          --nextclade_version "~{nextclade_version}" \
          --seq_run ${write_lines(seq_run)}\
          --assembler_version_table_list ${write_lines(assembler_version_non_empty)}
    }

    output {
        File sequencing_results_csv = "${seq_run[0]}_sequencing_results.csv"
        File sequence_assembly_metrics_csv = "${seq_run[0]}_sequence_assembly_metrics.csv"
        File wgs_horizon_report_csv = "${seq_run[0]}_wgs_horizon_report.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task transfer {
    input {
        Array[String] out_dir
        File cat_fastas
        File pangolin_lineage
        File nextclade_json
        File nextclade_csv
        File nextclade_clades_csv
        File nextclade_variants_csv
        File sequencing_results_csv
        File sequence_assembly_metrics_csv
        File wgs_horizon_report_csv
    }

    String outdir = '${out_dir[0]}'
    String outdirpath = sub(outdir, "/$", "")

    command <<<

        gsutil -m cp ~{cat_fastas} ~{outdirpath}/multifasta/
        gsutil -m cp ~{pangolin_lineage} ~{outdirpath}/pangolin_out/
        gsutil -m cp ~{nextclade_json} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_csv} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_clades_csv} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_variants_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{sequence_assembly_metrics_csv} ~{outdirpath}/
        gsutil -m cp ~{wgs_horizon_report_csv} ~{outdirpath}/summary_results/
    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
