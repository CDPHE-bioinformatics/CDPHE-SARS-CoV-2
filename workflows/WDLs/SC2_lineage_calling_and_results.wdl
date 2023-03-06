version 1.0

workflow SC2_lineage_calling_and_results {

    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File?] cov_out
        Array[File?] percent_cvg_csv
        Array[String] out_dir_array
        # Array[String] plate_name
        # Array[String] plate_sample_well
        # Array[String] primer_set
        Array[String] project_name_array
        # Array[String] tech_platform
        # Array[String] read_type
        Array[String?] assembler_version_array
        Array[File] workbook_path_array

        # python scripts
        File nextclade_json_parser_py
        File concat_seq_results_py



    }

    # secret variables - for static values convert from array to single entity
    String project_name = project_name_array[0]
    File workbook_path = workbook_path_array[0]
    String assembler_version = select_all(assembler_version_array)[0]
    String out_dir = out_dir_array[0]

    call concatenate {
        input:
            renamed_consensus = select_all(renamed_consensus)
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
        project_name = project_name,
        nextclade_json_parser_py = nextclade_json_parser_py,
        nextclade_json = nextclade.nextclade_json
    }

    call results_table {
      input:
        sample_name = sample_name,
        # plate_name =  plate_name,
        # plate_sample_well = plate_sample_well,
        # primer_set = primer_set,
        # tech_platform = tech_platform,
        # read_type = read_type,
        concat_seq_results_py = concat_seq_results_py,
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        pangolin_lineage_csv = pangolin.lineage,
        pangolin_version = pangolin.pangolin_version,
        nextclade_clades_csv = parse_nextclade.nextclade_clades_csv,
        nextclade_variants_csv = parse_nextclade.nextclade_variants_csv,
        nextclade_version = nextclade.nextclade_version,
        project_name = project_name,
        assembler_version= assembler_version

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
        File wgs_horizon_report_csv = results_table.wgs_horizon_report_csv
    }
}

task concatenate {

    input {
        Array[File] renamed_consensus
    }

    command <<<

        cat ~{sep=" " renamed_consensus} > concatenate_assemblies.fasta

    >>>

    output {

        File cat_fastas = "concatenate_assemblies.fasta"

    }

    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 10 SSD"
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
      File nextclade_json_parser_py
      File nextclade_json
      String project_name
    }

    command {
      python ~{nextclade_json_parser_py} \
          --nextclade_json ~{nextclade_json} \
          --project_name ~{project_name}

    }

    output {
      File nextclade_clades_csv = '${project_name}_nextclade_results.csv'
      File nextclade_variants_csv = '${project_name}_nextclade_variant_summary.csv'
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
      Array[String] sample_name
    #   Array[String] plate_name
    #   Array[String] plate_sample_well
    #   Array[String] primer_set
    #   Array[String] tech_platform
    #   Array[String] read_type
      File concat_seq_results_py
      Array[File] cov_out
      Array[File] percent_cvg_csv
      File pangolin_lineage_csv
      String pangolin_version
      File nextclade_clades_csv
      File nextclade_variants_csv
      String nextclade_version
      String project_name
      String assembler_version
      File workbook_path

    }

    command <<<
        python ~{concat_seq_results_py} \
            --sample_name_array ~{write_lines(sample_name)} \ 
            --workbook_path ~{workbook_path} \ 
            --cov_out_files ~{cov_out} \
            --percent_cvg_files ~{write_lines(percent_cvg_csv)} \
            --assembler_version ~{assembler_version} \
            --pangolin_lineage_csv ~{pangolin_lineage_csv} \
            --pangolin_version ~{pangolin_version} \
            --nextclade_clades_csv ~{nextclade_clades_csv} \
            --nextclade_version ~{nextclade_version} \ 
            --project_name ~{project_name} 

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv"
        File wgs_horizon_report_csv = "~{project_name}_wgs_horizon_report.csv"
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
        String out_dir
        File cat_fastas
        File pangolin_lineage
        File nextclade_json
        File nextclade_csv
        File nextclade_clades_csv
        File nextclade_variants_csv
        File sequencing_results_csv
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
        gsutil -m cp ~{wgs_horizon_report_csv} ~{outdirpath}/summary_results/
    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
