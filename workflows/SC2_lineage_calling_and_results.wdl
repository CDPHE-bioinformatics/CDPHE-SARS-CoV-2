version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture

workflow SC2_lineage_calling_and_results {

    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File?] cov_out
        Array[File?] percent_cvg_csv
        Array[String] out_dir_array
        Array[String] project_name_array
        Array[File] terra_data_table_path_array

        # workspace data
        File cdc_lineage_groups_json
        
        # python scripts
        File nextclade_json_parser_py
        File concat_seq_results_py
        File version_capture_lineage_calling_and_results_py

    }

    # secret variables - for static values convert from array to single entity
    String project_name = select_all(project_name_array)[0]
    File terra_data_table_path = select_all(terra_data_table_path_array)[0]
    String out_dir = select_all(out_dir_array)[0]
    String outdirpath = sub(out_dir, "/$", "")

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

    call version_capture.workflow_version_capture as workflow_version_capture{
        input:
    }

    call results_table {
      input:
        sample_name = sample_name,
        concat_seq_results_py = concat_seq_results_py,
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        pangolin_lineage_csv = pangolin.lineage,
        cdc_lineage_groups_json = cdc_lineage_groups_json,
        nextclade_clades_csv = parse_nextclade.nextclade_clades_csv,
        nextclade_variants_csv = parse_nextclade.nextclade_variants_csv,
        nextclade_version = nextclade.nextclade_version,
        project_name = project_name,
        terra_data_table_path = terra_data_table_path,
        workflow_version = workflow_version_capture.workflow_version,
        analysis_date = workflow_version_capture.analysis_date
    }




    call create_version_capture_file {
        input: 
            version_capture_lineage_calling_and_results_py = version_capture_lineage_calling_and_results_py,
            pangolin_version = pangolin.pangolin_version,
            nextclade_version = nextclade.nextclade_version,
            project_name = project_name,
            workflow_version = workflow_version_capture.workflow_version,
            analysis_date = workflow_version_capture.analysis_date,
            pangolin_lineage_csv = pangolin.lineage
    }

    call transfer {
      input:
          outdirpath = outdirpath,
          cat_fastas = concatenate.cat_fastas,
          pangolin_lineage = pangolin.lineage,
          nextclade_json = nextclade.nextclade_json,
          nextclade_csv = nextclade.nextclade_csv,
          nextclade_clades_csv = parse_nextclade.nextclade_clades_csv,
          nextclade_variants_csv = parse_nextclade.nextclade_variants_csv,
          sequencing_results_csv = results_table.sequencing_results_csv,
          version_capture_lineage_calling_and_results = create_version_capture_file.version_capture_lineage_calling_and_results
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
        File version_capture_lineage_calling_and_results = create_version_capture_file.version_capture_lineage_calling_and_results
        String transfer_date_lineage_calling = transfer.transfer_date_lineage_calling
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

    command <<<

        pangolin --version > VERSION

        pangolin --skip-scorpio --expanded-lineage --threads 32 \
            --outfile pangolin_lineage_report.csv ~{cat_fastas}

    >>>

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

    command <<<
      python ~{nextclade_json_parser_py} \
          --nextclade_json ~{nextclade_json} \
          --project_name ~{project_name}

    >>>

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
      File concat_seq_results_py
      Array[File] cov_out
      Array[File] percent_cvg_csv
      File pangolin_lineage_csv
      File cdc_lineage_groups_json
      File nextclade_clades_csv
      File nextclade_variants_csv
      String nextclade_version
      String project_name
      File terra_data_table_path
      String workflow_version
      String analysis_date
      

    }

    command <<<
    python ~{concat_seq_results_py} \
        --sample_name_array "~{write_lines(sample_name)}" \
        --cov_out_files "~{write_lines(cov_out)}" \
        --percent_cvg_files "~{write_lines(percent_cvg_csv)}" \
        --pangolin_lineage_csv "~{pangolin_lineage_csv}" \
        --cdc_lineage_groups_json "~{cdc_lineage_groups_json}" \
        --nextclade_variants_csv "~{nextclade_variants_csv}" \
        --nextclade_clades_csv "~{nextclade_clades_csv}" \
        --nextclade_version "~{nextclade_version}" \
        --project_name "~{project_name}" \
        --terra_data_table_path "~{terra_data_table_path}" \
        --workflow_version "~{workflow_version}" \
        --analysis_date "~{analysis_date}"

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results_v~{workflow_version}.csv"

    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task create_version_capture_file {
    meta {
        description: "generate version capture file for workflow"
    }

    input {
        File version_capture_lineage_calling_and_results_py
        String project_name
        String pangolin_version
        String nextclade_version
        String analysis_date
        String workflow_version
        File pangolin_lineage_csv

    }

    command <<<

        python ~{version_capture_lineage_calling_and_results_py} \
        --project_name "~{project_name}" \
        --pangolin_version "~{pangolin_version}" \
        --nextclade_version "~{nextclade_version}" \
        --analysis_date "~{analysis_date}" \
        --workflow_version "~{workflow_version}" \
        --pangolin_lineage_csv "~{pangolin_lineage_csv}"

    >>>

    output {
        File version_capture_lineage_calling_and_results = "version_capture_lineage_calling_and_results_~{project_name}_~{workflow_version}.csv"
    }

    runtime {

      docker: "mchether/py3-bio:v4"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }
}

task transfer {
    input {
        String outdirpath
        File cat_fastas
        File pangolin_lineage
        File nextclade_json
        File nextclade_csv
        File nextclade_clades_csv
        File nextclade_variants_csv
        File sequencing_results_csv
        File version_capture_lineage_calling_and_results
    }


    command <<<

        gsutil -m cp ~{cat_fastas} ~{outdirpath}/multifasta/
        gsutil -m cp ~{pangolin_lineage} ~{outdirpath}/pangolin_out/
        gsutil -m cp ~{nextclade_json} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_csv} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_clades_csv} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_variants_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{version_capture_lineage_calling_and_results} ~{outdirpath}/summary_results/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    
    >>>

    output {

        String transfer_date_lineage_calling = read_string("TRANSFER_DATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
