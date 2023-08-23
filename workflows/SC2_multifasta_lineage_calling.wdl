version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture

workflow SC2_multifasta_lineage_calling {

    input {
        File multifasta
        String sample_id
        String out_dir

        # workspace references
        File version_capture_multifasta_lineage_calling_py

        # secret variables
        String outdirpath = sub(out_dir, "/$", "")
    }

    call nextclade {
        input:
            multifasta = multifasta,
            sample_id = sample_id
    }
    
    call pangolin {
        input:
            multifasta = multifasta,
            sample_id = sample_id
    }

    call version_capture.workflow_version_capture as workflow_version_capture {
        input:
    }
    
    call_create_version_catpure_file {
        input:
            version_capture_multifasta_lineage_calling_py = version_capture_multifasta_lineage_calling_py,
            pangolin_version = pangolin.pangolin_version,
            nextclade_version = nextclade.nextclade_version, 
            analysis_date = workflow_version_capture.analysis_date,
            workflow_version = workflow_version_capture.workflow_version     

    } 

    call transfer {
      input:
          outdirpath = outdirpath,
          nextclade_json = nextclade.nextclade_json,
          auspice_json = nextclade.auspice_json,
          nextclade_csv = nextclade.nextclade_csv,
          pangolin_lineage = pangolin.lineage,
          version_capture_multifasta_lineage_calling = create_version_capture_file.version_capture_multifasta_lineage_calling
    }

    output {
        String nextclade_version = nextclade.nextclade_version
        File nextclade_json = nextclade.nextclade_json
        File auspice_json = nextclade.auspice_json
        File nextclade_csv = nextclade.nextclade_csv
        String pangolin_version = pangolin.pangolin_version
        File pangolin_lineage = pangolin.lineage
        File version_capture_multifasta_lineage_calling = create_version_capture_file.version_capture_multifasta_lineage_calling
        String transfer_date = transfer.transfer_date
    }
}

task nextclade {

    input {
        File multifasta
        String sample_id
    }

    command {

        nextclade --version > VERSION
        nextclade dataset get --name='sars-cov-2' --reference='MN908947' --output-dir='data/sars-cov-2'
        nextclade run --input-dataset data/sars-cov-2 --output-json ${sample_id}_nextclade.json --output-csv ${sample_id}_nextclade.csv --output-tree ${sample_id}_nextclade.auspice.json ${multifasta}
    }

    output {
        String nextclade_version = read_string("VERSION")
        File nextclade_json = "${sample_id}_nextclade.json"
        File auspice_json = "${sample_id}_nextclade.auspice.json"
        File nextclade_csv = "${sample_id}_nextclade.csv"
    }

    runtime {
        docker: "nextstrain/nextclade"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 500 HDD"
    }
}

task pangolin {

    input {
        File multifasta
        String sample_id
    }

    command {
        pangolin --version > VERSION
        pangolin --skip-scorpio --outfile ${sample_id}_pangolin_lineage_report.csv ${multifasta}
    }

    output {
        String pangolin_version = read_string("VERSION")
        File lineage = "${sample_id}_pangolin_lineage_report.csv"
    }

    runtime {
        docker: "staphb/pangolin"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 500 HDD"
    }
}

task create_version_capture_file {

    input {

        File version_capture_multifasta_lineage_calling_py
        String pangolin_version
        String nextclade_version
        String analysis_date
        String workflow_version
    }

    command <<<

        python ~{version_capture_multifasta_lineage_calling_py} \
        --pangolin_version "~{pangolin_version}" \
        --nextclade_version "~{nextclade_version" \
        --analysis_date "~{analysis_date}" \
        --workfow_version "~{workflow_version}"
     >>>

    output {
        File version_capture_multifasta_lineage_calling = "version_capture_multifasta_lineage_calling_~{project_name}_~{workflow_version}.csv"
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
        File auspice_json
        File nextclade_csv
        File nextclade_json
        File pangolin_lineage
        File version_capture_multifasta_lineage_calling
    }

    String outdir = sub(out_dir, "/$", "")

    command <<<

        gsutil -m cp ~{nextclade_json} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{auspice_json} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{nextclade_csv} ~{outdirpath}/nextclade_out/
        gsutil -m cp ~{pangolin_lineage} ~{outdirpath}/pangolin_out/
        gsutil -m cp ~{version_capture_multifasta_lineage_calling} ~{outdirpath}/summary_results/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE

    >>>

    output {
        transfer_date = read_string("TRANSFER_DATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
