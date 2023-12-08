version 1.0

workflow SC2_novel_mutations {

    input {
        Array[String] project_names_array
        Array[File] combined_mutations_array
        String covwwt_path 
        String historical_data_path
        String today

        # Reference files
        File historical_full
        File historical_unique
        File metadata
        File gff_mutations

        # Python scripts
        File novel_mutations_append_py
    }


    call append_new_mutations {
        input:
            project_names_array = project_names_array,
            combined_mutations_array = combined_mutations_array,
            novel_mutations_append_py = novel_mutations_append_py,
            historical_full = historical_full,
            historical_unique = historical_unique,
            metadata = metadata,
            gff_mutations = gff_mutations,
            today = today
    }

    scatter (project in append_new_mutations.project_unique_mutations) {
        call transfer_project_outputs {
            input:
                project_unique_mutations = project,
                covwwt_path = covwwt_path
        }    
    }

    call transfer_appended_outputs {
        input:
            historical_full_updated = append_new_mutations.historical_full_updated,
            historical_unique_updated = append_new_mutations.historical_unique_updated,
            historical_data_path = historical_data_path,
            recurrent_mutations = append_new_mutations.recurrent_mutations,
            novel_mutations = append_new_mutations.novel_mutations,
            covwwt_path = covwwt_path
    }

    output {
        Array[File] project_unique_mutations = append_new_mutations.project_unique_mutations
        File historical_full_updated = append_new_mutations.historical_full_updated
        File historical_unique_updated = append_new_mutations.historical_unique_updated
        Array[File]? project_missing_dates = append_new_mutations.project_missing_dates
        File? recurrent_mutations = append_new_mutations.recurrent_mutations
        File? novel_mutations = append_new_mutations.novel_mutations
    }
}

task append_new_mutations {

    meta {
        volatile: true
    }

    input {
        Array[String] project_names_array
        Array[File] combined_mutations_array
        String today
        File historical_full
        File historical_unique
        File metadata
        File gff_mutations
        File novel_mutations_append_py
    }

    command <<<
        date +"%d-%m-%Y" > today

        python ~{novel_mutations_append_py} --project_names ~{sep=' ' project_names_array}\
        --combined_mutations_files ~{sep=' ' combined_mutations_array} \
        --historical_full ~{historical_full} \
        --historical_unique ~{historical_unique} \
        --metadata ~{metadata} \
        --gff ~{gff_mutations} \
        --today ~{today}
    >>>

    output {
        Array[File]? project_missing_dates = glob("*_missing_dates.tsv")
        Array[File] project_unique_mutations = glob("*_unique_mutations.tsv")
        File historical_full_updated = "novel_mutations_historical_full.tsv"
        File historical_unique_updated = "novel_mutations_historical_unique.tsv"        
        File? recurrent_mutations = "recurrent_mutations_{today}.tsv"
        File? novel_mutations = "novel_mutations_{today}.tsv"
    }

    runtime {
        docker: "ariannaesmith/py3.10.9-bio:v1"
        memory: "1 GB"
        cpu: 4
        disks: "local-disk 10 SSD"
    }
}

task transfer_project_outputs {
    input {
        File project_unique_mutations
        String covwwt_path
    }

    String project_name = basename(project_unique_mutations, "_unique_mutations.tsv")

    command <<<
        gsutil -m cp ~{project_unique_mutations} ~{covwwt_path}/~{project_name}/novel_mutations/
    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}

task transfer_appended_outputs {

    input {
        File historical_full_updated
        File historical_unique_updated
        String historical_data_path
        File? recurrent_mutations
        File? novel_mutations
        String covwwt_path
    }

    command <<<
        gsutil -m cp ~{historical_full_updated} ~{historical_data_path}
        gsutil -m cp ~{historical_unique_updated} ~{historical_data_path}
        gsutil -m cp ~{recurrent_mutations} ~{covwwt_path}/recurrent_mutations
        gsutil -m cp ~{novel_mutations} ~{covwwt_path}/novel_mutations
    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}