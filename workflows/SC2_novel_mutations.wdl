version 1.1

workflow SC2_novel_mutations {

    input {
        Array[String] project_names_array
        Array[File] combined_mutations_array
        Array[String] out_dir_array 
        String historical_full_path
        String historical_unique_path


        # Reference files
        File historical_full
        File historical_unique

        # Python scripts
        File novel_mutations_append_py
    }
    # Secret variables
    String out_dir = out_dir_array[0]
    String out_dir_path = sub(out_dir, "/$", "")

    call append_new_mutations {
        input:
            project_names_array = project_names_array,
            combined_mutations_array = combined_mutations_array,
            novel_mutations_append_py = novel_mutations_append_py,
            historical_full = historical_full,
            historical_unique = historical_unique
    }

    scatter (project in append_new_mutations.project_unique_mutations) {
        call transfer_project_outputs {
            input:
                project_unique_mutations = project,
                out_dir_path = out_dir_path
        }    
    }

    call transfer_appended_outputs {
        input:
            historical_full_updated = append_new_mutations.historical_full_updated,
            historical_unique_updated = append_new_mutations.historical_unique_updated,
            historical_full_path = historical_full_path,
            historical_unique_path = historical_unique_path
    }

    output {
        Array[File] project_unique_mutations = append_new_mutations.project_unique_mutations
        File historical_full_updated = append_new_mutations.historical_full_updated
        File historical_unique_updated = append_new_mutations.historical_unique_updated
    }
}

task append_new_mutations {

    input {
        Array[String] project_names_array
        Array[File] combined_mutations_array
        File historical_full
        File historical_unique
        File novel_mutations_append_py
    }

    command <<<
        python ~{novel_mutations_append_py} --project_names ~{sep(',', project_names_array)}\
        --combined_mutations_files ~{sep(',', combined_mutations_array)} \
        --historical_full ~{historical_full} \
        --historical_unique ~{historical_unique}
    >>>

    output {
        File historical_full_updated = "historical_full_updated.csv"
        File historical_unique_updated = "historical_unique_updated.csv"
        Array[File] project_unique_mutations = glob("*_unique_mutations.csv")
    }

    runtime {
        docker: "mchether/py3-bio:v4"
        memory: "1 GB"
        cpu: 4
        disks: "local-disk 10 SSD"
    }
}

task transfer_project_outputs {
    input {
        File project_unique_mutations
        String out_dir_path
    }

    String project_name = basename(project_unique_mutations, "_unique_mutations.csv")

    command <<<
        gsutil -m cp ~{project_unique_mutations} ~{out_dir_path}/~{project_name}/novel_mutations/
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
        String historical_full_path
        String historical_unique_path
    }

    command <<<
        gsutil -m cp ~{historical_full_updated} ~{historical_full_path}
        gsutil -m cp ~{historical_unique_updated} ~{historical_unique_path}
    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}