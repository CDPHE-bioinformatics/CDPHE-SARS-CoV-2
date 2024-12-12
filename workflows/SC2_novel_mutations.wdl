version 1.0

import "../tasks/transfer_task.wdl" as transfer_task

workflow SC2_novel_mutations {

    input {
        Array[String] project_names_array
        Array[File] combined_mutations_array
        String covwwt_path 
        Boolean overwite_project_and_set
        String historical_data_path
        Boolean overwrite_historical
        String today
        Array[String] sites_to_drop

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
            today = today,
            sites_to_drop = sites_to_drop
    }


    Boolean recurrent_mutations_defined = defined(append_new_mutations.recurrent_mutations)
    Boolean novel_mutations_defined = defined(append_new_mutations.novel_mutations)
    

    scatter (project_file in append_new_mutations.project_unique_mutations) {

        String project_name = basename(project_file, "_unique_mutations.tsv")
        String project_out_dir = "~{covwwt_path}/~{project_name}/terra_outputs"

        SubdirsToFiles project_subdirs_to_files = object { subdirs_to_files: [
            (project_file, "novel_mutations")
        ]}

        call transfer_task.transfer as transfer_project_outputs {
            input:
                out_dir = project_out_dir,
                overwrite = overwite_project_and_set,
                subdirs_to_files = project_subdirs_to_files
        }    
    }

    String covwwt_novel_mutations_path = "~{covwwt_path}/novel_mutations/"

    if (recurrent_mutations_defined || novel_mutations_defined) {

        SubdirsToFiles set_subdirs_to_files = object { subdirs_to_files: [
            (append_new_mutations.recurrent_mutations, ""),
            (append_new_mutations.novel_mutations, "")
        ]}

        call transfer_task.transfer as transfer_set_outputs {
            input:
                out_dir = covwwt_novel_mutations_path,
                overwrite = overwite_project_and_set,
                subdirs_to_files = set_subdirs_to_files
        }
    }

    SubdirsToFiles historical_subdirs_to_files = object { subdirs_to_files: [
        (append_new_mutations.historical_full_updated, ""),
        (append_new_mutations.historical_unique_updated, "")
    ]}

    call transfer_task.transfer as transfer_appended_outputs {
        input:
            out_dir = historical_data_path,
            overwrite = overwrite_historical,
            subdirs_to_files = historical_subdirs_to_files
    }

    output {
        Array[File] project_unique_mutations = append_new_mutations.project_unique_mutations
        File historical_full_updated = append_new_mutations.historical_full_updated
        File historical_unique_updated = append_new_mutations.historical_unique_updated
        Array[File]? project_missing_dates = append_new_mutations.project_missing_dates
        File? recurrent_mutations = append_new_mutations.recurrent_mutations
        File? novel_mutations = append_new_mutations.novel_mutations
        String transfer_date_novel_mutations = transfer_appended_outputs.transfer_date
    }
}

task append_new_mutations {

    meta {
        volatile: true
    }

    input {
        Array[String] project_names_array
        Array[File] combined_mutations_array
        Array[String] sites_to_drop
        String today
        File historical_full
        File historical_unique
        File metadata
        File gff_mutations
        File novel_mutations_append_py
    }

    command <<<

        python ~{novel_mutations_append_py} --project_names ~{sep=' ' project_names_array} \
        --combined_mutations_files ~{sep=' ' combined_mutations_array} \
        --historical_full ~{historical_full} \
        --historical_unique ~{historical_unique} \
        --metadata ~{metadata} \
        --gff ~{gff_mutations} \
        --today ~{today} \
        --sites_to_drop ~{sep=' ' sites_to_drop}

    >>>

    output {
        Array[File]? project_missing_dates = glob("*_missing_dates.tsv")
        Array[File] project_unique_mutations = glob("*_unique_mutations.tsv")
        File historical_full_updated = "novel_mutations_historical_full.tsv"
        File historical_unique_updated = "novel_mutations_historical_unique.tsv"        
        File? recurrent_mutations = "recurrent_mutations_~{today}.tsv"
        File? novel_mutations = "novel_mutations_~{today}.tsv"
    }

    runtime {
        docker: "ariannaesmith/py3.10.9-bio:v1"
        memory: "1 GB"
        cpu: 4
        disks: "local-disk 10 SSD"
    }
}
