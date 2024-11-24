version 1.0

task transfer {
    input {
        String outdirpath
        Array[File] file_to_subdir
    }


    command <<<
        for file in ~{sep=' ' file_to_subdir}; do
            gsutil -m cp "$file" "~{outdirpath}/all_files/"
        done
    
        transferdate=$(date)
        echo "$transferdate" | tee TRANSFERDATE
    >>>


    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
