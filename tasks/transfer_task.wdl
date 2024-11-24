version 1.0

task transfer {
    input {
        String outdirpath
        Array[Array[String]] file_to_subdir
    }


    command <<<
        while IFS=$'\t' read -r file subdir; do
            gsutil -m cp "$file" "~{outdirpath}/${subdir}/"
        done < ~{write_tsv(file_to_subdir)}
    
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
