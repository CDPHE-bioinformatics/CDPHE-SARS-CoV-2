version 1.0

task transfer {
    input {
        String out_dir
        Boolean overwrite

        # Map[File, String] caused random failures on Terra
        # Use [select_first([file, ""]), "subdir"] for File? type
        Array[Array[String]] file_to_subdir
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<

        file_to_subdir_tsv=~{write_tsv(file_to_subdir)}

        # Check if files already exist at the destination
        declare -a destinations
        declare -a existing_files
        while IFS=$'\t' read -r file subdir; do
            if [[ -n "$file" ]]; then
                filename=$(basename "$file")
                destination="~{outdirpath}/${subdir}/${filename}"
                destinations+=( "$destination" )
            fi
        done < "$file_to_subdir_tsv"
        existing_files=( "$(gsutil ls "${destinations[@]}")" )

        # existing_files will contain one element even if 'empty'
        if [[ ~{overwrite} = true && ${#existing_files[@]} == 1 ]]; then
            echo "Error: overwrite set to true but no files at destination to overwrite" >&2
            exit 1
        fi
        if [[ ~{overwrite} = false && ${#existing_files[@]} != 1 ]]; then
            echo "Error: overwrite set to false but files exist at destination" >&2
            echo "Existing files: ${existing_files[@]}" >&2
            exit 1
        fi

        while IFS=$'\t' read -r file subdir; do
            if [[ -n "$file" ]]; then
                gsutil cp "$file" "~{outdirpath}/${subdir}/"
            fi
        done < "$file_to_subdir_tsv"

        transferdate=$(date)
        echo "$transferdate" | tee TRANSFERDATE

    >>>


    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}
