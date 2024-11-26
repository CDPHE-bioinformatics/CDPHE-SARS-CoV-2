version 1.0

# workaround cromwell bug with write_json of Array
# https://github.com/broadinstitute/cromwell/issues/4625
struct FilesToSubdirs {
    Array[Pair[File?, String]] files_to_subdirs
}

task transfer {
    input {
        String out_dir
        Boolean overwrite
 
        FilesToSubdirs files_to_subdirs
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<

        file_to_subdir_json=~{write_json(files_to_subdirs)}
        cat $file_to_subdir_json

        # if grep -q '|' "$file_to_subdir_tsv"; then
        #     echo "Error: filename or directory cannot contain '|' character" >&2
        #     exit 1
        # fi

        # # Check if files already exist at the destination
        # declare -a destinations
        # declare -a existing_files
        # while IFS='|' read -r file subdir; do
        #     if [[ -n "$file" ]]; then
        #         filename=$(basename "$file")
        #         destination="~{outdirpath}/${subdir}/${filename}"
        #         destinations+=( "$destination" )
        #     fi
        # done < <(tr '\t' '|' < "$file_to_subdir_tsv")
        # existing_files=( $(gsutil ls "${destinations[@]}") )

        # if [[ ~{overwrite} = true && ${#existing_files[@]} == 0 ]]; then
        #     echo "Error: overwrite set to true but no files at destination to overwrite" >&2
        #     exit 1
        # fi
        # if [[ ~{overwrite} = false && ${#existing_files[@]} != 0 ]]; then
        #     echo "Error: overwrite set to false but files exist at destination" >&2
        #     echo "Existing files: ${existing_files[@]}" >&2
        #     exit 1
        # fi

        # while IFS='|' read -r file subdir; do
        #     if [[ -n "$file" ]]; then
        #         gsutil cp "$file" "~{outdirpath}/${subdir}/"
        #     fi
        # done < <(tr '\t' '|' < "$file_to_subdir_tsv")

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
