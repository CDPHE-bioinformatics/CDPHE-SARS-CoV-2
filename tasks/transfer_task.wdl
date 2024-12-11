version 1.0

# workaround cromwell bug with write_json of Array
# https://github.com/broadinstitute/cromwell/issues/4625
struct FilesToSubdirs {
    Array[Pair[Array[File?], String]] files_to_subdirs
}

task transfer {
    input {
        String out_dir
        Boolean overwrite
        Int cpu = 1
        FilesToSubdirs files_to_subdirs
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        files_to_subdirs_json=~{write_json(files_to_subdirs)}

        if grep -q '|' "$files_to_subdirs_json"; then
            echo "Error: filename or directory cannot contain '|' character" >&2
            exit 1
        fi

        python3 <<CODE

        import json
        import os
        import subprocess
        import sys

        with open("${files_to_subdirs_json}", 'r') as infile:
            if '|' in infile.read():
                sys.exit('Error: filename or directory cannot contain "|" character')
            pairs = json.load(infile)['files_to_subdirs']
            for pair in pairs:
                # filename = pair['left'] if pair['left'] is not None else ''
                # outfile.write(filename + '|' + pair['right'] + '\n')
                filenames = pair['left']
                subdir = pair['right']
                destination = os.path.join(~{outdirpath}, subdir)
                if filenames:
                    command = ['gsutil', '-m', 'cp', filenames, destination]
                    subprocess.run(command)

        CODE

        # # Check if files already exist at the destination
        # declare -a destinations
        # declare -a existing_files
        # while IFS='|' read -r file subdir; do
        #     if [[ -n "$file" ]]; then
        #         filename=$(basename "$file")
        #         destination="~{outdirpath}/${subdir}/${filename}"
        #         destinations+=( "$destination" )
        #     fi
        # done < files_to_subdirs.txt
        # existing_files=( $(gsutil ls "${destinations[@]}") )

        # if [[ ~{overwrite} = true && ${#existing_files[@]} == 0 ]]; then
        #     echo "Warning: overwrite set to true but no files at destination to overwrite" >&2
        # fi
        # if [[ ~{overwrite} = false && ${#existing_files[@]} != 0 ]]; then
        #     echo "Error: overwrite set to false but files exist at destination" >&2
        #     echo "Existing files: ${existing_files[@]}" >&2
        #     exit 1
        # fi

        # while IFS='|' read -r file subdir; do
        #     if [[ -n "$file" ]]; then
        #         if [[ ~{overwrite} = true ]]; then
        #             gsutil -m cp "$file" "~{outdirpath}/${subdir}/"

        #         # Do not clobber in case of TOCTOU race condition
        #         else
        #             gsutil -m cp -n "$file" "~{outdirpath}/${subdir}/"
        #         fi
        #     fi
        # done < files_to_subdirs.txt

        transferdate=$(date)
        echo "$transferdate" | tee TRANSFERDATE

    >>>


    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: cpu
        disks: "local-disk 100 SSD"
    }
}
