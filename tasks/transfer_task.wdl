version 1.0

# workaround cromwell bug with write_json of Array
# https://github.com/broadinstitute/cromwell/issues/4625
struct SubdirsToFiles {
    Array[Pair[String, Array[File?]]] subdirs_to_files
}

task transfer {
    input {
        String out_dir
        Boolean overwrite
        Int cpu = 1
        SubdirsToFiles subdirs_to_files
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        subdirs_to_files_json=~{write_json(subdirs_to_files)}


        python3 <<CODE

        import json
        import os
        import subprocess
        import sys

        from collections import defaultdict

        overwrite = True if '~{overwrite}' == 'true' else False

        # map destination folders to corresponding source files
        with open('${subdirs_to_files_json}', 'r') as infile:
            pairs = json.load(infile)['subdirs_to_files']
            destinations = []
            destinations_dict = defaultdict(list)
            for pair in pairs:
                subdir = pair['left']
                sources = pair['right'] if not all(s is None for s in pair['right']) else []
                destination = os.path.join('~{outdirpath}', subdir)
                filenames = [os.path.basename(s) for s in sources]
                destination_files = [os.path.join(destination, f) for f in filenames]
                destinations.extend(destination_files)
                destinations_dict[destination].extend(sources)

        # check if files already exist at the destination
        command = ['gsutil', 'ls', *destinations]
        output = subprocess.run(command, capture_output=True, text=True)
        existing_files = output.stdout.splitlines()
        if overwrite and not existing_files:
            print('Warning: overwrite set to true but no files at destination to overwrite', file=sys.stderr)
        if not overwrite and existing_files:
            sys.exit(f'Error: overwrite set to false but files exist at the destination: {existing_files}')

        # copy files
        for destination in destinations_dict:
            sources = destinations_dict[destination]
            if filenames:
                clobber = '-n' if not overwrite else ''
                command = ['gsutil', '-m', 'cp', clobber, *sources, destination]
                subprocess.run(command)

        CODE

        transferdate=$(date)
        echo "$transferdate" | tee TRANSFERDATE

    >>>


    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.2"
        memory: "2 GB"
        cpu: cpu
        disks: "local-disk 100 SSD"
    }
}
