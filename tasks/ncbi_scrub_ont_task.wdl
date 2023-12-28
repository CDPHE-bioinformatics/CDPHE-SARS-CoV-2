# This file contains modified code from the ncbi_scrub_se task from Theiagen's
# Public Health Bioinformatics (PHB) repo.
# https://github.com/theiagen/public_health_bioinformatics/blob/e4cf6083c53757f6df972ab7503c1a3a08c008a2/tasks/quality_control/task_ncbi_scrub.wdl
# Task is modified to handle unconcatenated FASTQ files output by ONT GridIon

version 1.0

task ncbi_scrub_ont {
  input {
    Array[File] fastq_files
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE

    total_removed=0

    # dehost reads
    for fastq in ~{sep=' ' fastq_files}; do
        echo $fastq
        echo ${fastq}
        cat $fastq > test.fastq
        ls -l $fastq
        /opt/scrubber/scripts/scrub.sh ${fastq}
        /opt/scrubber/scripts/scrub.sh $fastq
        /opt/scrubber/scripts/scrub.sh test.fastq
        removed=$(/opt/scrubber/scripts/scrub.sh ${fastq} |& tail -n1 | awk -F" " '{print $1}')
        total_removed=$(echo "print(${total_removed} + ${removed})" | python3)
        name=$(basename ${fastq} .fastq)
        mv ./${name}.fastq.clean "${name}_dehosted.fastq"
        gzip "${name}_dehosted.fastq" -c > "${name}_dehosted.fastq.gz"
    done

    echo ${total_removed} > SPOTS_REMOVED
  >>>
  output {
    Array[File] fastq_files_dehosted  = glob("*_dehosted.fastq.gz")
    Int human_spots_removed = read_int("SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}