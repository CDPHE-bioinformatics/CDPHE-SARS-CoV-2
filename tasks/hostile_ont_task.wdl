# Modified from pull request from Theiagen's PHB repo
# https://github.com/theiagen/public_health_bioinformatics/pull/256/files
# Modified to only handle ONT data and takes a array of FASTQ files
# for each sample (as output on the machine), rather than a single FASTQ
# file as in hostile.task.wdl


version 1.0

task hostile_ont {
  input {
    Array[File] fastq_files
    String docker = "quay.io/biocontainers/hostile:0.3.0--pyhdfd78af_0"

    Int disk_size = 100
    Int cpu = 4
    Int mem = 16
  }
  command <<<

    # date and version control
    date | tee DATE
    hostile --version | tee VERSION

    total_removed=0
    total_proportion=0

    # dehost reads based on sequencing method
    # this task currently not used for ONT because task expects
    # concatenated reads. Instead using hostile_ont.wdl
    for fastq in ~{sep=' ' fastq_files}; do
      hostile clean \
        --fastq1 ${fastq} \
        --aligner "minimap2" \
        --threads ~{cpu} > decontamination-log.json
      name=$(basename ${fastq} .fastq)
      mv ./${name}.clean.fastq.gz "${name}_dehosted.fastq.gz"

      # extract the number of removed human reads
      removed=$(grep '"reads_removed":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}')
      proportion=$(grep '"reads_removed_proportion":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}')
      total_removed=$(echo "print(${total_removed} + ${removed})" | python)
      total_proportion=$(echo "print(${total_proportion} + ${proportion})" | python)
    done

    echo ${total_removed} > HUMANREADS
    echo ${total_proportion} > HUMANREADS_PROP
  >>>
  output {
    String hostile_version = read_string("VERSION")
    Array[File] fastq_files_dehosted = glob("*_dehosted.fastq.gz")
    Int? human_reads_removed = read_int("HUMANREADS")
    Float? human_reads_removed_proportion = read_float("HUMANREADS_PROP")
    String hostile_docker = docker
  }
  runtime {
      docker: docker
      memory: "~{mem} GB"
      cpu: cpu
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 3
  }
}