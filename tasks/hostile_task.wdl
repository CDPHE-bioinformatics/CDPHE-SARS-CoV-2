# Adapted from pull request from Theiagen's PHB repo:
# https://github.com/theiagen/public_health_bioinformatics/pull/256/files

version 1.0

import "../tasks/version_capture_task.wdl"

task hostile {
  input {
    File fastq1
    File? fastq2
    String seq_method

    # genome index may be multiple files like for bowtie2
    Array[File] genome_index
  }

  # make these variables private and constant for now but may move to input to allow customization
  String docker = "quay.io/biocontainers/hostile:1.0.0--pyhdfd78af_0"
  Int disk_size = 100
  Int cpu = 4
  Int mem = 16

  String base_name = basename(basename(basename(fastq1, ".gz"), ".fastq"), ".fq")
  String fastq1_scrubbed_name = base_name + "_scrubbed.fastq.gz"

  # workaround since can't use basename() on an optional file
  String fastq2_scrubbed_name = sub(fastq1_scrubbed_name, "1(?=_scrubbed)", "2")

  command <<<
    # date and version control
    date | tee DATE
    hostile --version | tee VERSION

    # dehost reads based on sequencing method
    if [[ "~{seq_method}" == "OXFORD_NANOPORE" ]]; then
      hostile clean \
        --fastq1 ~{fastq1} \
        --aligner "minimap2" \
        --threads ~{cpu} \
        --index ~{genome_index[0]} \
      | tee decontamination-log.json
      # rename scrubbed fastq
      mv ./*.clean.fastq.gz "~{fastq1_scrubbed_name}"
    else
      hostile clean \
        --fastq1 ~{fastq1} \
        --fastq2 ~{fastq2} \
        --aligner "bowtie2" \
        --threads ~{cpu} \
        --index ~{sub(genome_index[0], ".1.bt2", "")} \
      | tee decontamination-log.json
      # rename scrubbed fastqs
      mv ./*.clean_1.fastq.gz "~{fastq1_scrubbed_name}"
      mv ./*.clean_2.fastq.gz "~{fastq2_scrubbed_name}"
    fi

    # extract the number of removed human reads
    grep '"reads_removed":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS
    grep '"reads_removed_proportion":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS_PROP
  >>>
  output {
    File fastq1_scrubbed = "${fastq1_scrubbed_name}"
    File? fastq2_scrubbed = "${fastq2_scrubbed_name}"
    String human_reads_removed = read_string("HUMANREADS")
    String human_reads_removed_proportion = read_string("HUMANREADS_PROP")
    String hostile_docker = docker

    VersionInfo hostile_version_info = object {
      software: "hostile",
      docker: docker,
      version: read_string("VERSION")
    }
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