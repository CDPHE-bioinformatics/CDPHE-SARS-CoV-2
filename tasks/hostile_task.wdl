# Taken from pull request from Theiagen's PHB repo:
# https://github.com/theiagen/public_health_bioinformatics/pull/256/files

version 1.0

task hostile {
  input {
    File fastq1
    File? fastq2
    String docker = "quay.io/biocontainers/hostile:0.3.0--pyhdfd78af_0"
    String seq_method

    Int disk_size = 100
    Int cpu = 4
    Int mem = 16
  }

  String fastq1_scrubbed_name = select_first([basename(fastq1, ".fastq.gz"), basename(fastq1, ".fastq")]) + "_scrubbed.fastq.gz"

  # workaround to get basename() for optional file
  String fastq2_name = if defined(fastq2) then fastq2 else ""
  String fastq2_scrubbed_name = if defined(fastq2) then select_first([basename(fastq2_name, ".fastq.gz"), basename(fastq2_name, ".fastq")]) + "_scrubbed.fastq.gz" else ""

  command <<<
    # date and version control
    date | tee DATE
    hostile --version | tee VERSION

    # empty file for optional output fastq2
    touch NULL

    # dehost reads based on sequencing method
    if [[ "~{seq_method}" == "OXFORD_NANOPORE" ]]; then
      hostile clean \
        --fastq1 ~{fastq1} \
        --aligner "minimap2" \
        --threads ~{cpu} | tee decontamination-log.json

      # rename scrubbed fastq
      mv ./*.clean.fastq.gz "~{fastq1_scrubbed_name}"
    else
      hostile clean \
        --fastq1 ~{fastq1} \
        --fastq2 ~{fastq2} \
        --aligner "bowtie2" \
        --threads ~{cpu} | tee decontamination-log.json

      # rename scrubbed fastqs
      mv ./*.clean_1.fastq.gz "~{fastq1_scrubbed_name}"
      mv ./*.clean_2.fastq.gz "~{fastq2_scrubbed_name}"
    fi

    # extract the number of removed human reads
    grep '"reads_removed":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS
    grep '"reads_removed_proportion":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS_PROP
  >>>
  output {
    String hostile_version = read_string("VERSION")
    File fastq1_scrubbed = "${fastq1_scrubbed_name}"
    File? fastq2_scrubbed = if defined(fastq2) then "${fastq2_scrubbed_name}" else "NULL"
    String human_reads_removed = read_string("HUMANREADS")
    String human_reads_removed_proportion = read_string("HUMANREADS_PROP")
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