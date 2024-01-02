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
  command <<<
    # date and version control
    date | tee DATE
    hostile --version | tee VERSION

    fastq1_name=~{fastq1}
    fastq1_scrubbed_name="${fastq1_name%%.*}_scrubbed_fastq.gz"
    echo ${fastq1_scrubbed_name} > FASTQ1_SCRUBBED_NAME

    # dehost reads based on sequencing method
    if [[ "~{seq_method}" == "OXFORD_NANOPORE" ]]; then
      hostile clean \
        --fastq1 ~{fastq1} \
        --aligner "minimap2" \
        --threads ~{cpu} | tee decontamination-log.json

      # rename scrubbed fastq
      mv ./*.clean.fastq.gz "${fastq1_scrubbed_name}"
    else
      hostile clean \
        --fastq1 ~{fastq1} \
        --fastq2 ~{fastq2} \
        --aligner "bowtie2" \
        --threads ~{cpu} | tee decontamination-log.json

      # rename scrubbed fastqs
      mv ./*.clean_1.fastq.gz "${fastq1_scrubbed_name}"
      fastq2_name=~{fastq2}
      fastq2_scrubbed_name="${fastq2_name%%.*}_scrubbed_fastq.gz"
      echo ${fastq2_scrubbed_name} > FASTQ2_SCRUBBED_NAME
      mv ./*.clean_2.fastq.gz "${fastq2_scrubbed_name}"
    fi

    # extract the number of removed human reads
    grep '"reads_removed":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS
    grep '"reads_removed_proportion":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS_PROP
  >>>
  output {
    String hostile_version = read_string("VERSION")
    File fastq1_scrubbed = basename(read_string("FASTQ1_SCRUBBED_NAME"))
    File? fastq2_scrubbed = basename(read_string("FASTQ2_SCRUBBED_NAME"))
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