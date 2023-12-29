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
    fastq1_basename=${fastq1_name%%.*}  # get name without any extensions

    # dehost reads based on sequencing method
    if [[ "~{seq_method}" == "OXFORD_NANOPORE" ]]; then
      hostile clean \
        --fastq1 ~{fastq1} \
        --aligner "minimap2" \
        --threads ~{cpu} | tee decontamination-log.json

      # rename scrubbed fastq
      mv ./*.clean.fastq.gz "${fastq1_basename}_scrubbed.fastq.gz"
    else
      hostile clean \
        --fastq1 ~{fastq1} \
        --fastq2 ~{fastq2} \
        --aligner "bowtie2" \
        --threads ~{cpu} | tee decontamination-log.json

      # rename scrubbed fastqs
      mv ./*.clean_1.fastq.gz "${fastq1_basename}_scrubbed.fastq.gz"
      fastq2_name=~{fastq2}
      fastq2_basename=${fastq2_name%%.*}
      mv ./*.clean_2.fastq.gz "${fastq1_basename}_scrubbed.fastq.gz"
    fi

    # extract the number of removed human reads
    grep '"reads_removed":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS
    grep '"reads_removed_proportion":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS_PROP
  >>>
  output {
    String hostile_version = read_string("VERSION")

    # Illumina has R1/R2 so make sure to grab R1 file for fastq1_scrubbed
    File fastq1_scrubbed = if seq_method == "OXFORD_NANOPORE" then glob("*_scrubbed.fastq.gz")[0] else glob("*_R1_scrubbed.fastq.gz")[0]

    File? fastq2_scrubbed = glob("*_R2_scrubbed.fastq.gz")[0]
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