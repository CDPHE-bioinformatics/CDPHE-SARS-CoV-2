# Taken from pull request from Theiagen's PHB repo:
# https://github.com/theiagen/public_health_bioinformatics/pull/256/files

version 1.0

task hostile {
  input {
    File read1
    File? read2
    String samplename
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

    # dehost reads based on sequencing method
    # this task currently not used for ONT because task expects
    # concatenated reads. Instead using hostile_ont.wdl
    if [[ "~{seq_method}" == "OXFORD_NANOPORE" ]]; then
      hostile clean \
        --fastq1 ~{read1} \
        --aligner "minimap2" \
        --threads ~{cpu} > decontamination-log.json

      mv ./*.clean.fastq.gz "~{samplename}_R1_dehosted.fastq.gz"
    else
      hostile clean \
        --fastq1 ~{read1} \
        --fastq2 ~{read2} \
        --aligner "bowtie2" \
        --threads ~{cpu} > decontamination-log.json

      mv ./*.clean_1.fastq.gz "~{samplename}_R1_dehosted.fastq.gz"
      mv ./*.clean_2.fastq.gz "~{samplename}_R2_dehosted.fastq.gz"
    fi
    # extract the number of removed human reads
    grep '"reads_removed":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS
    grep '"reads_removed_proportion":' ./decontamination-log.json | awk -F': ' '{print $2}' | awk -F',' '{print $1}' > HUMANREADS_PROP
  >>>
  output {
    String hostile_version = read_string("VERSION")
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    File? read2_dehosted = "~{samplename}_R2_dehosted.fastq.gz"
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