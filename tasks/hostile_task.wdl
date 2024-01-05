# Taken from pull request from Theiagen's PHB repo:
# https://github.com/theiagen/public_health_bioinformatics/pull/256/files

version 1.0

task hostile {
  input {
    File fastq1
    File? fastq2
    String seq_method
    File genome_index

    String docker = "quay.io/biocontainers/hostile:0.3.0--pyhdfd78af_0"
    Int disk_size = 100
    Int cpu = 4
    Int mem = 16
  }

  String fastq1_scrubbed_name = select_first([basename(fastq1, ".fastq.gz"), basename(fastq1, ".fastq")]) + "_scrubbed.fastq.gz"
  String fastq2_scrubbed_name = sub(fastq1_scrubbed_name, "1(?=_scrubbed)", "2")

  command <<<
    # date and version control
    date | tee DATE
    hostile --version | tee VERSION

    echo "CPU ~{cpu}"

    # dehost reads based on sequencing method
    if [[ "~{seq_method}" == "OXFORD_NANOPORE" ]]; then
      hostile clean \
        --fastq1 ~{fastq1} \
        --aligner "minimap2" \
        --threads ~{cpu} \
        --index ~{genome_index} | tee decontamination-log.json
      # rename scrubbed fastq
      mv ./*.clean.fastq.gz "~{fastq1_scrubbed_name}"
    else
      tar xvf ~{genome_index}
      hostile clean \
        --fastq1 ~{fastq1} \
        --fastq2 ~{fastq2} \
        --aligner "bowtie2" \
        --threads ~{cpu} \
        --index ~{basename(genome_index, ".tar")} | tee decontamination-log.json
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
    File? fastq2_scrubbed = "${fastq2_scrubbed_name}"
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