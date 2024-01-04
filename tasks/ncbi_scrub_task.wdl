# This file copied from Thieagen's Public Health Bioinformatics (PHB) repo
# https://github.com/theiagen/public_health_bioinformatics/blob/e4cf6083c53757f6df972ab7503c1a3a08c008a2/tasks/quality_control/task_ncbi_scrub.wdl

version 1.0

task ncbi_scrub_pe {
  input {
    File fastq1
    File fastq2
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
    Int cpu = 4
  }

  String fastq1_scrubbed_name = select_first([basename(fastq1, ".fastq.gz"), basename(fastq1, ".fastq")]) + "_scrubbed.fastq.gz"
  String fastq2_scrubbed_name = sub(fastq1_scrubbed_name, "1(?=_scrubbed)", "2")

  command <<<
    # date and version control
    date | tee DATE

    # unzip read files as scrub tool does not take in .gz fastq files, and interleave them
    paste <(zcat ~{fastq1} | paste - - - -) <(zcat ~{fastq2} | paste - - - -) | tr '\t' '\n' > interleaved.fastq

    # dehost reads
    /opt/scrubber/scripts/scrub.sh -i interleaved.fastq |& tail -n1 | awk -F" " '{print $1}' > HUMANREADS

    # split interleaved reads and compress files
    paste - - - - - - - - < interleaved.fastq.clean \
      | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > ~{fastq1_scrubbed_name}) \
      | cut -f 5-8 | tr '\t' '\n' | gzip > ~{fastq2_scrubbed_name}
      
  >>>
  output {
    File fastq1_scrubbed = "~{fastq1_scrubbed_name}"
    File fastq2_scrubbed = "~{fastq2_scrubbed_name}"
    String human_reads_removed = read_string("HUMANREADS")
    String ncbi_scrub_docker = docker
  }
  runtime {
      docker: docker
      memory: "8 GB"
      cpu: cpu
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 3
  }
}

task ncbi_scrub_se {
  input {
    File fastq1
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
    Int cpu = 4
  }

  String fastq1_scrubbed_name = select_first([basename(fastq1, ".fastq.gz"), basename(fastq1, ".fastq")]) + "_scrubbed.fastq.gz"

  command <<<
    # date and version control
    date | tee DATE

    # unzip fwd file as scrub tool does not take in .gz fastq files
    if [[ "~{fastq1}" == *.gz ]]
    then
      gunzip -c ~{fastq1} > r1.fastq
      fastq1_unzip=r1.fastq
    else
      fastq1_unzip=~{fastq1}
    fi

    # dehost reads
    /opt/scrubber/scripts/scrub.sh -i ${fastq1_unzip} |& tail -n1 | awk -F" " '{print $1}' > HUMANREADS

    # gzip dehosted reads
    gzip $(basename ${fastq1_unzip}).clean -c > ~{fastq1_scrubbed_name}
  >>>
  output {
    File fastq1_scrubbed = "~{fastq1_scrubbed_name}"
    String human_reads_removed = read_string("HUMANREADS")
    String ncbi_scrub_docker = docker
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}