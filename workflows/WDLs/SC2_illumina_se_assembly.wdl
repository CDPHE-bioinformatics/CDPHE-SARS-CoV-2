version 1.0

workflow SC2_illumina_se_assembly {

    input {
        String    sample_name
        File    fastq
        File    primer_bed
        File    covid_genome
        File    covid_gff
        File    calc_percent_coverage_py
    }

    call trimmomatic {
        input:
            sample_name = sample_name,
            fastq = fastq
    }
    
    call bbduk_clean {
        input:
            sample_name = sample_name,
            fastq = trimmomatic.trimmed
    }
    
    call fastqc as fastqc_raw {
        input:
           fastq = fastq
    }

    call fastqc as fastqc_cleaned {
        input:
            fastq = bbduk_clean.cleaned
    }

    call align_reads {
        input:
            sample_name = sample_name,
            ref = covid_genome,
            fastq = bbduk_clean.cleaned
    }

    call ivar_trim {
        input:
            sample_name = sample_name,
            primers = primer_bed,
            bam = align_reads.out_bam
    }

    call ivar_var {
        input:
            sample_name = sample_name,
            ref = covid_genome,
            gff = covid_gff,
            bam = ivar_trim.trimsort_bam
    }

    call ivar_consensus {
        input:
            sample_name = sample_name,
            ref = covid_genome,
            bam = ivar_trim.trimsort_bam
    }

    call bam_stats {
        input:
            sample_name = sample_name,
            bam = ivar_trim.trimsort_bam
    }

    call rename_fasta {
        input:
            sample_name = sample_name,
            fasta = ivar_consensus.consensus_out
    }
    
    call calc_percent_cvg {
        input:
            sample_name = sample_name,
            fasta = rename_fasta.renamed_consensus,
            calc_percent_coverage_py = calc_percent_coverage_py
    }
    
    output {
        File trimmed_reads = trimmomatic.trimmed
        File trim_stats = trimmomatic.trim_stats
        File filtered_reads = bbduk_clean.cleaned
        File adapter_stats = bbduk_clean.adapter_stats
        File PhiX_stats = bbduk_clean.PhiX_stats
        File fastqc_raw_html = fastqc_raw.fastqc_html
        File fastqc_raw_zip = fastqc_raw.fastqc_zip
        File fastqc_clean_html = fastqc_cleaned.fastqc_html
        File fastqc_clean_zip = fastqc_cleaned.fastqc_zip
        File out_bam = align_reads.out_bam
        File trim_bam = ivar_trim.trim_bam
        File trimsort_bam = ivar_trim.trimsort_bam
        File trimsort_bamindex = ivar_trim.trimsort_bamindex
        File variants = ivar_var.var_out
        File consensus = ivar_consensus.consensus_out
        File flagstat_out = bam_stats.flagstat_out
        File stats_out = bam_stats.stats_out
        File covhist_out = bam_stats.covhist_out
        File cov_out = bam_stats.cov_out
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        String assembler_version = align_reads.assembler_version
    }
}

task trimmomatic {
    input {
        String sample_name
        File fastq
    }

    command <<<

        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 ~{fastq} ~{sample_name}_trimmed.fastq.gz SLIDINGWINDOW:4:30 MINLEN:25 > ~{sample_name}.trim.stats.txt

    >>>

    output {

        File trimmed = "${sample_name}_trimmed.fastq.gz"
        File trim_stats = "${sample_name}.trim.stats.txt"

    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk 10 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/trimmomatic:0.39"
    }
}

task bbduk_clean {
    input {
        String sample_name
        File fastq
    }

    command <<<

        bbduk.sh -Xmx"8g" in1=~{fastq} out1=~{sample_name}.rmadpt.fastq.gz ref=/bbmap/resources/adapters.fa stats=~{sample_name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
        bbduk.sh -Xmx"8g" in1=~{sample_name}.rmadpt.fastq.gz out1=~{sample_name}_cleaned.fastq.gz outm=~{sample_name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz minlength=25 k=31 hdist=1 stats=~{sample_name}.phix.stats.txt

    >>>

    output {

        File cleaned = "${sample_name}_cleaned.fastq.gz"
        File adapter_stats = "${sample_name}.adapters.stats.txt"
        File PhiX_stats = "${sample_name}.phix.stats.txt"

    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk 100 SSD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/bbtools:38.76"
    }
}

task fastqc {
    input {

        File fastq
    }
    
    String fastq_name = basename(basename(basename(fastq, ".gz"), ".fastq"), ".fq")

    command {

        fastqc --outdir $PWD ${fastq}

    }

    output {

        File fastqc_html = "${fastq_name}_fastqc.html"
        File fastqc_zip = "${fastq_name}_fastqc.zip"

    }

    runtime {
        cpu:    1
        memory:    "2 GiB"
        disks:    "local-disk 200 SSD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/fastqc:0.11.9"
    }
}

task align_reads {

    input {

        File fastq
        File ref
        String sample_name
    }
    
    Int disk_size = 3 * ceil(size(fastq, "GB"))

    command <<<
       
        echo bwa 0.7.17-r1188 > VERSION
        
        bwa index -p reference.fasta -a is ~{ref}
        bwa mem -t 2 reference.fasta ~{fastq} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ~{sample_name}_aln.sorted.bam
        samtools index ~{sample_name}_aln.sorted.bam

    >>>

    output {

        File out_bam = "${sample_name}_aln.sorted.bam"
        String assembler_version = read_string("VERSION")

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 200 SSD"
        bootDiskSizeGb:    100
        preemptible:    0
        maxRetries:    0
        docker:    "broadinstitute/viral-core:latest"
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_name
    }

    command {

        ivar trim -e -i ${bam} -b ${primers} -p ${sample_name}_trim.bam
        samtools sort ${sample_name}_trim.bam -o ${sample_name}_trim.sort.bam
        samtools index ${sample_name}_trim.sort.bam

    }

    output {

        File trim_bam = "${sample_name}_trim.bam"
        File trimsort_bam = "${sample_name}_trim.sort.bam"
        File trimsort_bamindex = "${sample_name}_trim.sort.bam.bai"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 200 SSD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_var {

    input {

        String sample_name
        File ref
        File gff
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_name}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}

    }

    output {

        File var_out = "${sample_name}_variants.tsv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 200 SSD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_consensus {

    input {

        String sample_name
        File ref
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar consensus -p ${sample_name}_consensus -q 20 -t 0.6 -m 10

    }

    output {

        File consensus_out = "${sample_name}_consensus.fa"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 10 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task bam_stats {

    input {

        String sample_name
        File bam
    }

    command {

        samtools flagstat ${bam} > ${sample_name}_flagstat.txt
        samtools stats ${bam} > ${sample_name}_stats.txt
        samtools coverage -m -o ${sample_name}_coverage_hist.txt ${bam}
        samtools coverage -o ${sample_name}_coverage.txt ${bam}

    }

    output {

        File flagstat_out  = "${sample_name}_flagstat.txt"
        File stats_out  = "${sample_name}_stats.txt"
        File covhist_out  = "${sample_name}_coverage_hist.txt"
        File cov_out  = "${sample_name}_coverage.txt"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 10 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.10"
    }
}

task rename_fasta {

    input {

        String sample_name
        File fasta
    }

    command <<<

        sed 's/>.*/>CO-CDPHE-~{sample_name}/' ~{fasta} > ~{sample_name}_consensus_renamed.fa

    >>>

    output {

        File renamed_consensus  = "${sample_name}_consensus_renamed.fa"

    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task calc_percent_cvg {

    input {

        File fasta
        String sample_name
        File calc_percent_coverage_py

    }

    command {
        python ~{calc_percent_coverage_py} \
          --sample_name ~{sample_name} \
          --fasta_file ~{fasta}
      }
    output {

      File percent_cvg_csv  = "${sample_name}_consensus_cvg_stats.csv"

    }

    runtime {

      docker: "mchether/py3-bio:v1"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }
}
