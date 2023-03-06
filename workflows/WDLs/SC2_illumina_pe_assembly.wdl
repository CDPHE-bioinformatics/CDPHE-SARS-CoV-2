version 1.0

workflow SC2_illumina_pe_assembly {

    input {
        String  sample_id
        File    fastq_1
        File    fastq_2
        File    primer_bed
        File    adapters_and_contaminants
        File    covid_genome
        File    covid_gff
        File    preprocess_python_script
        File    s_gene_amplicons
    }

    call seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_id = sample_id,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    call fastqc as fastqc_raw {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call fastqc as fastqc_cleaned {
        input:
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call align_reads {
        input:
            sample_id = sample_id,
            ref = covid_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call ivar_trim {
        input:
            sample_id = sample_id,
            primers = primer_bed,
            bam = align_reads.out_bam
    }

    call ivar_var {
        input:
            sample_id = sample_id,
            ref = covid_genome,
            gff = covid_gff,
            bam = ivar_trim.trimsort_bam
    }

    call ivar_consensus {
        input:
            sample_id = sample_id,
            ref = covid_genome,
            bam = ivar_trim.trimsort_bam
    }

    call bam_stats {
        input:
            sample_id = sample_id,
            bam = ivar_trim.trimsort_bam,
            bai = ivar_trim.trimsort_bamindex,
            s_gene_amplicons = s_gene_amplicons
    }

    call rename_fasta {
        input:
            sample_id = sample_id,
            fasta = ivar_consensus.consensus_out
    }

    call calc_percent_cvg {
        input:
            sample_id = sample_id,
            fasta = rename_fasta.renamed_consensus,
            preprocess_python_script = preprocess_python_script
    }

    output {
        File filtered_reads_1 = seqyclean.cleaned_1
        File filtered_reads_2 = seqyclean.cleaned_2
        File seqyclean_summary = seqyclean.seqyclean_summary
        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw1_zip = fastqc_raw.fastqc1_zip
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_raw2_zip = fastqc_raw.fastqc2_zip
        File fastqc_clean1_html = fastqc_cleaned.fastqc1_html
        File fastqc_clean1_zip = fastqc_cleaned.fastqc1_zip
        File fastqc_clean2_html = fastqc_cleaned.fastqc2_html
        File fastqc_clean2_zip = fastqc_cleaned.fastqc2_zip
        File out_bam = align_reads.out_bam
        File out_bamindex = align_reads.out_bamindex
        File trim_bam = ivar_trim.trim_bam
        File trimsort_bam = ivar_trim.trimsort_bam
        File trimsort_bamindex = ivar_trim.trimsort_bamindex
        File variants = ivar_var.var_out
        File consensus = ivar_consensus.consensus_out
        File flagstat_out = bam_stats.flagstat_out
        File stats_out = bam_stats.stats_out
        File covhist_out = bam_stats.covhist_out
        File cov_out = bam_stats.cov_out
        File cov_s_gene_out = bam_stats.cov_s_gene_out
        File cov_s_gene_amplicons_out = bam_stats.cov_s_gene_amplicons_out
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        String assembler_version = align_reads.assembler_version
    }
}

task seqyclean {
    input {
        File contam
        String sample_id
        File fastq_1
        File fastq_2
    }

    command {

        seqyclean -minlen 70 -qual 30 30 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contam} -o ${sample_id}_clean

    }

    output {

        File cleaned_1 = "${sample_id}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_id}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_id}_clean_SummaryStatistics.tsv"

    }

    runtime {
        cpu:    2
        memory:    "6 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/seqyclean:1.10.09"
    }
}

task fastqc {
    input {

        File fastq_1
        File fastq_2
    }

    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command {

        fastqc --outdir $PWD ${fastq_1} ${fastq_2}

    }

    output {

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"

    }

    runtime {
        cpu:    1
        memory:    "2 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/fastqc:0.11.9"
    }
}

task align_reads {

    input {

        File fastq_1
        File fastq_2
        File ref
        String sample_id
    }

    command {

        echo bwa 0.7.17-r1188 > VERSION
        
        bwa index -p reference.fasta -a is ${ref}
        bwa mem -t 2 reference.fasta ${fastq_1} ${fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./${sample_id}_aln.sorted.bam
        samtools index ./${sample_id}_aln.sorted.bam

    }

    output {

        File out_bam = "${sample_id}_aln.sorted.bam"
        File out_bamindex = "${sample_id}_aln.sorted.bam.bai"
        String assembler_version = read_string("VERSION")

    }

    runtime {
        cpu:    2
        memory:    "12 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "broadinstitute/viral-core:latest"
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_id
    }

    command {

        ivar trim -e -i ${bam} -b ${primers} -p ${sample_id}_trim.bam
        samtools sort ${sample_id}_trim.bam -o ${sample_id}_trim.sort.bam
        samtools index ${sample_id}_trim.sort.bam

    }

    output {

        File trim_bam = "${sample_id}_trim.bam"
        File trimsort_bam = "${sample_id}_trim.sort.bam"
        File trimsort_bamindex = "${sample_id}_trim.sort.bam.bai"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_var {

    input {

        String sample_id
        File ref
        File gff
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_id}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}

    }

    output {

        File var_out = "${sample_id}_variants.tsv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_consensus {

    input {
        String sample_id
        File ref
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar consensus -p ${sample_id}_consensus -q 20 -t 0.6 -m 10

    }

    output {
        File consensus_out = "${sample_id}_consensus.fa"
    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task bam_stats {

    input {
        String sample_id
        File bam
        File bai
        File s_gene_amplicons
    }

    command <<<

        samtools flagstat ~{bam} > ~{sample_id}_flagstat.txt
        samtools stats ~{bam} > ~{sample_id}_stats.txt
        samtools coverage -m -o ~{sample_id}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_id}_coverage.txt ~{bam}

        # Calculate depth of coverage over entire S gene
        echo "Calculating overall S gene depth"
        samtools coverage --region MN908947.3:21,563-25,384 \
            -o ~{sample_id}_S_gene_coverage.txt ~{bam}

        # Calculate depth of coverage over S gene amplicon regions (excludes overlapping regions with adjacent amplicons)
        echo "calculating depths for ~{s_gene_amplicons}"
        {
            s_gene_depths="~{sample_id}_S_gene_depths.tsv"

            # write header line to s_gene_depths output file
            read header
            echo -e "${header}\tdepth" | tee $s_gene_depths

            # write amplicon info and depths to output file
            IFS=$'\t'
            while read amplicon coords description; do
                line=$(echo -e "${amplicon}\t${coords}\t${description}\t")

                # extract mean amplicon depth from samtools coverage output
                line+=$(samtools coverage --region MN908947.3:${coords} ~{bam} \
                            | cut -f 7 | sed '2q;d')

                echo -e "$line" | tee -a $s_gene_depths
            done
        } < ~{s_gene_amplicons}

    >>>

    output {

        File flagstat_out  = "${sample_id}_flagstat.txt"
        File stats_out  = "${sample_id}_stats.txt"
        File covhist_out  = "${sample_id}_coverage_hist.txt"
        File cov_out  = "${sample_id}_coverage.txt"
        File cov_s_gene_out = "${sample_id}_S_gene_coverage.txt"
        File cov_s_gene_amplicons_out = "${sample_id}_S_gene_depths.tsv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.16"
    }
}

task rename_fasta {

    input {

        String sample_id
        File fasta
    }

    command <<<

        sed 's/>.*/>CO-CDPHE-~{sample_id}/' ~{fasta} > ~{sample_id}_consensus_renamed.fa

    >>>

    output {

        File renamed_consensus  = "${sample_id}_consensus_renamed.fa"

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
        String sample_id
        File preprocess_python_script

    }

    command {
        python ~{preprocess_python_script} \
          --sample_id ~{sample_id} \
          --fasta_file ~{fasta}
      }
    output {

      File percent_cvg_csv  = "${sample_id}_consensus_cvg_stats.csv"

    }

    runtime {

      docker: "mchether/py3-bio:v1"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }

}