version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture
import "../tasks/hostile_task.wdl" as hostile_task
import "../tasks/transfer_task.wdl" as transfer_task

workflow SC2_illumina_pe_assembly {

    input {
        String    sample_name
        File    fastq_1
        File    fastq_2
        File    primer_bed
        File    adapters_and_contaminants
        File    covid_genome
        File    covid_gff
        Boolean scrub_reads
        Array[File]? scrub_genome_index
        String  project_name
        String out_dir
        Boolean overwrite

        # python scripts
        File    calc_percent_coverage_py
        File    s_gene_amplicons
        File    version_capture_py
    }

    if (scrub_reads) {
        call hostile_task.hostile as hostile {
            input:
                fastq1 = fastq_1,
                fastq2 = fastq_2,
                genome_index = select_first([scrub_genome_index]),
                seq_method = "ILLUMINA"
        }
    }

    call seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_name = sample_name,
            fastq_1 = select_first([hostile.fastq1_scrubbed, fastq_1]),
            fastq_2 = select_first([hostile.fastq2_scrubbed, fastq_2])
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
            sample_name = sample_name,
            ref = covid_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
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
            bam = ivar_trim.trimsort_bam,
            bai = ivar_trim.trimsort_bamindex,
            s_gene_amplicons = s_gene_amplicons,
            primer_bed = primer_bed
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

    call version_capture.workflow_version_capture  as workflow_version_capture{
        input:
    }

    Array[VersionInfo] version_array = [
        seqyclean.seqyclean_version_info,
        fastqc_cleaned.fastqc_version_info,
        align_reads.bwa_version_info,
        align_reads.samtools_version_info,
        ivar_consensus.ivar_version_info,
        ivar_consensus.samtools_version_info,
        bam_stats.samtools_version_info
    ]
    if (scrub_reads) {
        Array[VersionInfo] version_array_with_hostile = flatten([version_array, select_all([hostile.hostile_version_info])])
    }
    call version_capture.task_version_capture as task_version_capture {
        input:
            version_array = select_first([version_array_with_hostile, version_array]),
            workflow_name = "SC2_illumina_pe_assembly",
            workflow_version_path = workflow_version_capture.workflow_version_path,
            project_name = project_name,
            analysis_date = workflow_version_capture.analysis_date,
            version_capture_py = version_capture_py
    }
    
    FilesToSubdirs files_to_subdirs = object { files_to_subdirs: [
            (hostile.fastq1_scrubbed, "fastq_scrubbed"),
            (hostile.fastq2_scrubbed, "fastq_scrubbed"),
            (seqyclean.seqyclean_summary, "seqyclean"),
            (fastqc_raw.fastqc1_html, "fastqc"),
            (fastqc_raw.fastqc1_zip, "fastqc"),
            (fastqc_raw.fastqc2_html, "fastqc"),
            (fastqc_raw.fastqc2_zip, "fastqc"),
            (fastqc_cleaned.fastqc1_html, "fastqc"),
            (fastqc_cleaned.fastqc1_zip, "fastqc"),
            (fastqc_cleaned.fastqc2_html, "fastqc"),
            (fastqc_cleaned.fastqc2_zip, "fastqc"),
            (ivar_trim.trimsort_bam, "alignments"),
            (ivar_trim.trimsort_bamindex, "alignments"),
            (ivar_var.var_out, "variants"),
            (bam_stats.flagstat_out, "bam_stats"),
            (bam_stats.stats_out, "bam_stats"),
            (bam_stats.covhist_out, "bam_stats"),
            (bam_stats.cov_out, "bam_stats"),
            (bam_stats.depth_out, "bam_stats"),
            (bam_stats.cov_s_gene_out, "bam_stats"),
            (bam_stats.cov_s_gene_amplicons_out, "bam_stats"),
            (rename_fasta.renamed_consensus, "assemblies"),
            (task_version_capture.version_capture_file, "sample_version_capture")
    ]}

    call transfer_task.transfer {
        input:
            out_dir = out_dir,
            overwrite = overwrite,
            files_to_subdirs = files_to_subdirs
    }

    output {
        Int? human_reads_removed = hostile.human_reads_removed
        Float? human_reads_removed_proportion = hostile.human_reads_removed_proportion
        File? fastq1_scrubbed = hostile.fastq1_scrubbed
        File? fastq2_scrubbed = hostile.fastq2_scrubbed
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
        File depth_out = bam_stats.depth_out
        File cov_s_gene_out = bam_stats.cov_s_gene_out
        File cov_s_gene_amplicons_out = bam_stats.cov_s_gene_amplicons_out
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv

        File version_capture_illumina_pe_assembly = task_version_capture.version_capture_file
        String transfer_date_assembly = transfer.transfer_date

    }
}

task seqyclean {
    input {
        File contam
        String sample_name
        File fastq_1
        File fastq_2
    }

    String docker = "staphb/seqyclean:1.10.09"

    command <<<

        seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_1} -2 ~{fastq_2} -c ~{contam} -o ~{sample_name}_clean

        # grab seqyclean version 
        seqyclean -h | awk '/Version/ {print $2}' | tee VERSION
    >>>

    output {

        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
        
        VersionInfo seqyclean_version_info = object {
            software: "seqyclean",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    2
        memory:    "6 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task fastqc {
    input {

        File fastq_1
        File fastq_2
    }

    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    String docker = "staphb/fastqc:0.11.9"

    command <<<

        fastqc --outdir $PWD ~{fastq_1} ~{fastq_2}

        # grab version 
        fastqc --version | awk '/FastQC/ {print $2}' | tee VERSION

    >>>

    output {

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"

        VersionInfo fastqc_version_info = object {
            software: "fastqc",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    1
        memory:    "2 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task align_reads {

    input {

        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    String docker = "quay.io/broadinstitute/viral-core:2.2.3"

    command <<<

        # echo bwa 0.7.17-r1188 > VERSION
        # grab version bwa and samtools versions
        bwa 2>&1 | awk '/Version/{print $2}' | tee VERSION_bwa
        samtools --version | awk '/samtools/ {print $2}' |tee VERSION_samtools
        
        bwa index -p reference.fasta -a is ~{ref}
        bwa mem -t 2 reference.fasta ~{fastq_1} ~{fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./~{sample_name}_aln.sorted.bam
        samtools index ./~{sample_name}_aln.sorted.bam

    >>>

    output {

        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"

        VersionInfo bwa_version_info = object {
            software: "bwa",
            docker: docker,
            version: read_string("VERSION_bwa")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_samtools")
        }
    }

    runtime {
        cpu:    2
        memory:    "12 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_name
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<

        ivar trim -e -i ~{bam} -b ~{primers} -p ~{sample_name}_trim.bam
        samtools sort ~{sample_name}_trim.bam -o ~{sample_name}_trim.sort.bam
        samtools index ~{sample_name}_trim.sort.bam

    >>>

    output {

        File trim_bam = "${sample_name}_trim.bam"
        File trimsort_bam = "${sample_name}_trim.sort.bam"
        File trimsort_bamindex = "${sample_name}_trim.sort.bam.bai"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task ivar_var {

    input {

        String sample_name
        File ref
        File gff
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<

        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{ref} ~{bam} | \
        ivar variants -p ~{sample_name}_variants -q 20 -t 0.6 -m 10 -r ~{ref} -g ~{gff}

    >>>

    output {

        File var_out = "${sample_name}_variants.tsv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task ivar_consensus {

    input {

        String sample_name
        File ref
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<

        # grab ivar and samtools versions
        ivar version | awk '/version/ {print $3}' | tee VERSION_ivar
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_samtools

        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{ref} ~{bam} | \
        ivar consensus -p ~{sample_name}_consensus -q 20 -t 0.6 -m 10

    >>>

    output {

        File consensus_out = "${sample_name}_consensus.fa"

        VersionInfo ivar_version_info = object {
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_ivar")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string ("VERSION_samtools")
        }
    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task bam_stats {

    input {

        String sample_name
        File bam
        File bai
        File s_gene_amplicons
        File primer_bed
    }

    String docker = "staphb/samtools:1.16"

    command <<<

        # grab version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION

        samtools flagstat ~{bam} > ~{sample_name}_flagstat.txt
        samtools stats ~{bam} > ~{sample_name}_stats.txt
        samtools coverage -m -o ~{sample_name}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_name}_coverage.txt ~{bam}
        samtools depth -a -o ~{sample_name}_depth.txt ~{bam}


        # Calculate depth of coverage over entire S gene
        echo "Calculating overall S gene depth"
        samtools coverage --region MN908947.3:21,563-25,384 \
            -o ~{sample_name}_S_gene_coverage.txt ~{bam}

        # Calculate depth of coverage over S gene amplicon regions (excludes overlapping regions with adjacent amplicons)
        echo "calculating depths for ~{s_gene_amplicons}"
        {
            s_gene_depths="~{sample_name}_S_gene_depths.tsv"

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

        File flagstat_out  = "${sample_name}_flagstat.txt"
        File stats_out  = "${sample_name}_stats.txt"
        File covhist_out  = "${sample_name}_coverage_hist.txt"
        File cov_out  = "${sample_name}_coverage.txt"
        File depth_out = "${sample_name}_depth.txt"
        File cov_s_gene_out = "${sample_name}_S_gene_coverage.txt"
        File cov_s_gene_amplicons_out = "${sample_name}_S_gene_depths.tsv"

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION")
        }

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
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

    command <<<
        python ~{calc_percent_coverage_py} \
          --sample_name ~{sample_name} \
          --fasta_file ~{fasta}
      >>>

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
