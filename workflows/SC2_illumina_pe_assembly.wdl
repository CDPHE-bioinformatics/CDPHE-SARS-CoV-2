version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture



workflow SC2_illumina_pe_assembly {

    input {
        String    sample_name
        File    fastq_1
        File    fastq_2
        File    primer_bed
        File    adapters_and_contaminants
        File    covid_genome
        File    covid_gff
        String  project_name
        String out_dir

        # python scripts
        File    calc_percent_coverage_py
        File    s_gene_amplicons
        File    version_capture_illumina_pe_assembly_py
    }

    # secrete variables
    String outdirpath = sub(out_dir, "/$", "")

    call seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_name = sample_name,
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
            s_gene_amplicons = s_gene_amplicons
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

    call create_version_capture_file {
        input:
            version_capture_illumina_pe_assembly_py = version_capture_illumina_pe_assembly_py,
            project_name = project_name,
            seqyclean_version = seqyclean.seqyclean_version,
            fastqc_version = fastqc_raw.fastqc_version,
            bwa_version = align_reads.bwa_version,
            samtools_version_broadinstitute = align_reads.samtools_version_broadinstitute,
            ivar_version = ivar_consensus.ivar_version,
            samtools_version_andersenlabapps = ivar_consensus.samtools_version_andersenlabapps,
            samtools_version_staphb = bam_stats.samtools_version_staphb,
            analysis_date = workflow_version_capture.analysis_date,
            workflow_version = workflow_version_capture.workflow_version

            
    }
    call transfer {
        input:
            outdirpath = outdirpath,
            seqyclean_summary = seqyclean.seqyclean_summary,
            fastqc_raw1_html = fastqc_raw.fastqc1_html,
            fastqc_raw1_zip = fastqc_raw.fastqc1_zip,
            fastqc_raw2_html = fastqc_raw.fastqc2_html,
            fastqc_raw2_zip = fastqc_raw.fastqc2_zip,
            fastqc_clean1_html = fastqc_cleaned.fastqc1_html,
            fastqc_clean1_zip = fastqc_cleaned.fastqc1_zip,
            fastqc_clean2_html = fastqc_cleaned.fastqc2_html,
            fastqc_clean2_zip = fastqc_cleaned.fastqc2_zip,
            trimsort_bam = ivar_trim.trimsort_bam,
            trimsort_bamindex = ivar_trim.trimsort_bamindex,
            variants = ivar_var.var_out,
            flagstat_out = bam_stats.flagstat_out,
            stats_out = bam_stats.stats_out,
            covhist_out = bam_stats.covhist_out,
            cov_out = bam_stats.cov_out,
            cov_s_gene_out = bam_stats.cov_s_gene_out,
            cov_s_gene_amplicons_out = bam_stats.cov_s_gene_amplicons_out,
            renamed_consensus = rename_fasta.renamed_consensus,
            version_capture_illumina_pe_assembly = create_version_capture_file.version_capture_illumina_pe_assembly
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

        File version_capture_illumina_pe_assembly = create_version_capture_file.version_capture_illumina_pe_assembly
        String transfer_date_assembly = transfer.transfer_date_assembly

    }
}

task seqyclean {
    input {
        File contam
        String sample_name
        File fastq_1
        File fastq_2
    }

    command <<<

        seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_1} -2 ~{fastq_2} -c ~{contam} -o ~{sample_name}_clean

        # grab seqyclean version 
        seqyclean -h | awk '/Version/ {print $2}' | tee VERSION
    >>>

    output {

        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
        String seqyclean_version = read_string("VERSION")

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
        String fastqc_version = read_string("VERSION")

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
        String sample_name
    }

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
        String bwa_version = read_string("VERSION_bwa")
        String samtools_version_broadinstitute = read_string("VERSION_samtools")

    }

    runtime {
        cpu:    2
        memory:    "12 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "broadinstitute/viral-core:2.0.21"
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_name
    }

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
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_consensus {

    input {

        String sample_name
        File ref
        File bam
    }

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
        String ivar_version = read_string("VERSION_ivar")
        String samtools_version_andersenlabapps = read_string("VERSION_samtools")

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

        String sample_name
        File bam
        File bai
        File s_gene_amplicons
    }

    command <<<

        # grab version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION

        samtools flagstat ~{bam} > ~{sample_name}_flagstat.txt
        samtools stats ~{bam} > ~{sample_name}_stats.txt
        samtools coverage -m -o ~{sample_name}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_name}_coverage.txt ~{bam}

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
        File cov_s_gene_out = "${sample_name}_S_gene_coverage.txt"
        File cov_s_gene_amplicons_out = "${sample_name}_S_gene_depths.tsv"
        String samtools_version_staphb = read_string("VERSION")

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

task create_version_capture_file {
    meta {
        description: "generate version capture file for workflow"
    }

    input {
        File version_capture_illumina_pe_assembly_py
        String project_name
        String seqyclean_version
        String fastqc_version
        String bwa_version
        String samtools_version_broadinstitute
        String ivar_version
        String samtools_version_andersenlabapps
        String samtools_version_staphb
        String analysis_date
        String workflow_version

    }

    command <<<

        python ~{version_capture_illumina_pe_assembly_py} \
        --project_name "~{project_name}" \
        --seqyclean_version "~{seqyclean_version}" \
        --fastqc_version "~{fastqc_version}" \
        --bwa_version "~{bwa_version}" \
        --samtools_version_broadinstitute "~{samtools_version_broadinstitute}" \
        --ivar_version "~{ivar_version}" \
        --samtools_version_andersenlabapps "~{samtools_version_andersenlabapps}" \
        --samtools_version_staphb "~{samtools_version_staphb}" \
        --analysis_date "~{analysis_date}" \
        --workflow_version "~{workflow_version}" \

    >>>

    output {
        File version_capture_illumina_pe_assembly = "version_capture_illumina_pe_asembly_~{project_name}_v~{workflow_version}.csv"
    }

    runtime {

      docker: "mchether/py3-bio:v4"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }
}

task transfer {
    input {
        String outdirpath
        File seqyclean_summary 
        File fastqc_raw1_html
        File fastqc_raw1_zip
        File fastqc_raw2_html
        File fastqc_raw2_zip
        File fastqc_clean1_html
        File fastqc_clean1_zip
        File fastqc_clean2_html
        File fastqc_clean2_zip
        File trimsort_bam
        File trimsort_bamindex
        File variants
        File flagstat_out
        File stats_out
        File covhist_out
        File cov_out
        File cov_s_gene_out
        File cov_s_gene_amplicons_out
        File renamed_consensus
        File version_capture_illumina_pe_assembly       

    }

    command <<<

        gsutil -m cp ~{seqyclean_summary} ~{outdirpath}/seqyclean/
        gsutil -m cp ~{fastqc_raw1_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_raw1_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_raw2_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_raw2_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_clean1_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_clean1_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_clean2_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{fastqc_clean2_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{trimsort_bam} ~{outdirpath}/alignments/
        gsutil -m cp ~{trimsort_bamindex} ~{outdirpath}/alignments/
        gsutil -m cp ~{variants} ~{outdirpath}/variants/
        gsutil -m cp ~{cov_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{cov_s_gene_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{cov_s_gene_amplicons_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{flagstat_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{stats_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{renamed_consensus} ~{outdirpath}/assemblies/
        gsutil -m cp ~{version_capture_illumina_pe_assembly} ~{outdirpath}/summary_results/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE

    >>>


    output {
        String transfer_date_assembly = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}