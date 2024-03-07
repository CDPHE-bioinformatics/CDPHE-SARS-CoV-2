version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture
import "../tasks/hostile_task.wdl" as hostile_task

workflow SC2_ont_assembly {

    input {
        String    gcs_fastq_dir
        String    sample_name
        String    index_1_id
        String    primer_set
        String    barcode_kit
        String    medaka_model
        Boolean  scrub_reads
        File?    scrub_genome_index
        File    covid_genome
        File    primer_bed
        File    s_gene_primer_bed
        File    s_gene_amplicons
        String  project_name
        String  out_dir

        # python scripts
        File    calc_percent_coverage_py
        File    version_capture_ont_assembly_py
    }

    # secret variables
    String outdirpath = sub(out_dir, "/$", "")

    call ListFastqFiles {
        input:
            gcs_fastq_dir = gcs_fastq_dir
    }

    call Demultiplex {
        input:
            fastq_files = ListFastqFiles.fastq_files,
            index_1_id = index_1_id,
            barcode_kit = barcode_kit
    }
    if (scrub_reads) {
        call concatenate_fastqs {
            input:
                sample_name = sample_name,
                fastq_files = Demultiplex.guppy_demux_fastq
        }
        call hostile_task.hostile as hostile {
            input:
                fastq1 = concatenate_fastqs.concatenated_fastq,
                seq_method = "OXFORD_NANOPORE",
                genome_index = [hostile_genome_index],
        }
    }

    call Read_Filtering {
        input:
            fastq_files = select_first([hostile.fastq1_scrubbed, Demultiplex.guppy_demux_fastq]),
            index_1_id = index_1_id,
            sample_name = sample_name,
            primer_set = primer_set
    }

    call Medaka{
        input:
            filtered_reads = Read_Filtering.guppyplex_fastq,
            sample_name = sample_name,
            index_1_id = index_1_id,
            medaka_model = medaka_model
    }

    Boolean consensus_defined = defined(Medaka.consensus)
    if (consensus_defined) {
        Float consensus_size = size(Medaka.consensus)
        Boolean empty_fasta = if (consensus_size < 30) then true else false

        if (empty_fasta) {
            String exit_reason = "Empty fasta"
            call exit_wdl {
                input:
                    exit_reason = exit_reason
            }
        }
    }

    call Bam_stats {
        input:
            bam = Medaka.trimsort_bam,
            bai = Medaka.trimsort_bai,
            sample_name = sample_name,
            index_1_id = index_1_id,
            s_gene_amplicons = s_gene_amplicons,
            primer_bed = primer_bed

    }

    call Scaffold {
        input:
            sample_name = sample_name,
            index_1_id = index_1_id,
            ref = covid_genome,
            fasta = Medaka.consensus
    }

    call rename_fasta {
        input:
            sample_name = sample_name,
            fasta = Scaffold.scaffold_consensus
    }

    call calc_percent_cvg {
        input:
            sample_name = sample_name,
            fasta = rename_fasta.renamed_consensus,
            calc_percent_coverage_py = calc_percent_coverage_py
    }

    call get_primer_site_variants {
        input:
            variants = Medaka.variants,
            variants_index = Medaka.variants_index,
            sample_name = sample_name,
            s_gene_primer_bed = s_gene_primer_bed
    }

    call version_capture.workflow_version_capture as workflow_version_capture{
        input:
    }

    call create_version_capture_file {
        input:
            project_name = project_name,
            guppy_version = Demultiplex.guppy_version,
            artic_version = Medaka.artic_version,
            medaka_version = Medaka.medaka_version,
            samtools_version = Bam_stats.samtools_version,
            pyScaf_version = Scaffold.pyScaf_version,
            bcftools_version = get_primer_site_variants.bcftools_version,
            analysis_date = workflow_version_capture.analysis_date,
            workflow_version = workflow_version_capture.workflow_version,
            version_capture_ont_assembly_py = version_capture_ont_assembly_py
    }

    call transfer {
        input:
        outdirpath = outdirpath,
        fastq_scrubbed = hostile.fastq1_scrubbed,
        trimsort_bam = Medaka.trimsort_bam,
        trimsort_bai = Medaka.trimsort_bai,
        flagstat_out = Bam_stats.flagstat_out,
        samstats_out = Bam_stats.stats_out,
        covhist_out = Bam_stats.covhist_out,
        cov_out = Bam_stats.cov_out,
        depth_out = Bam_stats.depth_out,
        cov_s_gene_out = Bam_stats.cov_s_gene_out,
        cov_s_gene_amplicons_out = Bam_stats.cov_s_gene_amplicons_out,
        variants = Medaka.variants,
        renamed_consensus = rename_fasta.renamed_consensus,
        primer_site_variants = get_primer_site_variants.primer_site_variants,
        version_capture_ont_assembly = create_version_capture_file.version_capture_ont_assembly
    }

    output {
        File index_1_id_summary = Demultiplex.index_1_id_summary
        Array[File] guppy_demux_fastq = Demultiplex.guppy_demux_fastq
        File fastq_files_scrubbed = hostile.fastq1_scrubbed
        Int human_reads_removed = hostile.human_reads_removed
        Float human_reads_removed_proportion = hostile.human_reads_removed_proportion
        String? hostile_version = hostile.hostile_version
        String? hostile_docker = hostile.hostile_docker
        File filtered_fastq = Read_Filtering.guppyplex_fastq
        File sorted_bam = Medaka.sorted_bam
        File trimsort_bam = Medaka.trimsort_bam
        File trimsort_bai = Medaka.trimsort_bai
        File flagstat_out = Bam_stats.flagstat_out
        File samstats_out = Bam_stats.stats_out
        File covhist_out = Bam_stats.covhist_out
        File cov_out = Bam_stats.cov_out
        File depth_out = Bam_stats.depth_out
        File cov_s_gene_out = Bam_stats.cov_s_gene_out
        File cov_s_gene_amplicons_out = Bam_stats.cov_s_gene_amplicons_out
        File variants = Medaka.variants
        File consensus = Medaka.consensus
        File scaffold_consensus = Scaffold.scaffold_consensus
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        File primer_site_variants = get_primer_site_variants.primer_site_variants

        File version_capture_ont_assembly = create_version_capture_file.version_capture_ont_assembly
        String transfer_date_assembly = transfer.transfer_date_assembly
    }
}

task ListFastqFiles {
    input {
        String gcs_fastq_dir
    }

    String indir = sub(gcs_fastq_dir, "/$", "")

    command <<<
        gsutil ls ~{indir}/* > fastq_files.txt
    >>>

    output {
        Array[File] fastq_files = read_lines("fastq_files.txt")
    }

    runtime {
        cpu:    1
        memory:    "1 GB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
}

task Demultiplex {
    input {
        Array[File] fastq_files
        String index_1_id
        String barcode_kit
    }

    Int disk_size = 3 * ceil(size(fastq_files, "GB"))

    command <<<
        guppy_barcoder -v | awk '/Version/ {print $13}' | tee VERSION
        set -e
        mkdir fastq_files
        ln -s ~{sep=' ' fastq_files} fastq_files
        ls -alF fastq_files
        guppy_barcoder --require_barcodes_both_ends --barcode_kits ~{barcode_kit} --fastq_out -i fastq_files -s demux_fastq
        ls -alF demux_fastq
    >>>

    output {
        Array[File] guppy_demux_fastq = glob("demux_fastq/${index_1_id}/*.fastq")
        File index_1_id_summary = "demux_fastq/barcoding_summary.txt"
        String guppy_version = read_string("VERSION")
    }

    runtime {
        cpu:    8
        memory:    "16 GB"
        disks:    "local-disk 100 SSD"
        preemptible:    0
        maxRetries:    3
        docker:    "genomicpariscentre/guppy:6.4.6"
    }
}

task concatenate_fastqs {
    input {
        String sample_name
        Array[File] fastq_files
    }

    command <<<
        cat ~{sep=" " fastq_files} > ~{sample_name}.fastq
    >>>

    output {
        File concatenated_fastq = "~{sample_name}.fastq"
    }

    runtime {
        cpu: 1
        memory: "6 GB"
        docker: "ubuntu:jammy"
    }
}


task Read_Filtering {
    input {
        Array[File] fastq_files
        String index_1_id
        String sample_name
        String primer_set
    }

    Int max_length = if primer_set == "Midnight" then 1500 else 700

    command <<<
        set -e
        mkdir fastq_files
        ln -s ~{sep=' ' fastq_files} fastq_files
        ls -alF fastq_files

        artic guppyplex --min-length 400 --max-length ~{max_length} --directory fastq_files --output ~{sample_name}_~{index_1_id}.fastq

    >>>

    output {
        File guppyplex_fastq = "${sample_name}_${index_1_id}.fastq"
    }

    runtime {
        docker: "quay.io/staphb/artic-ncov2019:1.3.0"
        memory: "16 GB"
        cpu: 8
        disks: "local-disk 100 SSD"
        preemptible: 0
    }
}

task Medaka {
    input {
        String index_1_id
        String sample_name
        File filtered_reads
        String medaka_model
    }

    command <<<

        artic minion --medaka --medaka-model ~{medaka_model} --normalise 20000 --threads 8 --read-file ~{filtered_reads} nCoV-2019 ~{sample_name}_~{index_1_id}

        artic -v > VERSION_artic
        medaka --version | tee VERSION_medaka

    >>>

    output {
        File consensus = "${sample_name}_${index_1_id}.consensus.fasta"
        File sorted_bam = "${sample_name}_${index_1_id}.trimmed.rg.sorted.bam"
        File trimsort_bam = "${sample_name}_${index_1_id}.primertrimmed.rg.sorted.bam"
        File trimsort_bai = "${sample_name}_${index_1_id}.primertrimmed.rg.sorted.bam.bai"
        File variants = "${sample_name}_${index_1_id}.pass.vcf.gz"
        File variants_index = "${sample_name}_${index_1_id}.pass.vcf.gz.tbi"
        String artic_version = read_string("VERSION_artic")
        String medaka_version = read_string("VERSION_medaka")
    }

    runtime {
        docker: "staphb/artic:1.2.4-1.11.1"
        memory: "16 GB"
        cpu: 8
        disks: "local-disk 100 SSD"
        preemptible: 0
    }
}

task exit_wdl {
    input {
        String exit_reason
    }
    
    command <<<
        echo "~{exit_reason}"
        exit 1
    >>>

    runtime {
        return_codes: 0
        docker: "ubuntu:latest"
    }
}

task Bam_stats {
    input {
        String sample_name
        String index_1_id
        File bam
        File bai
        File s_gene_amplicons
        File primer_bed
    }

    Int disk_size = 3 * ceil(size(bam, "GB"))

    command <<<

        # grab samtools version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION

        samtools flagstat ~{bam} > ~{sample_name}_~{index_1_id}_flagstat.txt

        samtools stats ~{bam} > ~{sample_name}_~{index_1_id}_stats.txt

        samtools coverage -m -o ~{sample_name}_~{index_1_id}_coverage_hist.txt ~{bam}


        samtools coverage -o ~{sample_name}_~{index_1_id}_coverage.txt ~{bam}

        samtools depth -a -o ~{sample_name}_~{index_1_id}_depth.txt ~{bam}


        # Calculate depth of coverage over entire S gene
        echo "Calculating overall S gene depth"
        samtools coverage --region MN908947.3:21,563-25,384 \
            -o ~{sample_name}_~{index_1_id}_S_gene_coverage.txt ~{bam}

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

        File flagstat_out  = "${sample_name}_${index_1_id}_flagstat.txt"
        File stats_out  = "${sample_name}_${index_1_id}_stats.txt"
        File covhist_out  = "${sample_name}_${index_1_id}_coverage_hist.txt"
        File cov_out  = "${sample_name}_${index_1_id}_coverage.txt"
        File depth_out = "${sample_name}_${index_1_id}_depth.txt"
        File cov_s_gene_out = "${sample_name}_${index_1_id}_S_gene_coverage.txt"
        File cov_s_gene_amplicons_out = "${sample_name}_S_gene_depths.tsv"
        String samtools_version = read_string("VERSION")
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.16"
    }
}

task Scaffold {
    input {
        String sample_name
        String index_1_id
        File fasta
        File ref
    }

    Int disk_size = 3 * ceil(size(fasta, "GB"))

    command <<<

        # grab version
        pyScaf.py --version > VERSION 2>&1 # writes version to stderr instead of stdout

        pyScaf.py -f ~{fasta} -o ~{sample_name}_~{index_1_id}_consensus_scaffold.fa -r ~{ref}

    >>>

    output {
        File scaffold_consensus = "${sample_name}_${index_1_id}_consensus_scaffold.fa"
        String pyScaf_version = read_string("VERSION")
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "chrishah/pyscaf-docker"
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

      docker: "mchether/py3-bio:v4"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }

}

task get_primer_site_variants {

    meta {
        description: "Get variants at primer-binding sites on the S gene"
    }
    input {
        File variants
        File variants_index
        String sample_name
        File s_gene_primer_bed
    }

    command <<<

        # grab version
        bcftools --version | awk '/bcftools/ {print $2}' | tee VERSION

        bcftools view --regions-file ~{s_gene_primer_bed} --no-header ~{variants} \
            | tee ~{sample_name}_S_gene_primer_variants.txt

    >>>

    output {
        File primer_site_variants = "${sample_name}_S_gene_primer_variants.txt"
        String bcftools_version = read_string("VERSION")
    }

    runtime {
        docker: "staphb/bcftools:1.16"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task create_version_capture_file {


    input {
        File version_capture_ont_assembly_py
        String project_name
        String guppy_version
        String artic_version
        String medaka_version
        String samtools_version
        String pyScaf_version
        String bcftools_version
        String? hostile_version
        String analysis_date
        String workflow_version
    }

    command <<<

        python ~{version_capture_ont_assembly_py} \
        --project_name "~{project_name}" \
        --guppy_version "~{guppy_version}" \
        --artic_version "~{artic_version}" \
        --medaka_version "~{medaka_version}" \
        --samtools_version "~{samtools_version}" \
        --pyScaf_version "~{pyScaf_version}" \
        --bcftools_version "~{bcftools_version}" \
        --hostile_version = "~{hostile_version}" \
        --analysis_date "~{analysis_date}" \
        --workflow_version "~{workflow_version}"

    >>>

    output {
        File version_capture_ont_assembly = 'version_capture_ont_asembly_~{project_name}_~{workflow_version}.csv'
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
        File? fastq_scrubbed
        File trimsort_bam
        File trimsort_bai
        File flagstat_out
        File samstats_out
        File covhist_out
        File cov_out
        File depth_out
        File cov_s_gene_out
        File cov_s_gene_amplicons_out
        File variants
        File renamed_consensus
        File primer_site_variants
        File version_capture_ont_assembly
    }

    command <<<

        gsutil -m cp ~{fastq_scrubbed} ~{outdirpath}/fastq_scrubbed/
        gsutil -m cp ~{trimsort_bam} ~{outdirpath}/alignments/
        gsutil -m cp ~{trimsort_bai} ~{outdirpath}/alignments/
        gsutil -m cp ~{flagstat_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{samstats_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{cov_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{depth_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{cov_s_gene_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{cov_s_gene_amplicons_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{variants} ~{outdirpath}/variants/
        gsutil -m cp ~{primer_site_variants} ~{outdirpath}/primer_site_variants/
        gsutil -m cp ~{renamed_consensus} ~{outdirpath}/assemblies/
        gsutil -m cp ~{version_capture_ont_assembly} ~{outdirpath}/summary_results/


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
