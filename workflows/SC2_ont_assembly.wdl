version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture
import "../tasks/hostile_task.wdl" as hostile_task
import "../tasks/transfer_task.wdl" as transfer_task

workflow SC2_ont_assembly {

    input {
        String    gcs_fastq_dir
        String    sample_name
        String    index_1_id
        String    primer_set
        String    barcode_kit
        String?    model
        Boolean  scrub_reads
        File?    scrub_genome_index
        File    covid_genome
        File    primer_bed
        File    s_gene_primer_bed
        File    s_gene_amplicons
        String  project_name
        String  out_dir
        Boolean overwrite

        # python scripts
        File    calc_percent_coverage_py
        File    version_capture_py
    }

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
                genome_index = select_all([scrub_genome_index]),
        }
        call Read_Filtering as Scrubbed_Read_Filtering {
            input:
                fastq_files = [hostile.fastq1_scrubbed],
                index_1_id = index_1_id,
                sample_name = sample_name,
                primer_set = primer_set
        }
    }

    if (!scrub_reads) { # no 'else' in WDL yet
        call Read_Filtering {
            input:
                fastq_files = Demultiplex.guppy_demux_fastq,
                index_1_id = index_1_id,
                sample_name = sample_name,
                primer_set = primer_set
        }
    }

    call call_consensus_artic{
        input:
            filtered_reads = select_first([Read_Filtering.guppyplex_fastq, Scrubbed_Read_Filtering.guppyplex_fastq]),
            fastq_file = Demultiplex.guppy_demux_fastq[0],
            sample_name = sample_name,
            index_1_id = index_1_id,
            model = model,
            ref = covid_genome,
            bed = primer_bed,
    }

    Boolean consensus_defined = defined(call_consensus_artic.consensus)
    if (consensus_defined) {
        Float consensus_size = size(call_consensus_artic.consensus)
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
            bam = call_consensus_artic.trimsort_bam,
            bai = call_consensus_artic.trimsort_bai,
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
            fasta = call_consensus_artic.consensus
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
            variants = call_consensus_artic.variants,
            variants_index = call_consensus_artic.variants_index,
            sample_name = sample_name,
            s_gene_primer_bed = s_gene_primer_bed
    }

    call version_capture.workflow_version_capture as workflow_version_capture{
        input:
    }

    Array[VersionInfo] version_array = [
        Demultiplex.guppy_version_info,
        call_consensus_artic.artic_version_info,
        call_consensus_artic.clair3_version_info,
        call_consensus_artic.model_version_info,
        Bam_stats.samtools_version_info,
        Scaffold.pyscaf_version_info,
        get_primer_site_variants.bcftools_version_info,
    ]
    if (scrub_reads) {
        Array[VersionInfo] version_array_with_hostile = flatten([version_array, select_all([hostile.hostile_version_info])])
    }

    call  version_capture.task_version_capture as task_version_capture {
        input:
            version_array = select_first([version_array_with_hostile, version_array]),
            workflow_name = "SC2_ont_assembly",
            workflow_version_path = workflow_version_capture.workflow_version_path,
            project_name = project_name,
            sample_name = sample_name,
            analysis_date = workflow_version_capture.analysis_date,
            version_capture_py = version_capture_py
    }

    SubdirsToFiles subdirs_to_files = object { subdirs_to_files: [
        ("fastq_scrubbed", [hostile.fastq1_scrubbed]),
        ("alignments",
            [call_consensus_artic.trimsort_bam,
             call_consensus_artic.trimsort_bai]),
        ("bam_stats",
            [Bam_stats.flagstat_out, Bam_stats.stats_out, Bam_stats.covhist_out,
             Bam_stats.cov_out, Bam_stats.depth_out, Bam_stats.cov_s_gene_out,
             Bam_stats.cov_s_gene_amplicons_out]),
        ("variants", [call_consensus_artic.variants]),
        ("primer_site_variants", [get_primer_site_variants.primer_site_variants]),
        ("assemblies", [rename_fasta.renamed_consensus]),
        ("sample_version_capture", [task_version_capture.version_capture_file])
    ]}

    call transfer_task.transfer {
        input:
            out_dir = out_dir,
            overwrite = overwrite,
            subdirs_to_files = subdirs_to_files
    }

    output {
        File index_1_id_summary = Demultiplex.index_1_id_summary
        Array[File] guppy_demux_fastq = Demultiplex.guppy_demux_fastq
        File? fastq_files_scrubbed = hostile.fastq1_scrubbed
        Int? human_reads_removed = hostile.human_reads_removed
        Float? human_reads_removed_proportion = hostile.human_reads_removed_proportion
        File filtered_fastq = select_first([Read_Filtering.guppyplex_fastq, Scrubbed_Read_Filtering.guppyplex_fastq])
        File sorted_bam = call_consensus_artic.sorted_bam
        File trimsort_bam = call_consensus_artic.trimsort_bam
        File trimsort_bai = call_consensus_artic.trimsort_bai
        File flagstat_out = Bam_stats.flagstat_out
        File samstats_out = Bam_stats.stats_out
        File covhist_out = Bam_stats.covhist_out
        File cov_out = Bam_stats.cov_out
        File depth_out = Bam_stats.depth_out
        File cov_s_gene_out = Bam_stats.cov_s_gene_out
        File cov_s_gene_amplicons_out = Bam_stats.cov_s_gene_amplicons_out
        File variants = call_consensus_artic.variants
        File consensus = call_consensus_artic.consensus
        File scaffold_consensus = Scaffold.scaffold_consensus
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        File primer_site_variants = get_primer_site_variants.primer_site_variants
        File version_capture_ont_assembly = task_version_capture.version_capture_file
        String transfer_date_assembly = transfer.transfer_date
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
        Array[String] fastq_files = read_lines("fastq_files.txt")
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

    String docker = "genomicpariscentre/guppy:6.4.6"
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

        VersionInfo guppy_version_info = object {
            software: "guppy",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    8
        memory:    "16 GB"
        disks:    "local-disk 100 SSD"
        preemptible:    0
        maxRetries:    3
        docker:    docker
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

task call_consensus_artic {
    input {
        String index_1_id
        String sample_name
        String? model
        File filtered_reads
        File fastq_file  # need unzipped FASTQ with original header for model detection
        File ref
        File bed
    }

    String docker = "sambaird/artic-fieldbioinformatics:1.5.3"

    command <<<

        PATH=$PATH:/opt/conda/bin
        export CONDA_PREFIX=/opt/conda

        # Auto-detect model from FASTQ if not provided
        if [[ -z "~{model}" ]]; then
            model=$(
                python -c \
                "import artic.utils; \
                model_dict = artic.utils.choose_model(\"~{fastq_file}\"); \
                print(model_dict['name'])"
            )
        else
            model="~{model}"
        fi

        artic minion \
            --model "${model}" \
            --normalise 20000 \
            --read-file "~{filtered_reads}" \
            --ref "~{ref}" \
            --bed "~{bed}" \
            "~{sample_name}_~{index_1_id}"

        artic -v > VERSION_artic
        run_clair3.sh --version | tee VERSION_clair3
        echo "${model}" > VERSION_model

    >>>

    output {
        File consensus = "${sample_name}_${index_1_id}.consensus.fasta"
        File sorted_bam = "${sample_name}_${index_1_id}.trimmed.rg.sorted.bam"
        File trimsort_bam = "${sample_name}_${index_1_id}.primertrimmed.rg.sorted.bam"
        File trimsort_bai = "${sample_name}_${index_1_id}.primertrimmed.rg.sorted.bam.bai"
        File variants = "${sample_name}_${index_1_id}.pass.vcf.gz"
        File variants_index = "${sample_name}_${index_1_id}.pass.vcf.gz.tbi"
        
        VersionInfo artic_version_info = object {
            software: "artic",
            docker: docker,
            version: read_string("VERSION_artic")
        }

        VersionInfo clair3_version_info = object {
            software: "clair3",
            docker: docker,
            version: read_string("VERSION_clair3")
        }

        VersionInfo model_version_info = object {
            software: "model",
            docker: docker,
            version: read_string("VERSION_model")
        }
    }

    runtime {
        docker: docker
        memory: "16 GB"
        cpu: 8
        disks: "local-disk 100 SSD"
        bootDiskSizeGb: 15  # Since Terra allocating 0 GB for docker container for some reason
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

    String docker = "staphb/samtools:1.16"
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
        
        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task Scaffold {
    input {
        String sample_name
        String index_1_id
        File fasta
        File ref
    }

    String docker = "chrishah/pyscaf-docker"
    Int disk_size = 3 * ceil(size(fasta, "GB"))

    command <<<

        # grab version
        pyScaf.py --version > VERSION 2>&1 # writes version to stderr instead of stdout

        pyScaf.py -f ~{fasta} -o ~{sample_name}_~{index_1_id}_consensus_scaffold.fa -r ~{ref}

    >>>

    output {
        File scaffold_consensus = "${sample_name}_${index_1_id}_consensus_scaffold.fa"

        VersionInfo pyscaf_version_info = object {
            software: "pyScaf",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        cpu:    4
        memory:    "8 GiB"
        disks:    "local-disk " + disk_size + " HDD"
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

    String docker = "staphb/bcftools:1.16"

    command <<<

        # grab version
        bcftools --version | awk '/bcftools/ {print $2}' | tee VERSION

        bcftools view --regions-file ~{s_gene_primer_bed} --no-header ~{variants} \
            | tee ~{sample_name}_S_gene_primer_variants.txt

    >>>

    output {
        File primer_site_variants = "${sample_name}_S_gene_primer_variants.txt"
        String bcftools_version = read_string("VERSION")

        VersionInfo bcftools_version_info = object {
            software: "bcftools",
            docker: docker,
            version: read_string("VERSION")
        }
    }

    runtime {
        docker: docker
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

