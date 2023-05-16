version 1.0

workflow SC2_ont_assembly {

    input {
        String    gcs_fastq_dir
        String    sample_name
        String    index_1_id
        String    primer_set
        File    covid_genome
        File    primer_bed
        File    s_gene_primer_bed
        File    s_gene_amplicons
        String project_name

        # python scripts
        File    calc_percent_coverage_py
        File    concat_assembly_software_ont_py
    }
    
    call ListFastqFiles {
        input:
            gcs_fastq_dir = gcs_fastq_dir
    }
    call Demultiplex {
        input:
            fastq_files = ListFastqFiles.fastq_files,
            index_1_id = index_1_id
    }
    call Read_Filtering {
        input:
            fastq_files = Demultiplex.guppy_demux_fastq,
            index_1_id = index_1_id,
            sample_name = sample_name,
            primer_set = primer_set
    }
    call Medaka {
        input:
            filtered_reads = Read_Filtering.guppyplex_fastq,
            sample_name = sample_name,
            index_1_id = index_1_id,
            primer_bed = primer_bed
    }
    call Bam_stats {
        input:
            bam = Medaka.trimsort_bam,
            bai = Medaka.trimsort_bai,
            sample_name = sample_name,
            index_1_id = index_1_id,
            s_gene_amplicons = s_gene_amplicons
            
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

    call create_software_assembly_file {
        input:
            concat_assembly_software_ont_py = concat_assemby_software_py,
            guppy_version = Demultiplex.guppy_version,
            medaka_version = Medaka.assembler_version,
            project_name = project_name
    }

    output {
        File index_1_id_summary = Demultiplex.index_1_id_summary
        Array[File] guppy_demux_fastq = Demultiplex.guppy_demux_fastq
        File filtered_fastq = Read_Filtering.guppyplex_fastq
        File sorted_bam = Medaka.sorted_bam
        File trimsort_bam = Medaka.trimsort_bam
        File trimsort_bai = Medaka.trimsort_bai
        File flagstat_out = Bam_stats.flagstat_out
        File samstats_out = Bam_stats.stats_out
        File covhist_out = Bam_stats.covhist_out
        File cov_out = Bam_stats.cov_out
        File cov_s_gene_out = Bam_stats.cov_s_gene_out
        File cov_s_gene_amplicons_out = Bam_stats.cov_s_gene_amplicons_out
        File variants = Medaka.variants
        File consensus = Medaka.consensus
        File scaffold_consensus = Scaffold.scaffold_consensus
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        File primer_site_variants = get_primer_site_variants.primer_site_variants
        String assembler_version = Medaka.assembler_version
        File assembly_software_file = create_software_assembly_file.assembly_software_file
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
    }

    Int disk_size = 3 * ceil(size(fastq_files, "GB"))

    command <<<
        guppy_barcoder -v | awk '/Version/ {print $13}' | tee VERSION
        set -e
        mkdir fastq_files
        ln -s ~{sep=' ' fastq_files} fastq_files
        ls -alF fastq_files
        guppy_barcoder --require_barcodes_both_ends --barcode_kits "EXP-NBD196" --fastq_out -i fastq_files -s demux_fastq
        ls -alF demux_fastq
    >>>

    output {
        Array[File] guppy_demux_fastq = glob("demux_fastq/${index_1_id}/*.fastq")
        File index_1_id_summary = "demux_fastq/barcoding_summary.txt"
        String guppy_version = read_string('VERSION')
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
        File primer_bed
    }

    command <<<
    
        mkdir -p ./primer-schemes/nCoV-2019/Vuser
        cp /primer-schemes/nCoV-2019/V3/nCoV-2019.reference.fasta ./primer-schemes/nCoV-2019/Vuser/nCoV-2019.reference.fasta
        cp ~{primer_bed} ./primer-schemes/nCoV-2019/Vuser/nCoV-2019.scheme.bed
        
        artic -v > VERSION
        artic minion --medaka --medaka-model r941_min_high_g360 --normalise 20000 --threads 8 --scheme-directory ./primer-schemes --read-file ~{filtered_reads} nCoV-2019/Vuser ~{sample_name}_~{index_1_id}

    >>>

    output {
        File consensus = "${sample_name}_${index_1_id}.consensus.fasta"
        File sorted_bam = "${sample_name}_${index_1_id}.trimmed.rg.sorted.bam"
        File trimsort_bam = "${sample_name}_${index_1_id}.primertrimmed.rg.sorted.bam"
        File trimsort_bai = "${sample_name}_${index_1_id}.primertrimmed.rg.sorted.bam.bai"
        File variants = "${sample_name}_${index_1_id}.pass.vcf.gz"
        File variants_index = "${sample_name}_${index_1_id}.pass.vcf.gz.tbi"
        String assembler_version = read_string("VERSION")
    }

    runtime {
        docker: "quay.io/staphb/artic-ncov2019:1.3.0"
        memory: "16 GB"
        cpu: 8
        disks: "local-disk 100 SSD"
        preemptible: 0
    }
}

task Bam_stats {
    input {
        String sample_name
        String index_1_id
        File bam
        File bai
        File s_gene_amplicons
    }

    Int disk_size = 3 * ceil(size(bam, "GB"))

    command <<<

        samtools flagstat ~{bam} > ~{sample_name}_~{index_1_id}_flagstat.txt

        samtools stats ~{bam} > ~{sample_name}_~{index_1_id}_stats.txt

        samtools coverage -m -o ~{sample_name}_~{index_1_id}_coverage_hist.txt ~{bam}


        samtools coverage -o ~{sample_name}_~{index_1_id}_coverage.txt ~{bam}


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
        File cov_s_gene_out = "${sample_name}_${index_1_id}_S_gene_coverage.txt"
        File cov_s_gene_amplicons_out = "${sample_name}_S_gene_depths.tsv"
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

    command {

        pyScaf.py -f ${fasta} -o ${sample_name}_${index_1_id}_consensus_scaffold.fa -r ${ref}

    }

    output {
        File scaffold_consensus = "${sample_name}_${index_1_id}_consensus_scaffold.fa"
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

    command {

        sed 's/>.*/>CO-CDPHE-~{sample_name}/' ~{fasta} > ~{sample_name}_consensus_renamed.fa

    }

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

        bcftools view --regions-file ~{s_gene_primer_bed} --no-header ~{variants} \
            | tee ~{sample_name}_S_gene_primer_variants.txt

    >>>

    output {
        File primer_site_variants = "${sample_name}_S_gene_primer_variants.txt"
    }

    runtime {
        docker: "staphb/bcftools:1.16"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task create_software_assembly_file {
    meta {
        description: "pull assembly software into a sinlge tsv file"
    }

    input {
        File concat_assembly_software_ont_py
        String guppy_version
        String medaka_version
        String project_name
    }

    command <<<

        python ~{concat_assembly_software_ont_py} \
        --project_name "~{project_name}" \
        --bwa_version "~{guppy_version}" \
        --ivar_version "~{medaka_version}"

    >>>

    output {
        File assemlby_software_file = '~{project_name}_assembly_software.tsv'
    }
}
