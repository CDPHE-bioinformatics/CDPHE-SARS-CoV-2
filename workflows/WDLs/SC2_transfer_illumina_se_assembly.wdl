version 1.0

workflow SC2_transfer_illumina_se_assembly {

    input {
        Array[File] fastqc_raw_html
        Array[File] fastqc_raw_zip
        Array[File] fastqc_clean_html
        Array[File] fastqc_clean_zip
        Array[File] adapter_stats
        Array[File] PhiX_stats
        Array[File] filtered_reads
        Array[File] trimsort_bam
        Array[File] trimsort_bamindex
        Array[File] consensus
        Array[File] variants
        Array[File] cov_out
        Array[File] covhist_out
        Array[File] flagstat_out
        Array[File] stats_out
        Array[File] renamed_consensus
        Array[String] out_dir
    }

    call transfer_outputs {
        input:
            fastqc_raw_html = fastqc_raw_html,
            fastqc_raw_zip = fastqc_raw_zip,
            fastqc_clean_html = fastqc_clean_html,
            fastqc_clean_zip = fastqc_clean_zip,
            adapter_stats = adapter_stats,
            PhiX_stats = PhiX_stats,
            filtered_reads = filtered_reads,
            trimsort_bam = trimsort_bam,
            trimsort_bamindex = trimsort_bamindex,
            consensus = consensus,
            variants = variants,
            cov_out = cov_out,
            covhist_out = covhist_out,
            flagstat_out = flagstat_out,
            stats_out = stats_out,
            renamed_consensus = renamed_consensus,
            out_dir = out_dir
    }
    
    output {
        String transfer_date = transfer_outputs.transfer_date
    }
}    

task transfer_outputs {
    input {
        Array[String] out_dir
        Array[File] fastqc_raw_html
        Array[File] fastqc_raw_zip
        Array[File] fastqc_clean_html
        Array[File] fastqc_clean_zip
        Array[File] adapter_stats
        Array[File] PhiX_stats
        Array[File] filtered_reads
        Array[File] trimsort_bam
        Array[File] trimsort_bamindex
        Array[File] consensus
        Array[File] variants
        Array[File] cov_out
        Array[File] covhist_out
        Array[File] flagstat_out
        Array[File] stats_out
        Array[File] renamed_consensus
    }
    
    String outdir = '${out_dir[0]}'
    String outdirpath = sub(outdir, "/$", "")

    command <<<
        
        gsutil -m cp ~{sep=' ' fastqc_raw_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_raw_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_clean_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_clean_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' adapter_stats} ~{outdirpath}/filtered_reads/
        gsutil -m cp ~{sep=' ' PhiX_stats} ~{outdirpath}/filtered_reads/
        gsutil -m cp ~{sep=' ' filtered_reads} ~{outdirpath}/filtered_reads/
        gsutil -m cp ~{sep=' ' trimsort_bam} ~{outdirpath}/alignments/
        gsutil -m cp ~{sep=' ' trimsort_bamindex} ~{outdirpath}/alignments/
        gsutil -m cp ~{sep=' ' consensus} ~{outdirpath}/assemblies/
        gsutil -m cp ~{sep=' ' variants} ~{outdirpath}/variants/
        gsutil -m cp ~{sep=' ' cov_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' covhist_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' flagstat_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' stats_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' renamed_consensus} ~{outdirpath}/assemblies/
        
        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 500 SSD"
    }
}
