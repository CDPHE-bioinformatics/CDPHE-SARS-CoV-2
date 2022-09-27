version 1.0

workflow SC2_transfer_illumina_pe_assembly {

    input {
        Array[File] fastqc_raw1_html
        Array[File] fastqc_raw1_zip
        Array[File] fastqc_raw2_html
        Array[File] fastqc_raw2_zip
        Array[File] fastqc_clean1_html
        Array[File] fastqc_clean1_zip
        Array[File] fastqc_clean2_html
        Array[File] fastqc_clean2_zip
        Array[File] seqyclean_summary
        Array[File] filtered_reads_1
        Array[File] filtered_reads_2
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
            fastqc_raw1_html = fastqc_raw1_html,
            fastqc_raw1_zip = fastqc_raw1_zip,
            fastqc_raw2_html = fastqc_raw2_html,
            fastqc_raw2_zip = fastqc_raw2_zip,
            fastqc_clean1_html = fastqc_clean1_html,
            fastqc_clean1_zip = fastqc_clean1_zip,
            fastqc_clean2_html = fastqc_clean2_html,
            fastqc_clean2_zip = fastqc_clean2_zip,
            seqyclean_summary = seqyclean_summary,
            filtered_reads_1 = filtered_reads_1,
            filtered_reads_2 = filtered_reads_2,
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
        Array[File] fastqc_raw1_html
        Array[File] fastqc_raw1_zip
        Array[File] fastqc_raw2_html
        Array[File] fastqc_raw2_zip
        Array[File] fastqc_clean1_html
        Array[File] fastqc_clean1_zip
        Array[File] fastqc_clean2_html
        Array[File] fastqc_clean2_zip
        Array[File] seqyclean_summary
        Array[File] filtered_reads_1
        Array[File] filtered_reads_2
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
        
        gsutil -m cp ~{sep=' ' fastqc_raw1_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_raw1_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_raw2_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_raw2_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_clean1_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_clean1_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_clean2_html} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' fastqc_clean2_zip} ~{outdirpath}/fastqc/
        gsutil -m cp ~{sep=' ' seqyclean_summary} ~{outdirpath}/seqyclean/
        gsutil -m cp ~{sep=' ' filtered_reads_1} ~{outdirpath}/seqyclean/
        gsutil -m cp ~{sep=' ' filtered_reads_2} ~{outdirpath}/seqyclean/
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
        disks: "local-disk 100 SSD"
    }
}
