version 1.0

workflow SC2_transfer_illumina_pe_assembly {

    input {
        File fastqc_raw1_html
        File fastqc_raw1_zip
        File fastqc_raw2_html
        File fastqc_raw2_zip
        File fastqc_clean1_html
        File fastqc_clean1_zip
        File fastqc_clean2_html
        File fastqc_clean2_zip
        File seqyclean_summary
        File filtered_reads_1
        File filtered_reads_2
        File trimsort_bam
        File trimsort_bamindex
        File consensus
        File variants
        File cov_out
        File covhist_out
        File flagstat_out
        File stats_out
        File renamed_consensus
        String out_dir
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
        String out_dir
        File fastqc_raw1_html
        File fastqc_raw1_zip
        File fastqc_raw2_html
        File fastqc_raw2_zip
        File fastqc_clean1_html
        File fastqc_clean1_zip
        File fastqc_clean2_html
        File fastqc_clean2_zip
        File seqyclean_summary
        File filtered_reads_1
        File filtered_reads_2
        File trimsort_bam
        File trimsort_bamindex
        File consensus
        File variants
        File cov_out
        File covhist_out
        File flagstat_out
        File stats_out
        File renamed_consensus
    }

    command <<<
        
        gsutil -m cp ~{fastqc_raw1_html} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_raw1_zip} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_raw2_html} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_raw2_zip} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_clean1_html} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_clean1_zip} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_clean2_html} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_clean2_zip} ~{out_dir}/fastqc/
        gsutil -m cp ~{seqyclean_summary} ~{out_dir}/seqyclean/
        gsutil -m cp ~{filtered_reads_1} ~{out_dir}/seqyclean/
        gsutil -m cp ~{filtered_reads_2} ~{out_dir}/seqyclean/
        gsutil -m cp ~{trimsort_bam} ~{out_dir}/alignments/
        gsutil -m cp ~{trimsort_bamindex} ~{out_dir}/alignments/
        gsutil -m cp ~{consensus} ~{out_dir}/assemblies/
        gsutil -m cp ~{variants} ~{out_dir}/variants/
        gsutil -m cp ~{cov_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{flagstat_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{stats_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{renamed_consensus} ~{out_dir}/assemblies/
        
        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
