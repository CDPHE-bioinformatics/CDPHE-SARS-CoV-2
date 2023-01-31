version 1.0

workflow SC2_transfer_illumina_se_assembly {

    input {
        File fastqc_raw_html
        File fastqc_raw_zip
        File fastqc_clean_html
        File fastqc_clean_zip
        File adapter_stats
        File PhiX_stats
        File filtered_reads
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
        String out_dir
        File fastqc_raw_html
        File fastqc_raw_zip
        File fastqc_clean_html
        File fastqc_clean_zip
        File adapter_stats
        File PhiX_stats
        File filtered_reads
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
        
        gsutil -m cp ~{fastqc_raw_html} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_raw_zip} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_clean_html} ~{out_dir}/fastqc/
        gsutil -m cp ~{fastqc_clean_zip} ~{out_dir}/fastqc/
        gsutil -m cp ~{adapter_stats} ~{out_dir}/filtered_reads/
        gsutil -m cp ~{PhiX_stats} ~{out_dir}/filtered_reads/
        gsutil -m cp ~{filtered_reads} ~{out_dir}/filtered_reads/
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
        memory: "1 GB"
        cpu: 4
        disks: "local-disk 500 SSD"
    }
}
