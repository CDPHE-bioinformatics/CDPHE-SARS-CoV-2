version 1.0

workflow SC2_transfer_ont_assembly {

    input {
        File filtered_fastq
        File trimsort_bam
        File flagstat_out
        File samstats_out
        File covhist_out
        File cov_out
        File variants
        File scaffold_consensus
        File renamed_consensus
        File cov_s_gene_out
        File cov_s_gene_amplicons_out
        File primer_site_variants
        String out_dir
    }

    call transfer_outputs {
        input:
            filtered_fastq = filtered_fastq,
            trimsort_bam = trimsort_bam,
            flagstat_out = flagstat_out,
            samstats_out = samstats_out,
            covhist_out = covhist_out,
            cov_out = cov_out,
            cov_s_gene_out = cov_s_gene_out,
            cov_s_gene_amplicons_out = cov_s_gene_amplicons_out,
            variants = variants,
            scaffold_consensus = scaffold_consensus,
            renamed_consensus = renamed_consensus,
            primer_site_variants = primer_site_variants,
            out_dir = out_dir
    }
    
    output {
        String transfer_date = transfer_outputs.transfer_date
    }
}    

task transfer_outputs {
    input {
        File filtered_fastq
        File trimsort_bam
        File flagstat_out
        File samstats_out
        File covhist_out
        File cov_out
        File cov_s_gene_out
        File cov_s_gene_amplicons_out
        File variants
        File scaffold_consensus
        File renamed_consensus
        File primer_site_variants
        String out_dir
    }

    command <<<
        
        gsutil -m cp ~{filtered_fastq} ~{out_dir}/filtered_fastq/
        gsutil -m cp ~{trimsort_bam} ~{out_dir}/alignments/
        gsutil -m cp ~{flagstat_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{samstats_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{cov_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{cov_s_gene_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{cov_s_gene_amplicons_out} ~{out_dir}/bam_stats/
        gsutil -m cp ~{scaffold_consensus} ~{out_dir}/assemblies/
        gsutil -m cp ~{variants} ~{out_dir}/variants/
        gsutil -m cp ~{primer_site_variants} ~{out_dir}/primer_site_variants/
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
        disks: "local-disk 100 SSD"
    }
}
