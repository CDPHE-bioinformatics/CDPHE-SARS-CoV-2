version 1.0

workflow SC2_transfer_ont_assembly {

    input {
        Array[File] filtered_fastq
        Array[File] trimsort_bam
        Array[File] flagstat_out
        Array[File] samstats_out
        Array[File] covhist_out
        Array[File] cov_out
        Array[File] variants
        Array[File] scaffold_consensus
        Array[File] renamed_consensus
        Array[File] cov_s_gene_out
        Array[File] cov_s_gene_amplicons_out
        Array[File] primer_site_variants
        Array[String] out_dir
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
        Array[File] filtered_fastq
        Array[File] trimsort_bam
        Array[File] flagstat_out
        Array[File] samstats_out
        Array[File] covhist_out
        Array[File] cov_out
        Array[File] cov_s_gene_out
        Array[File] cov_s_gene_amplicons_out
        Array[File] variants
        Array[File] scaffold_consensus
        Array[File] renamed_consensus
        Array[File] primer_site_variants
        Array[String] out_dir
    }
    
    String outdir = '${out_dir[0]}'
    String outdirpath = sub(outdir, "/$", "")

    command <<<
        
        gsutil -m cp ~{sep=' ' filtered_fastq} ~{outdirpath}/filtered_fastq/
        gsutil -m cp ~{sep=' ' trimsort_bam} ~{outdirpath}/alignments/
        gsutil -m cp ~{sep=' ' flagstat_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' samstats_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' covhist_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' cov_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' cov_s_gene_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' cov_s_gene_amplicons_out} ~{outdirpath}/bam_stats/
        gsutil -m cp ~{sep=' ' scaffold_consensus} ~{outdirpath}/assemblies/
        gsutil -m cp ~{sep=' ' variants} ~{outdirpath}/variants/
        gsutil -m cp ~{sep=' ' primer_site_variants} ~{outdirpath}/primer_site_variants/
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
