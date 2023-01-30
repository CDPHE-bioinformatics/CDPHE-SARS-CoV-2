version 1.0

workflow SC2_wastewater_variant_calling {

    input {
        Array[File] trimsort_bam
        File covid_genome
        File voc_bed
        File voc_annotations
        Array[String] sample_id
        String out_dir
        Array[String] seq_run_array
        File dashboard_formatting_py
    }

    scatter (id_bam in zip(sample_id, trimsort_bam)) {
        call add_RG {
            input:
                sample_id = id_bam.left,
                bam = id_bam.right
        }
        call variant_calling {
            input:
                bam = add_RG.rgbam,
                ref = covid_genome,
                sample_id = id_bam.left
        }
        call freyja_demix {
            input:
                variants = variant_calling.variants,
                depth = variant_calling.depth,
                sample_id = id_bam.left
        }
        call fill_NA {
            input:
                variants = variant_calling.variants,
                sample_id = id_bam.left,
                voc_bed = voc_bed
        }
        call reformat_tsv {
            input:
                tsv = fill_NA.fill_NA_tsv,
                sample_id = id_bam.left
        }
        call summary_prep {
            input:
                tsv = reformat_tsv.reformatted_tsv,
                sample_id = id_bam.left,
                voc_annotations = voc_annotations
        }
        call demix_reformat {
            input:
                tsv = freyja_demix.demix,
                sample_id = id_bam.left
        }
    }

    call freyja_aggregate {
        input:
            demix = freyja_demix.demix
    }

    call format_freyja_aggregate_for_dashboard{
        input:
            dashboard_formatting_py = dashboard_formatting_py,
            freyja_aggregate_file = freyja_aggregate.demix_aggregated,
            seq_run_array = seq_run_array
    }

    call combine_tsv {
        input:
            tsv = summary_prep.sample_voc_tsv_summary,
            tsv_counts = summary_prep.sample_voc_tsv_counts,
            voc_annotations = voc_annotations,
            demix_tsv = demix_reformat.demix_reformatted
    }
    call summary_tsv {
        input:
            tsv = combine_tsv.voc_summary_temp
    }
    call transfer_outputs {
        input:
            variants = variant_calling.variants,
            depth = variant_calling.depth,
            demix = freyja_demix.demix,
            sample_voc_tsv_summary = summary_prep.sample_voc_tsv_summary,
            sample_voc_tsv_counts = summary_prep.sample_voc_tsv_counts,
            voc_counts = combine_tsv.voc_counts,
            voc_summary = summary_tsv.voc_summary,
            demix_aggregated = freyja_aggregate.demix_aggregated,
            demix_summary = combine_tsv.demix_summary,
            wwt_variant_abundances_long_format = format_freyja_aggregate_for_dashboard.wwt_variant_abundances_long_format,
            wwt_variant_abundances_long_format_mean = format_freyja_aggregate_for_dashboard.wwt_variant_abundances_long_format_mean,
            out_dir = out_dir
    }

    output {
        Array[File] addrg_bam = add_RG.rgbam
        Array[File] variants = variant_calling.variants
        Array[File] depth = variant_calling.depth
        Array[File] demix = freyja_demix.demix
        Array[File] fill_NA_tsv = fill_NA.fill_NA_tsv
        Array[File] reformatted_tsv = reformat_tsv.reformatted_tsv
        Array[File] sample_voc_tsv_summary = summary_prep.sample_voc_tsv_summary
        Array[File] sample_voc_tsv_counts = summary_prep.sample_voc_tsv_counts
        Array[File] demix_reformatted = demix_reformat.demix_reformatted
        File demix_aggregated = freyja_aggregate.demix_aggregated
        File voc_summary_temp = combine_tsv.voc_summary_temp
        File voc_counts = combine_tsv.voc_counts
        File voc_summary = summary_tsv.voc_summary
        File demix_summary = combine_tsv.demix_summary
        String transfer_date = transfer_outputs.transfer_date
        File wwt_variant_abundances_long_format = format_freyja_aggregate_for_dashboard.wwt_variant_abundances_long_format
        File wwt_variant_abundances_long_format_mean = format_freyja_aggregate_for_dashboard.wwt_variant_abundances_long_format_mean

    }
}

task add_RG {
    input {
        String sample_id
        File bam
    }

    command <<<

        samtools addreplacerg -r ID:~{sample_id} -r LB:L1 -r SM:~{sample_id} -o ~{sample_id}_addRG.bam ~{bam}

    >>>

    output {
        File rgbam = "${sample_id}_addRG.bam"
    }

    runtime {
        docker: "staphb/samtools:1.10"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}

task variant_calling {
    input {
        String sample_id
        File bam
        File ref

    }

    command <<<

        freyja variants ~{bam} --variants ~{sample_id}_variants.tsv --depths ~{sample_id}_depth.tsv --ref ~{ref}

    >>>

    output {
        File variants = "${sample_id}_variants.tsv"
        File depth = "${sample_id}_depth.tsv"
    }

    runtime {
        docker: "staphb/freyja"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
    }
}

task freyja_demix {
    input {
        String sample_id
        File variants
        File depth

    }

    command <<<

        #get updated lineages for demixing
        mkdir ./freyja_db
        freyja update --outdir ./freyja_db
        
        #creates a temp file with the same name as the intended output file that will get output in case of failure or overwritten in case of sucess
        echo -e "\t~{sample_id}\nsummarized\tLowCov\nlineages\tLowCov\nabundances\tLowCov\nresid\tLowCov\ncoverage\tLowCov" > ~{sample_id}_demixed.tsv
        
        freyja demix --eps 0.01 --covcut 10 --barcodes ./freyja_db/usher_barcodes.csv --meta ./freyja_db/curated_lineages.json --confirmedonly ~{variants} ~{depth} --output ~{sample_id}_demixed.tsv

    >>>

    output {
        File demix = "${sample_id}_demixed.tsv"
    }

    runtime {
        docker: "staphb/freyja"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
        continueOnReturnCode: [0, 1]
    }
}

task freyja_aggregate {
    input {
        Array[File] demix
    }

    command <<<

        mkdir demix_outputs
        mv ~{sep=' ' demix} demix_outputs/
        freyja aggregate demix_outputs/ --output demix_aggregated.tsv

    >>>

    output {
        File demix_aggregated = "demix_aggregated.tsv"
    }

    runtime {
        docker: "staphb/freyja"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
    }
}

task format_freyja_aggregate_for_dashboard {
    input {
        File dashboard_formatting_py
        File freyja_aggregate_file
        Array[String] seq_run_array

    }

    command <<<
    python ~{dashboard_formatting_py} \
        --freyja_aggregate_path ~{freyja_aggregate_file} \
        --seq_run ~{write_lines(seq_run_array)}
    >>>

    output{
        File wwt_variant_abundances_long_format = "wwt_variant_abundances_long_format_${seq_run_array[0]}.csv"
        File wwt_variant_abundances_long_format_mean = "wwt_variant_abundances_long_format_mean_${seq_run_array[0]}.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task fill_NA {
    input {
        File variants
        String sample_id
        File voc_bed
    }

    command <<<

    #input is the freyja variants tsv, we first need to cut and order the columns we want
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $12 "\t" $5 "\t" $8 "\t" $11}' ~{variants} | awk 'NR==1; NR > 1 {print $0 | "sort -n -k 2,2"}' > ~{sample_id}_voc_mutations_temp1.tsv

    #filter tsv with mutations bed file to get the voc associated sites
    grep -f ~{voc_bed} ~{sample_id}_voc_mutations_temp1.tsv > ~{sample_id}_voc_mutations_temp2.tsv

    #fix the header names
    echo -e "CHROM\tPOS\tREF\t~{sample_id}_ALT\t~{sample_id}_DP\t~{sample_id}_RefDP\t~{sample_id}_AltDP\t~{sample_id}_AltFREQ" | cat - ~{sample_id}_voc_mutations_temp2.tsv | sed 's/\t/_/' | sort -t $'\t' -k1,1 > ~{sample_id}_voc_mutations_temp3.tsv

    #generate key file from the voc mutations bed file
    cat ~{voc_bed} | cut -f 1,2 | tr "\t" "_" | sort | uniq > keys.txt

    #get the NA-filled columns we want
    join -t $'\t' -e NA -a 1 -1 1 -2 1 -o "1.1,2.2,2.3,2.4,2.5,2.6,2.7" keys.txt "~{sample_id}_voc_mutations_temp3.tsv" > ~{sample_id}_voc_fill_NA.tsv

    >>>

    output {
         File fill_NA_tsv = "${sample_id}_voc_fill_NA.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task reformat_tsv {
    input {
        File tsv
        String sample_id
    }

    command <<<

        # combine the rows based on matching nucl location

        awk '{f2[$1]=f2[$1] sep[$1] $2;
                f3[$1]=f3[$1] sep[$1] $3;
                f4[$1]=f4[$1] sep[$1] $4;
                f5[$1]=f5[$1] sep[$1] $5;
                f6[$1]=f6[$1] sep[$1] $6;
                f7[$1]=f7[$1] sep[$1] $7;
                sep[$1]="|"}
        END {for(k in f2) print k,f2[k],f3[k],f4[k],f5[k],f6[k],f7[k]}' ~{tsv} > ~{sample_id}_voc_mutations_temp4.tsv

        #fix delimiters and add a column containing the sample id

        sed 's/ /\t/g' ~{sample_id}_voc_mutations_temp4.tsv | awk 'NF=NF+1{$NF="~{sample_id}"}1' > ~{sample_id}_voc_mutations_temp5.tsv

        # fix the column headers, convert from space to tab delimited and then sort by col1
        echo -e "CHROMPOS ~{sample_id}_REF ~{sample_id}_ALT ~{sample_id}_DP ~{sample_id}_RefDP ~{sample_id}_AltDP ~{sample_id}_AltFEQ sample_id" | cat - ~{sample_id}_voc_mutations_temp5.tsv | sed 's/ /\t/g' | sort -t $'\t' -k 1,1 -V > ~{sample_id}_voc_reformat.tsv

    >>>

    output {
         File reformatted_tsv = "${sample_id}_voc_reformat.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task summary_prep {
    input {
        File tsv
        String sample_id
        File voc_annotations
    }

    command <<<

      # cut the columns we want for the results summary (Alt allele and frequency) and make output file
      cut -f3,7 ~{tsv} > ~{sample_id}_voc_mutations_forsummary.tsv

      # cut the columns we want for the counts summary
      awk '{print $8 "\t" $3 "\t" $4 "\t" $6}' ~{tsv} > ~{sample_id}_voc_mutations_temp6.tsv

      # add annotations to the counts summary, reorder the columns, fix the column headers and make output file
      paste ~{voc_annotations} ~{sample_id}_voc_mutations_temp6.tsv | awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7}' | awk 'BEGIN{FS=OFS="\t"; print "sample_id", "AA_change", "Nucl_change", "Lineages", "ALT", "Total_count", "ALT_count"} NR>1{print $1, $2, $3, $4, $5, $6, $7}' > ~{sample_id}_voc_mutations_counts.tsv

   >>>

    output {
         File sample_voc_tsv_summary = "${sample_id}_voc_mutations_forsummary.tsv"
         File sample_voc_tsv_counts = "${sample_id}_voc_mutations_counts.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task demix_reformat {
    input {
        File tsv
        String sample_id
    }

    command <<<

        #datamash to transpose the tsv
        datamash -H transpose < ~{tsv} > ~{sample_id}_demixed_temp1.tsv

        #remove unnecessary characters, change spaces to commas
        tr -d \[ < ~{sample_id}_demixed_temp1.tsv | tr -d \] | tr -d \( | tr -d \) | tr -d \, | tr -d \' | tr ' ' ',' > ~{sample_id}_demixed_temp2.tsv

        # cut and keep the lineage and abundance columns (col3 and col4), fix column headers and delimiters
        awk '{print $3 "\t" $4}' ~{sample_id}_demixed_temp2.tsv | sed -e '1s/abundances/~{sample_id}_lineages/' -e '1s/resid/~{sample_id}_abundances/' | sed 's/ /\t/g' > ~{sample_id}_demixed_temp3.tsv

        # transpose
        datamash -H transpose < ~{sample_id}_demixed_temp3.tsv > ~{sample_id}_demixed_temp4.tsv

        #convert commas to tabs
        tr ',' '\t' < ~{sample_id}_demixed_temp4.tsv > ~{sample_id}_demixed_reformatted.tsv

   >>>

    output {
         File demix_reformatted = "${sample_id}_demixed_reformatted.tsv"
    }

    runtime {
        docker: "rapatsky/debian"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task combine_tsv {
    input {
        Array[File] tsv
        Array[File] tsv_counts
        Array[File] demix_tsv
        File voc_annotations
    }

    command <<<

      # concatenate the tsvs and make the counts tsv summary output
      awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' tsv_counts} >> voc_mutations_counts.tsv

      # fix delimiters in annotations file
      sed 's/ /\t/g' ~{voc_annotations} > voc_annotations.tsv

      # concatentate tsvs for allele frequency summary file and make output
      paste voc_annotations.tsv ~{sep=' ' tsv} > voc_mutations_summary_temp.tsv
      
      #concatenate reformatted demix tsv files into lineages and abundances summary
      cat ~{sep=' ' demix_tsv} > lineage_abundances.tsv

    >>>

    output {
        File voc_summary_temp = "voc_mutations_summary_temp.tsv"
        File voc_counts = "voc_mutations_counts.tsv"
        File demix_summary = "lineage_abundances.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}

task summary_tsv {
    input {
        File tsv
    }

    command <<<

        # datamash to tranpose results summary
        datamash -H transpose < ~{tsv} > voc_mutations_summary.tsv

    >>>

    output {
        File voc_summary = "voc_mutations_summary.tsv"
    }

    runtime {
        docker: "rapatsky/debian"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}

task transfer_outputs {
    input {
        Array[File] variants
        Array[File] depth
        Array[File] demix
        Array[File] sample_voc_tsv_summary
        Array[File] sample_voc_tsv_counts
        File voc_summary
        File voc_counts
        File demix_aggregated
        File demix_summary
        File wwt_variant_abundances_long_format
        File wwt_variant_abundances_long_format_mean
        String out_dir

    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<

        gsutil -m cp ~{sep=' ' variants} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' depth} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' demix} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' sample_voc_tsv_summary} ~{outdirpath}/waste_water_variant_calling/sample_variants/
        gsutil -m cp ~{sep=' ' sample_voc_tsv_counts} ~{outdirpath}/waste_water_variant_calling/sample_variants/
        gsutil -m cp ~{voc_summary} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{voc_counts} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{demix_aggregated} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{demix_summary} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{wwt_variant_abundances_long_format} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{wwt_variant_abundances_long_format_mean} ~{outdirpath}/waste_water_variant_calling/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}
