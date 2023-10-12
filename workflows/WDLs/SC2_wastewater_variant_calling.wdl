version 1.0

workflow SC2_wastewater_variant_calling {

    input {

        Array[File] trimsort_bam
        Array[String] sample_name
        Array[String] out_dir_array
        # Array[String] project_name_array

        # reference files/workspace data
        File covid_genome
        File covid_gff
        File covid_gff3

        # Python script
        File ivar_variants_to_vcf

    }
    # secret variables
    # String project_name = project_name_array[0]
    String out_dir = out_dir_array[0]


    scatter (id_bam in zip(sample_name, trimsort_bam)) {
        call add_RG {
            input:
                sample_name = id_bam.left,
                bam = id_bam.right
        }


        call variant_calling {
            input:
                bam = add_RG.rgbam,
                ref = covid_genome,
                ref_gff = covid_gff,
                sample_name = id_bam.left

        }

        call ivar_to_vcf {
            input:
                variants_ivar = variant_calling.variants_ivar,
                ivar_variants_to_vcf = ivar_variants_to_vcf,
                sample_name = id_bam.left
        }

        call csq_and_mutations {
            input:
                vcf = ivar_to_vcf.vcf,
                ref = covid_genome,
                ref_gff3 = covid_gff3,
                sample_name = id_bam.left
        }

        call freyja_demix {
            input:
                variants = csq_and_mutations.variants_bcftools,
                depth = variant_calling.depth,
                sample_name = id_bam.left
        }
    }

    call freyja_aggregate {
        input:
            demix = freyja_demix.demix
    }

    call combine_mutations_tsv {
        input:
            mutations_tsv = csq_and_mutations.mutations_tsv
    }
    
    call transfer_outputs {
        input:
            variants = csq_and_mutations.variants_bcftools,
            depth = variant_calling.depth,
            demix = freyja_demix.demix,
            combined_mutations_tsv = combine_mutations_tsv.combined_mutations_tsv,
            demix_aggregated = freyja_aggregate.demix_aggregated,
            out_dir = out_dir
    }

    output {
        Array[File] addrg_bam = add_RG.rgbam
        Array[File] variants = csq_and_mutations.variants_bcftools
        Array[File] depth = variant_calling.depth
        Array[File] demix = freyja_demix.demix
        File demix_aggregated = freyja_aggregate.demix_aggregated
        File combined_mutations_tsv = combine_mutations_tsv.combined_mutations_tsv
        String transfer_date = transfer_outputs.transfer_date
    }
}

task add_RG {
    input {
        String sample_name
        File bam
    }

    command <<<

        samtools addreplacerg -r ID:~{sample_name} -r LB:L1 -r SM:~{sample_name} -o ~{sample_name}_addRG.bam ~{bam}

    >>>

    output {
        File rgbam = "${sample_name}_addRG.bam"
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
        File bam
        File ref
        File ref_gff
        String sample_name
    }

    command <<<

    samtools mpileup -A -aa -d 600000 -B -Q 20 -q 0 -f ~{ref} ~{bam} | tee >(cut -f1-4 > ~{sample_name}_depth.tsv) | \
    ivar variants -p ~{sample_name}_variants.tsv -q 20 -t 0.0 -r ~{ref} -g ~{ref_gff}
    
    >>>

    output {
        File variants_ivar = "~{sample_name}_variants.tsv"
        File depth = "~{sample_name}_depth.tsv"

    }

     runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/ivar:1.4.2"
    } 
}

task ivar_to_vcf {
    input {
        File variants_ivar
        File ivar_variants_to_vcf
        String sample_name
    }

    command <<<

        python ~{ivar_variants_to_vcf} ~{variants_ivar} ~{sample_name}.vcf

    >>>

    output {
        File vcf = "${sample_name}.vcf"
    }

    runtime {

      docker: "mchether/py3-bio:v1"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }
}

task csq_and_mutations {
    input {
        File vcf
        File ref
        File ref_gff3
        String sample_name
    }

    command <<<

        bcftools csq -f ~{ref} -g ~{ref_gff3} ~{vcf} -Ov -o ~{sample_name}_variants.vcf
        bcftools query -f '[%CHROM\t%SAMPLE\t%POS\t%REF\t%ALT\t%REF_DP\t%ALT_DP\t%AF\t%DP\t%FILTER\t%INFO/GFF_FEATURE\t%INFO/REF_CODON\t%INFO/REF_AA\t%INFO/ALT_CODON\t%INFO/ALT_AA\t%INFO/POS_AA\t%TBCSQ\n]' ~{sample_name}_variants.vcf > ~{sample_name}_variants_temp.tsv
        echo -e "ref_genome\tsample_name\tposition\tref_nucl\talt_nucl\tref_dp\talt_dp\talt_freq\ttotal_dp\tpass\tgff_feature\tref_codon\tref_aa\talt_codon\talt_aa\tposition_aa\tbcsq" | cat - ~{sample_name}_variants_temp.tsv > ~{sample_name}_mutations.tsv

    >>>

    output {
        File variants_bcftools = "${sample_name}_variants.vcf"
        File mutations_tsv = "${sample_name}_mutations.tsv"
    }

    runtime {

        docker: "staphb/bcftools:1.18"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"

    }
}


task freyja_demix {
    input {
        String sample_name
        File variants
        File depth
    }

    command <<<

        #get updated lineages for demixing
        mkdir ./freyja_db
        freyja update --outdir ./freyja_db
        
        #creates a temp file with the same name as the intended output file that will get output in case of failure or overwritten in case of sucess
        echo -e "\t~{sample_name}\nsummarized\tLowCov\nlineages\tLowCov\nabundances\tLowCov\nresid\tLowCov\ncoverage\tLowCov" > ~{sample_name}_demixed.tsv
        
        freyja demix --eps 0.01 --covcut 10 --barcodes ./freyja_db/usher_barcodes.csv --meta ./freyja_db/curated_lineages.json --confirmedonly ~{variants} ~{depth} --output ~{sample_name}_demixed.tsv

    >>>

    output {
        File demix = "${sample_name}_demixed.tsv"
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

task combine_mutations_tsv {
    input {
        Array[File] mutations_tsv
    }

    command <<<
        
        # combine the counts and frequency files for all samples into one
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' mutations_tsv} >> combined_mutations.tsv
    
    >>>

    output {
        File combined_mutations_tsv = "combined_mutations.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 500 HDD"
    }
}

task transfer_outputs {
    input {
        Array[File] variants
        Array[File] depth
        Array[File] demix
        File demix_aggregated
        File combined_mutations_tsv
        String out_dir

    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<

        gsutil -m cp ~{sep=' ' variants} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' depth} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' demix} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{demix_aggregated} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{combined_mutations_tsv} ~{outdirpath}/waste_water_variant_calling/
    

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
