version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture

workflow SC2_wastewater_variant_calling {

    input {

        Array[File] trimsort_bam
        Array[String] sample_name
        Array[String] out_dir_array
        Array[String] project_name_array

        # reference files/workspace data
        File covid_genome
        File covid_gff

        # python scripts
        File version_capture_wwt_variant_calling_py

    }
    # secret variables
    String project_name = project_name_array[0]
    String out_dir = out_dir_array[0]
    String outdirpath = sub(out_dir, "/$", "")


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

        call freyja_demix {
            input:
                variants = variant_calling.variants,
                depth = variant_calling.depth,
                sample_name = id_bam.left
        }
        
        call mutations_tsv {
            input:
                variants = variant_calling.variants,
                sample_name = id_bam.left
        }
    }

    call freyja_aggregate {
        input:
            demix = freyja_demix.demix
    }

    call combine_mutations_tsv {
        input:
            mutations_tsv = mutations_tsv.mutations_tsv
    }
    
    call version_capture.workflow_version_capture as workflow_version_capture {
        input:
    }
    
    call create_version_capture_file {
        input:
            version_capture_wwt_variant_calling_py = version_capture_wwt_variant_calling_py
            project_name = project_name,
            samtools_version_staphb = select_all(add_RG.samtools_version_staphb)[0],
            samtools_version_andersenlabapps = select_all(variant_calling.samtools_version_andersenlabapps)[0],
            ivar_version = select_all(variant_calling.ivar_version)[0],
            freyja_version = select_all(freyja_demix.freyja_version)[0],
            analysis_date = workflow_version_capture.analysis_date,
            workflow_version = workflow_version_capture.workflow_version
            

    }

    call transfer_outputs {
        input:
            variants = variant_calling.variants,
            depth = variant_calling.depth,
            demix = freyja_demix.demix,
            combined_mutations_tsv = combine_mutations_tsv.combined_mutations_tsv,
            demix_aggregated = freyja_aggregate.demix_aggregated,
            outdirpath = outdirpath,
            version_capture_wwt_variant_calling = create_version_capture_file.version_capture_wwt_variant_calling
    }

    output {
        Array[File] addrg_bam = add_RG.rgbam
        Array[File] variants = variant_calling.variants
        Array[File] depth = variant_calling.depth
        Array[File] demix = freyja_demix.demix
        File demix_aggregated = freyja_aggregate.demix_aggregated
        File combined_mutations_tsv = combine_mutations_tsv.combined_mutations_tsv
        File version_capture_wwt_variant_calling = create_version_capture_file.version_capture_wwt_variant_calling
        String transfer_date_wwt_variant_calling = transfer_outputs.transfer_date_wwt_variant_calling
    }
}

task add_RG {
    input {
        String sample_name
        File bam
    }

    command <<<

        # grab samtools version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION

        samtools addreplacerg -r ID:~{sample_name} -r LB:L1 -r SM:~{sample_name} -o ~{sample_name}_addRG.bam ~{bam}

    >>>

    output {
        File rgbam = "${sample_name}_addRG.bam"
        String samtools_version_staphb = read_string("VERSION")
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

    # grab ivar and samtools versions
    ivar version | awk '/version/ {print $3}' | tee VERSION_ivar
    samtools --version | awk '/samtools/ {print $2}' | tee VERSION_samtools

    samtools mpileup -A -aa -d 600000 -B -Q 20 -q 0 -f ~{ref} ~{bam} | tee >(cut -f1-4 > ~{sample_name}_depth.tsv) | \
    ivar variants -p ~{sample_name}_variants.tsv -q 20 -t 0.0 -r ~{ref} -g ~{ref_gff}
    
    >>>

    output {
        File variants = "~{sample_name}_variants.tsv"
        File depth = "~{sample_name}_depth.tsv"
        String samtools_version_andersenlabapps = read_string("VERSION_samtools")
        String ivar_version = read_string("VERSION_ivar")

    }

     runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    } 
}

task freyja_demix {
    input {
        String sample_name
        File variants
        File depth
    }

    command <<<

        freyja --version | tee VERSION

        #get updated lineages for demixing
        mkdir ./freyja_db
        freyja update --outdir ./freyja_db
        
        #creates a temp file with the same name as the intended output file that will get output in case of failure or overwritten in case of sucess
        echo -e "\t~{sample_name}\nsummarized\tLowCov\nlineages\tLowCov\nabundances\tLowCov\nresid\tLowCov\ncoverage\tLowCov" > ~{sample_name}_demixed.tsv
        
        freyja demix --eps 0.01 --covcut 10 --barcodes ./freyja_db/usher_barcodes.csv --meta ./freyja_db/curated_lineages.json --confirmedonly ~{variants} ~{depth} --output ~{sample_name}_demixed.tsv

    >>>

    output {
        File demix = "${sample_name}_demixed.tsv"
        String freyja_version = read_string("VERSION")
    }

    runtime {
        docker: "staphb/freyja:latest"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
        continueOnReturnCode: [0, 1]
    }
}

task mutations_tsv {
    input {
        String sample_name
        File variants
    }

    command <<<
        
        #add a column that is filled in with the sample name
        sed 's/$/\t~{sample_name}/' ~{variants} > ~{sample_name}_variants_temp.tsv

        #cut columns needed for allele counts and frequency output and then rename the column headers
        awk '{print($20,"\t",$1,"\t",$2,"\t",$3,"\t",$4,"\t",$17,"\t",$19,"\t",$5,"\t",$8,"\t",$11)}' ~{sample_name}_variants_temp.tsv | sed -e '1s/~{sample_name}/sample_name/' -e '1s/REGION/ref_genome/' -e '1s/POS/position/' -e '1s/REF/ref_nucl/' -e '1s/ALT/alt_nucl/' -e '1s/REF_AA/ref_aa/' -e '1s/ALT_AA/alt_aa/' -e '1s/REF_DP/ref_depth/' -e '1s/ALT_DP/alt_depth/' -e '1s/ALT_FREQ/alt_freq/' > ~{sample_name}_mutations.tsv

    >>>

    output {
        File mutations_tsv = "${sample_name}_mutations.tsv"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 500 HDD"
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
        docker: "staphb/freyja:latest"
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
        
        # combine the coutns and frequency files for all samples into one
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

task create_version_capture_file {
    input {
        File version_capture_wwt_variant_calling_py
        String project_name
        String samtools_version_staphb
        String samtools_version_andersenlabapps
        String ivar_version
        String freyja_version
        String analysis_date
        String workflow_version
    }

    command <<<
    
        python ~{version_capture_wwt_variant_calling_py} \
        --project_name "~{project_name}" \
        --samtools_version_staphb "~{samtools_version_staphb}" \
        --samtools_version_andersenlabapps "~{samtools_version_andersenlabapps}" \
        --ivar_version "~{ivar_version}" \
        --freyja_version "~{freyja_version}" \
        --analysis_date "~{analysis_date}" \
        --workflow_version "~{workflow_version}"

    >>>

    output {
        File version_capture_wwt_variant_calling = 'version_capture_wwt_variant_calling_~{project_name}_v~{workflow_version}.csv'
    }

    runtime {

      docker: "mchether/py3-bio:v4"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }
}

task transfer_outputs {
    input {
        Array[File] variants
        Array[File] depth
        Array[File] demix
        File demix_aggregated
        File combined_mutations_tsv
        File version_capture_wwt_variant_calling
        String outdirpath

    }


    command <<<

        gsutil -m cp ~{sep=' ' variants} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' depth} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{sep=' ' demix} ~{outdirpath}/waste_water_variant_calling/freyja/
        gsutil -m cp ~{demix_aggregated} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{combined_mutations_tsv} ~{outdirpath}/waste_water_variant_calling/
        gsutil -m cp ~{version_capture_wwt_variant_calling} ~{outdirpath}/summary_results/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date_wwt_variant_calling = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 50 SSD"
    }
}
