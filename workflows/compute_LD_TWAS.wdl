version 1.0 


task ComputeLD {
    input {
        File DoseMatrix
        String PhenotypeID
        }

    command <<<
    Rscript TWAS_compute_LD.R \
        --PhenotypeID ~{PhenotypeID} \
        --DoseMatrix ~{DoseMatrix}
    >>>

    output {
        File MatrixLD = "~{PhenotypeID}.LD.rds" 
        }


    runtime {
        docker: "quay.io/biocontainers/bcftools:1.22--h3a4d415_0"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk 100 HDD"
    }
}


workflow PreprocessLD {
    input {
        File DoseMatrix
        String PhenotypeID
    }
    call ComputeLD {
        input:
            DoseMatrix = DoseMatrix,
            PhenotypeID = PhenotypeID
    }


}
