version 1.0 


task ComputeLD {
    input {
        File DoseMatrix
        File VariantList
        String PhenotypeID
        Int NumPrempt
        Int Memory
        }

    command <<<
    Rscript /tmp/compute_LD_TWAS.R \
        --PhenotypeID ~{PhenotypeID} \
        --DoseMatrix ~{DoseMatrix} \
        --VariantList ~{VariantList}
    >>>

    output {
        File MatrixLD = "~{PhenotypeID}.LD.rds" 
        }


    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        cpu: "4"
        preemptible: "${NumPrempt}"
        memory: "${Memory} GB"
        disks: "local-disk 100 HDD"
    }
}


workflow PreprocessLD {
    input {
        File DoseMatrix
        File VariantList
        String PhenotypeID
        Int NumPrempt
        Int Memory
    }
    call ComputeLD {
        input:
            DoseMatrix = DoseMatrix,
            VariantList = VariantList,
            PhenotypeID = PhenotypeID,
            NumPrempt = NumPrempt,
            Memory = Memory
    }
    output {
        File MatrixLD = ComputeLD.MatrixLD

    }

}
