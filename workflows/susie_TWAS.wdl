version 1.0 


task SusieTWAS {
    input {
        File LDMatrix 
        Array[File] SumStats 
        Array[File] SumStatsIndex 
        File FineMapping
        String PhenotypeID
        Int NumPrempt
        }

    command <<<
    Rscript /tmp/susie_TWAS.R \
        --PhenotypeID ~{PhenotypeID} \
        --LDMatrix ~{LDMatrix} \
        --SusieRes ~{FineMapping} \
        --SummaryStats ~{sep="," SumStats}
    >>>

    output {
        File OutTWAS = "~{PhenotypeID}.TWAS.txt" 
        }


    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        preemptible: "~{NumPrempt}"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk 100 HDD"
    }
}


workflow TWAS {
    input {
        File LDMatrix 
        Array[File] SumStats
        Array[File] SumStatsIndex
        File FineMapping
        String PhenotypeID
        Int NumPrempt

    }
    call SusieTWAS {
        input:
            LDMatrix = LDMatrix,
            PhenotypeID = PhenotypeID,
            FineMapping = FineMapping,
            SumStats = SumStats,
            SumStatsIndex = SumStatsIndex,
            NumPrempt = NumPrempt
    }
    output {
        File ResTWAS = SusieTWAS.OutTWAS 
    }

}
