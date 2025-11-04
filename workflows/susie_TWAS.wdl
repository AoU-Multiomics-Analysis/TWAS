version 1.0 


task SusieTWAS {
    input {
        File LDMatrix 
        File SumStats 
        File SumStatsIndex 
        File FineMapping
        String PhenotypeID
        }

    command <<<
    Rscript /tmp/susie_TWAS.R \
        --PhenotypeID ~{PhenotypeID} \
        --LDMatrix ~{LDMatrix} \
        --SusieRes ~{FineMapping} \
        --SummaryStats ~{SumStats}
    >>>

    output {
        File OutTWAS = "~{PhenotypeID}.TWAS.txt" 
        }


    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk 100 HDD"
    }
}


workflow TWAS {
    input {
        File LDMatrix 
        File Sumstats
        File SumsatsIndex
        File FineMapping
        String PhenotypeID
    }
    call ComputeLD {
        input:
            LDMatrix = LDMatrix,
            PhenotypeID = PhenotypeID,
            SumStats = SumStats,
            SumStatsIndex = SumStatsIndex, 
    }
    output {
        File ResTWAS = SusieTWAS.OutTWAS 
    }

}
