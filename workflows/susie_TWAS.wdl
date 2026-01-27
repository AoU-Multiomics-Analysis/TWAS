version 1.0 


task SusieTWAS {
    input {
        File LDMatrix 
        File SumStats 
        File SumStatsIndex
        String NameGWAS
        File FineMapping
        Int NumPrempt
        Int Memory
        }

    command <<<
    Rscript /tmp/susie_TWAS.R \
        --LDMatrix ~{LDMatrix} \
        --SusieRes ~{FineMapping} \
        --OutputPrefix ~{NameGWAS} \
        --SummaryStats ~{SumStats}
    >>>

    output {
        File OutTWAS = "~{NameGWAS}.TWAS.txt" 
        }


    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        preemptible: "~{NumPrempt}"
        cpu: "4"
        memory: "~{Memory} GB"
        disks: "local-disk 100 HDD"
    }
}


workflow TWAS {
    input {
        File LDMatrix 
        File SumStats
        File SumStatsIndex
        File FineMapping
        String NameGWAS
        Int NumPrempt
        Int Memory

    }
    call SusieTWAS {
        input:
            LDMatrix = LDMatrix,
            FineMapping = FineMapping,
            SumStats = SumStats,
            SumStatsIndex = SumStatsIndex,
            NumPrempt = NumPrempt,
            NameGWAS = NameGWAS,
            Memory = Memory
    }
    output {
        File ResTWAS = SusieTWAS.OutTWAS 
    }

}
