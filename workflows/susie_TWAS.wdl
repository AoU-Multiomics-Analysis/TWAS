version 1.0 


task SusieTWAS {
    input {
        File LDMatrix 
        File SumStats 
        File SumStatsIndex
        String NameGWAS
        File FineMapping
        File? GeneFilter
        String ChosenLabel = ""
        String QTLType = ""
        Int NumPrempt
        Int Memory
        }

    command <<<
    Rscript /tmp/susie_TWAS.R \
        --LDMatrix ~{LDMatrix} \
        --SusieRes ~{FineMapping} \
        --OutputPrefix ~{NameGWAS} \
        --SummaryStats ~{SumStats} \
        --ChosenLabel "~{ChosenLabel}" \
        --QTLType "~{QTLType}" \
        --GeneFilter '~{if defined(GeneFilter) then select_first([GeneFilter]) else ""}'
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
        File? GeneFilter
        String ChosenLabel = ""
        String QTLType = ""
        Int NumPrempt
        Int Memory

    }
    call SusieTWAS {
        input:
            LDMatrix = LDMatrix,
            FineMapping = FineMapping,
            SumStats = SumStats,
            SumStatsIndex = SumStatsIndex,
            GeneFilter = GeneFilter,
            ChosenLabel = ChosenLabel,
            QTLType = QTLType,
            NumPrempt = NumPrempt,
            NameGWAS = NameGWAS,
            Memory = Memory
    }
    output {
        File ResTWAS = SusieTWAS.OutTWAS 
    }

}
