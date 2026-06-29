version 1.0


task ConditionalRareTWAS {
    input {
        File LDMatrix
        File SumStats
        File SumStatsIndex
        String NameGWAS
        File FineMapping
        String RareColumn = "rare"
        File? GeneFilter
        String ChosenLabel = ""
        String QTLType = ""
        Int NumPrempt
        Int Memory
    }

    command <<<
    Rscript /tmp/conditional_rare_TWAS.R \
        --LDMatrix ~{LDMatrix} \
        --SusieRes ~{FineMapping} \
        --OutputPrefix ~{NameGWAS} \
        --SummaryStats ~{SumStats} \
        --RareColumn ~{RareColumn} \
        --ChosenLabel "~{ChosenLabel}" \
        --QTLType "~{QTLType}" \
        --GeneFilter '~{if defined(GeneFilter) then select_first([GeneFilter]) else ""}'
    >>>

    output {
        File OutTWAS = "~{NameGWAS}.TWAS.two_predictor.txt"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        preemptible: "~{NumPrempt}"
        cpu: "4"
        memory: "~{Memory} GB"
        disks: "local-disk 100 HDD"
    }
}


workflow ConditionalRareTWASAnalysis {
    input {
        File LDMatrix
        File SumStats
        File SumStatsIndex
        File FineMapping
        String NameGWAS
        String RareColumn = "rare"
        File? GeneFilter
        String ChosenLabel = ""
        String QTLType = ""
        Int NumPrempt
        Int Memory
    }

    call ConditionalRareTWAS {
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
            RareColumn = RareColumn,
            Memory = Memory
    }

    output {
        File ResTWAS = ConditionalRareTWAS.OutTWAS
    }
}
