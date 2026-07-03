version 1.0

task PrepareTWASManifest {
    input {
        File QTLManifest
        File GWASManifest
        String OutputPrefix
    }

    command <<<
    Rscript /tmp/prepare_TWAS_manifest.R \
        --QTLManifest "~{QTLManifest}" \
        --GWASManifest "~{GWASManifest}" \
        --OutputPrefix "~{OutputPrefix}"
    >>>

    output {
        File PairManifest = "twas_manifest_pairs.tsv"
        Array[String] PairIds = read_lines("pair_id.txt")
        Array[String] SusiePaths = read_lines("susie_path.txt")
        Array[String] LDMatrixPaths = read_lines("ld_matrix_path.txt")
        Array[String] QTLTypes = read_lines("qtl_type.txt")
        Array[String] QTLTypeSafe = read_lines("qtl_type_safe.txt")
        Array[String] SummaryStatsPaths = read_lines("summary_stats_path.txt")
        Array[String] SummaryStatsIndexPaths = read_lines("summary_stats_index_path.txt")
        Array[String] ChosenLabels = read_lines("chosen_label.txt")
        Array[String] ChosenLabelSafe = read_lines("chosen_label_safe.txt")
        Array[String] OutputPrefixes = read_lines("output_prefix.txt")
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        cpu: "1"
        memory: "4 GB"
        disks: "local-disk 10 HDD"
    }
}

task ManifestTWASRun {
    input {
        String PairId
        String FineMappingPath
        String LDMatrixPath
        String QTLType
        String QTLTypeSafe
        String SummaryStatsPath
        String SummaryStatsIndexPath
        String ChosenLabel
        String ChosenLabelSafe
        String OutputPrefix
        String AnalysisType = "conditional_rare"
        String RareColumn = "rare"
        File? GeneFilter
        Int NumPrempt
        Int Memory
        Int DiskGB
    }

    command <<<
    set -euo pipefail

    localize_one() {
        local src="$1"
        local dest="$2"

        if [[ "${src}" == gs://* ]]; then
            gsutil cp "${src}" "${dest}"
        elif [[ "${src}" != "${dest}" ]]; then
            cp "${src}" "${dest}"
        fi
    }

    fine_mapping_src='~{FineMappingPath}'
    ld_matrix_src='~{LDMatrixPath}'
    sumstats_src='~{SummaryStatsPath}'
    sumstats_index_src='~{SummaryStatsIndexPath}'

    fine_mapping_local="$(basename "${fine_mapping_src}")"
    ld_matrix_local="$(basename "${ld_matrix_src}")"
    sumstats_local="$(basename "${sumstats_src}")"
    sumstats_index_local="${sumstats_local}.tbi"

    localize_one "${fine_mapping_src}" "${fine_mapping_local}"
    localize_one "${ld_matrix_src}" "${ld_matrix_local}"
    localize_one "${sumstats_src}" "${sumstats_local}"
    localize_one "${sumstats_index_src}" "${sumstats_index_local}"

    analysis_type='~{AnalysisType}'

    if [[ "${analysis_type}" == "conditional_rare" ]]; then
        Rscript /tmp/conditional_rare_TWAS.R \
            --LDMatrix "${ld_matrix_local}" \
            --SusieRes "${fine_mapping_local}" \
            --OutputPrefix "~{OutputPrefix}" \
            --SummaryStats "${sumstats_local}" \
            --RareColumn "~{RareColumn}" \
            --ChosenLabel "~{ChosenLabel}" \
            --QTLType "~{QTLType}" \
            --GeneFilter '~{if defined(GeneFilter) then select_first([GeneFilter]) else ""}'

        cp "~{OutputPrefix}.TWAS.two_predictor.txt" result.TWAS.tsv
    elif [[ "${analysis_type}" == "standard" ]]; then
        Rscript /tmp/susie_TWAS.R \
            --LDMatrix "${ld_matrix_local}" \
            --SusieRes "${fine_mapping_local}" \
            --OutputPrefix "~{OutputPrefix}" \
            --SummaryStats "${sumstats_local}" \
            --ChosenLabel "~{ChosenLabel}" \
            --QTLType "~{QTLType}" \
            --GeneFilter '~{if defined(GeneFilter) then select_first([GeneFilter]) else ""}'

        cp "~{OutputPrefix}.TWAS.txt" result.TWAS.tsv
    else
        echo "AnalysisType must be 'standard' or 'conditional_rare'. Received: ${analysis_type}" >&2
        exit 1
    fi

    cat > result.metadata.tsv <<'METADATA'
pair_id	qtl_type	qtl_type_safe	chosen_label	chosen_label_safe	output_prefix	analysis_type	susie_path	ld_matrix_path	summary_stats_path	summary_stats_index_path
~{PairId}	~{QTLType}	~{QTLTypeSafe}	~{ChosenLabel}	~{ChosenLabelSafe}	~{OutputPrefix}	~{AnalysisType}	~{FineMappingPath}	~{LDMatrixPath}	~{SummaryStatsPath}	~{SummaryStatsIndexPath}
METADATA
    >>>

    output {
        File OutTWAS = "result.TWAS.tsv"
        File Metadata = "result.metadata.tsv"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        preemptible: "~{NumPrempt}"
        cpu: "4"
        memory: "~{Memory} GB"
        disks: "local-disk ~{DiskGB} HDD"
    }
}

task AggregateManifestTWAS {
    input {
        Array[File] TWASResults
        Array[File] MetadataFiles
        String OutputPrefix
        Int Memory
        Int NumThreads
        Int DiskGB
    }

    command <<<
    set -euo pipefail

    printf '%s\n' ~{sep=' ' TWASResults} > twas_files.txt
    printf '%s\n' ~{sep=' ' MetadataFiles} > metadata_files.txt

    Rscript /tmp/aggregate_manifest_TWAS.R \
        --TWASFiles twas_files.txt \
        --MetadataFiles metadata_files.txt \
        --OutputPrefix "~{OutputPrefix}"
    >>>

    output {
        File AllResults = "~{OutputPrefix}.all_TWAS.tsv"
        Array[File] ByQTLTypeResults = glob("by_qtl_type/*.tsv")
        File RunMetadata = "~{OutputPrefix}.manifest_run_metadata.tsv"
    }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        cpu: "~{NumThreads}"
        memory: "~{Memory} GB"
        disks: "local-disk ~{DiskGB} HDD"
    }
}

workflow ManifestTWAS {
    input {
        File QTLManifest
        File GWASManifest
        String OutputPrefix
        String AnalysisType = "conditional_rare"
        String RareColumn = "rare"
        File? GeneFilter
        Int NumPrempt
        Int TWASMemory
        Int TWASDiskGB = 500
        Int AggregateMemory
        Int AggregateThreads
        Int AggregateDiskGB = 500
    }

    call PrepareTWASManifest {
        input:
            QTLManifest = QTLManifest,
            GWASManifest = GWASManifest,
            OutputPrefix = OutputPrefix
    }

    scatter (pair_index in range(length(PrepareTWASManifest.PairIds))) {
        call ManifestTWASRun {
            input:
                PairId = PrepareTWASManifest.PairIds[pair_index],
                FineMappingPath = PrepareTWASManifest.SusiePaths[pair_index],
                LDMatrixPath = PrepareTWASManifest.LDMatrixPaths[pair_index],
                QTLType = PrepareTWASManifest.QTLTypes[pair_index],
                QTLTypeSafe = PrepareTWASManifest.QTLTypeSafe[pair_index],
                SummaryStatsPath = PrepareTWASManifest.SummaryStatsPaths[pair_index],
                SummaryStatsIndexPath = PrepareTWASManifest.SummaryStatsIndexPaths[pair_index],
                ChosenLabel = PrepareTWASManifest.ChosenLabels[pair_index],
                ChosenLabelSafe = PrepareTWASManifest.ChosenLabelSafe[pair_index],
                OutputPrefix = PrepareTWASManifest.OutputPrefixes[pair_index],
                AnalysisType = AnalysisType,
                RareColumn = RareColumn,
                GeneFilter = GeneFilter,
                NumPrempt = NumPrempt,
                Memory = TWASMemory,
                DiskGB = TWASDiskGB
        }
    }

    call AggregateManifestTWAS {
        input:
            TWASResults = ManifestTWASRun.OutTWAS,
            MetadataFiles = ManifestTWASRun.Metadata,
            OutputPrefix = OutputPrefix,
            Memory = AggregateMemory,
            NumThreads = AggregateThreads,
            DiskGB = AggregateDiskGB
    }

    output {
        File PairManifest = PrepareTWASManifest.PairManifest
        Array[File] PairResults = ManifestTWASRun.OutTWAS
        Array[File] PairMetadata = ManifestTWASRun.Metadata
        File AllResults = AggregateManifestTWAS.AllResults
        Array[File] ByQTLTypeResults = AggregateManifestTWAS.ByQTLTypeResults
        File RunMetadata = AggregateManifestTWAS.RunMetadata
    }
}
