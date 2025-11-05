version 1.0 

task AggregateTWAS {
    input {
        File TWASResFOFN
        Int Memory
        String OutputPrefix
        Int NumThreads
    }
    command <<< 
    export GSUTIL_PARALLEL_PROCESS_COUNT=32
    export GSUTIL_PARALLEL_THREAD_COUNT=8

    awk '{print $1}' ~{TWASResFOFN} | grep -v '^$' > file_paths.txt 

    mkdir -p localized
    gsutil -m cp -I localized/ < file_paths.txt 

    # Write the new local file paths into filelist.txt
    ls -1 "$(pwd)/localized/*" > filelist.txt
    Rscript /tmp/aggregate_TWAS.R --FilePaths filelist.txt  --OutputPrefix ~{OutputPrefix}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        disks: "local-disk 500 SSD"
        memory: "~{Memory}GB"
        cpu: "~{NumThreads}"
    }
 


    output {
        File mergedTWASRes = "~{OutputPrefix}_TWAS.tsv" 
    } 
}


workflow AggregateTWASWorkflow {
    input {
        File TWASResFOFN
        Int Memory
        String OutputPrefix
        Int NumThreads
    }
    call AggregateTWAS {
        input:
            TWASResFOFN = TWASResFOFN,
            Memory = Memory,
            OutputPrefix = OutputPrefix,
            NumThreads = NumThreads
    }
    output {
        File AggregatedTWAS = AggregateTWAS.mergedTWASRes
    }
}


