version 1.0 



task ComputeHeritabilityPlink {
    input {
        File pvar
        File pgen
        File psam
        File ExpressionBed 
        File Covars 
        String PhenotypeID
        Int Memory
        Int NumPrempt
    }

    command <<<

    Rscript /tmp/CreatePhenoFile.R \
        --BedFile ~{ExpressionBed} \
        --PhenotypeID ~{PhenotypeID}
    
    IFS=$'\t' read CHR START END < <(awk 'BEGIN{OFS="\t"} NR==2 {print $1, $2, $3}' gene_region.tsv)
    # subset pfiles to cis window for each gene  
    plink2 \
        --pgen ~{pgen} \
        --pvar ~{pvar} \
        --psam ~{psam} \
        --chr $CHR  \
        --from-bp $START \
        --to-bp $END \
        --make-bed \
        --out ~{PhenotypeID}
    
    # create GRM for cis window 
    gcta64 \
        --bfile ~{PhenotypeID} \
        --make-grm \
        --out ~{PhenotypeID} 
    
    # compute heritability estimate
    gcta64 \
        --grm ~{PhenotypeID} \
        --pheno pheno.txt \
        --reml \
        --qcovar ~{Covars} \
        --out ~{PhenotypeID}

    >>>
    
    output {
        File H2Estimate = "~{PhenotypeID}.hsq"
    }
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/twas:main"
        cpu: "4"
        preemptible: "${NumPrempt}"
        memory: "${Memory} GB"
        disks: "local-disk 100 HDD"
    }
}


workflow EstimateHeritability {
    input {
        File pvar
        File pgen
        File psam
        File ExpressionBed 
        File Covars 
        String PhenotypeID
        Int Memory 
        Int NumPrempt
    }
    
    call ComputeHeritabilityPlink {
        input: 
            pvar = pvar,
            pgen = pgen,
            psam = psam,
            ExpressionBed = ExpressionBed,
            Covars = Covars,
            PhenotypeID = PhenotypeID,
            Memory = Memory,
            NumPrempt = NumPrempt
    } 

    output {
        File HeritabilityEstimate = ComputeHeritabilityPlink.H2Estimate 
    }

}
