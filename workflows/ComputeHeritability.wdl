version 1.0 



task ComputeHeritabilityPlink {
    input {
        File pvar
        File pgen
        File psam
        File ExpressionBed 
        File Covars 
        File PhenotypeID
    }

    command <<<

    Rscript /tmp/CreatePhenoFile.R \
        --BedFile ~{ExpressionBed} \
        --PhenotypeID ~{PhenotypeID}
    
    IFS=$'\t' read CHR START END < <(awk 'NR==2 {print $1, $2, $3}' your_file.txt)   

    # subset pfiles to cis window for each gene  
    plink2 \
        --pgen ~{pgen} \
        --pvar ~{pvar} \
        --psam ~{psam} \
        --chr $CHR  \
        --from-bp $START \
        --to-bp $END \
        --pheno pheno.txt \ 
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
        --qcovars ~{Covars} \
        --out ~{PhenotypeID}

    >>>
    
    output {

    }
    
    runtime {

    }
}


workflow EstimateHeritability {
    input {
        File pvar
        File pgen
        File psam
        File ExpressionBed 
        File Covars 
        File PhenotypeID
    }
    
    call ComputeHeritabilityPlink {
        input: 
            pvar = pvar,
            pgen = pgen,
            psam = psam,
            ExpressionBed = ExpressionBed,
            Covars = Covars,
            PhenotypeID = PhenotypeID
    } 

    output {

    }

}
