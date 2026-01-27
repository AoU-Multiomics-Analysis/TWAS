suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(arrow)
    library(bedr)
    })

######### FUNCTIONS ######################

load_susie_data <- function(path) {
susie_data <- arrow::read_parquet(path) %>% 
        select(variant,posterior_mean,chromosome,position,ref,alt,pip,cs_id)
susie_data
    
}




extract_interval <- function(finemapping_data) {
    start <- min(as.numeric(finemapping_data$position))
    end <- max(as.numeric(finemapping_data$position))
    chromosome <- unique(finemapping_data$chromosome)[1]
    interval <- paste0(chromosome,':',start,'-',end)
    interval
}
# take variants from finemapping results and query GWAS file
# returns merged dataframe with fine-mapped variants and 
# gwas summary statistics
extract_gwas_data <- function(input_susie,input_gwas) {
    require(bedr)
    header <- strsplit(readLines(input_gwas, n = 1),"\t")[[1]]
    gwas_dat <- tabix(
                    extract_interval(input_susie),
                    input_gwas,
                    check.chr = FALSE
                )
    colnames(gwas_dat) <- header
    cleaned_gwas_dat <- gwas_dat %>% 
                mutate(
                        POS = as.numeric(POS),
                        Pvalue = as.numeric(Pvalue),
                        BETA = as.numeric(BETA),
                        SE = as.numeric(SE)
                    )  
    output <- input_susie %>% 
                mutate(position = as.numeric(position)) %>% 
                left_join(cleaned_gwas_dat,by = c('position' = 'POS','chromosome' ='CHR')) %>% 
                mutate(allele_match = case_when(ref == REF & alt == ALT ~ TRUE,
                                              ref == REF & alt == ALT ~ FALSE))  
    output
}

calculate_TWAS_Z <- function(variant_df,LD_matrix) {
    subset_LD <- LD_matrix[variant_df$variant,variant_df$variant]
    denom <- t(variant_df$posterior_mean) %*% data.matrix(subset_LD) %*%  variant_df$posterior_mean
    stat <- t(variant_df$posterior_mean) %*%  variant_df$z

    zscore <- stat/sqrt(denom)
    pvalue <- pchisq(zscore * zscore,1,lower.tail = FALSE)
    output <- data.frame(zscore = zscore,pvalue = pvalue,stat = stat) 
    output
}

empty_twas_row <- function(gene_id, gwas_name = NA_character_) {
  tibble(
    gene = gene_id,
    GWAS = gwas_name,
    zscore = NA_real_,
    pvalue = NA_real_
  )
}
TWAS <- function(GWAS_path,SusieData,LD) {
message('Extracting GWAS summary statistics for:')
message(GWAS_path)
PhenotypeID <- SusieData %>% distinct(geneID) %>% pull(geneID)
message(PhenotypeID)
  GWAS <- tryCatch(
    extract_gwas_data(SusieData, GWAS_path),
    error = function(e) {
      message("GWAS extraction failed for ", PhenotypeID, " : ", e$message)
      return(NULL)
    }
  )

  if (is.null(GWAS) || nrow(GWAS) == 0) {
    return(
      empty_twas_row(
        gene_id = PhenotypeID,
        gwas_name = tools::file_path_sans_ext(basename(GWAS_path))
      )
    )
  }
GeneLD <- LD[[PhenotypeID]]
    
# filters gwas data based on number of measurements 
# for each variant that are present in gwas data. 
# Sometimes variants have multiple measurements in a gwas 
# but im not qutie sure why
FilteredGWAS <- GWAS %>%
    filter(allele_match == TRUE) %>% 
    mutate(z = BETA/SE) %>% 
    group_by(variant) %>% 
    filter(dplyr::n() == 1) %>% 
    ungroup()

  if (nrow(FilteredGWAS) == 0) {
    return(
      empty_twas_row(
        gene_id = PhenotypeID,
        gwas_name = tools::file_path_sans_ext(basename(GWAS_path))
      )
    )
  }
#message('Computing TWAS Z')
ResTWAS <- FilteredGWAS %>% 
            calculate_TWAS_Z(GeneLD) %>% 
            mutate(gene = PhenotypeID,
                   GWAS = tools::file_path_sans_ext(basename(GWAS_path))
                    )

ResTWAS
}
######### PARSE COMMAND LINE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--LDMatrix"), type="character", default=NULL,
                        help="Path to RDS object that is a list and contains a matrix for each gene in analysis", metavar = "type"),
    optparse::make_option(c("--SummaryStats"), type="character", default=NULL,
                        help="path to indexed summary stats", metavar = "type"),
    optparse::make_option(c("--SusieRes"), type="character", default=NULL,
                        help="Path to finemapping data for a gene", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Path to finemapping data for a gene", metavar = "type")

    )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
MatrixLD <- opt$LDMatrix
FineMappingRes <- opt$SusieRes
SummaryStats <- opt$SummaryStats
OutFileName <- paste0(opt$OutputPrefix,'.TWAS.txt')
PhenotypeID <- opt$PhenotypeID


##################### LOAD DATA ##########################
# loading finemapping data 
message('Loading fine mapping')
susie_dat <- fread(FineMappingRes) %>% 
                mutate(variant = str_replace(variant,'chrchr','chr'))

# Loads LD matrix
message('Loading LD matrix')
LD <- readRDS(MatrixLD)

############# COMPUTE TWAS Z SCORE ########################
# computes TWAS Z statistic 
#ResTWAS <-  susie_dat %>%
            #TWAS(SummaryStats,LD,PhenotypeID)

# loop over summary stats and perform TWAS 
# for each phenotype 
#ResTWAS <-  SummaryStatsList %>%
                #map_dfr(~TWAS(.x, SusieData = susie_dat,LD = LD,PhenotypeID = PhenotypeID))

ResTWAS <- susie_dat %>% 
    filter(molecular_trait_id %in% names(LDobj)) %>%    
    mutate(geneID = molecular_trait_id) %>% 
    group_by(molecular_trait_id) %>% 
    group_modify(~TWAS(SummaryStats,.,LDobj))
ResTWAS %>% write_tsv(OutFileName)
