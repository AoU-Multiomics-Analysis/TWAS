library(tidyverse)
library(data.table)
library(arrow)
library(bedr)

load_susie_data <- function(path) {
susie_data <- arrow::read_parquet(path) %>% 
        select(variant,posterior_mean,chromosome,position,ref,alt,pip,cs_id)
susie_data
    
}



extract_interval <- function(finemapping_data) {
    start <- min(finemapping_data$position)
    end <- max(finemapping_data$position)
    chromosome <- unique(finemapping_data$chromosome)[1]
    interval <- paste0(chromosome,':',start,'-',end)
    interval
}

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
                        base_pair_location = as.numeric(base_pair_location),
                        p_value = as.numeric(p_value),
                        beta = as.numeric(beta),
                        standard_error = as.numeric(standard_error)
                    )  
    output <- input_susie %>% 
                mutate(position = as.numeric(position)) %>% 
                left_join(cleaned_gwas_dat,by = c('position' = 'base_pair_location','chromosome')) %>% 
                mutate(allele_match = case_when(ref == other_allele & alt == effect_allele ~ TRUE,
                                              ref == effect_allele & alt == other_allele ~ FALSE))  
    output
}

run_TWAS <- function(variant_df,LD_matrix) {
    subset_LD <- LD_matrix[variant_df$variant,variant_df$variant]
    denom <- t(variant_df$posterior_mean) %*% data.matrix(subset_LD) %*%  variant_df$posterior_mean
    stat <- t(variant_df$posterior_mean) %*%  variant_df$z

    zscore <- stat/sqrt(denom)
    pvalue <- pchisq(zscore * zscore,1,lower.tail = FALSE)
    output <- data.frame(zscore = zscore,pvalue = pvalue,stat = stat) 
    output
}

######### PARSE COMMAND LINE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--LDMatrix"), type="character", default=NULL,
                        help="Path to dose matrix", metavar = "type"),
    optparse::make_option(c("--PhenotypeID"), type="character", default=NULL,
                        help="Name of gene to compute LD on  ", metavar = "type"),
    optparse::make_option(c("--SummaryStats"), type="character", default=NULL,
                        help="Path to summary stats for GWAS", metavar = "type"),
    optparse::make_option(c("--SusieRes"), type="character", default=NULL,
                        help="Path to finemapping data for a gene", metavar = "type")
    #optparse::make_option(c("--AlleleFrequencies"), type="character", default=NULL,
    #                    help="Path to finemapping data for a gene", metavar = "type")
    )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
MatrixLD <- opt$LDMatrix
FineMappingRes <- opt$SusieRes
SummaryStats <- opt$SummaryStats
OutFileName <- paste0(opt$PhenotypeID,'.TWAS.txt')
PhenotypeID <- opt$PhenotypeID
#PathAlleleFrequencies <- opt$AlleleFrequencies


##################### LOAD DATA ##########################
# loading finemapping data 
susie_dat <- load_susie_data(basename(FineMappingRes)) %>% 
                mutate(variant = str_replace(variant,'chrchr','chr'))

# take variants from finemapping results and query GWAS file
# returns merged dataframe with fine-mapped variants and 
# gwas summary statistics
GWAS <- extract_gwas_data(
                    susie_dat,
                    SummaryStats
                    )

# filters gwas data based on number of measurements 
# for each variant that are present in gwas data. 
# Sometimes variants have multiple measurements in a gwas 
# but im not qutie sure why
FilteredGWAS <- GWAS %>%
    filter(allele_match == TRUE) %>% 
    mutate(z = beta/standard_error) %>% 
    group_by(variant) %>% 
    filter(dplyr::n() == 1) %>% 
    ungroup()

# reads in LD matrix and filters to the variants that are 
# present in both studies 
LD <- readRDS(MatrixLD)
#FilteredLD <- LD[FilteredGWAS$variant,FilteredGWAS$variant]


############# COMPUTE TWAS Z SCORE ########################
ResTWAS <- FilteredGWAS %>% 
            run_TWAS(LD) %>% 
            mutate(gene = PhenotypeID,
                   GWAS = file_path_sans_ext(basename(SummaryStats))
                    )
ResTWAS %>% write_tsv(OutFileName)
