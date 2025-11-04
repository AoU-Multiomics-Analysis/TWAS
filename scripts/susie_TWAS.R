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

calculate_TWAS_Z <- function(variant_df,LD_matrix) {
    subset_LD <- LD_matrix[variant_df$variant,variant_df$variant]
    denom <- t(variant_df$posterior_mean) %*% data.matrix(subset_LD) %*%  variant_df$posterior_mean
    stat <- t(variant_df$posterior_mean) %*%  variant_df$z

    zscore <- stat/sqrt(denom)
    pvalue <- pchisq(zscore * zscore,1,lower.tail = FALSE)
    output <- data.frame(zscore = zscore,pvalue = pvalue,stat = stat) 
    output
}

TWAS <- function(GWAS_path,SusieData,LD,PhenotypeID) {
GWAS <- extract_gwas_data(
                    SusieData,
                    GWAS_path 
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

ResTWAS <- FilteredGWAS %>% 
            calculate_TWAS_Z(LD) %>% 
            mutate(gene = PhenotypeID,
                   GWAS = tools::file_path_sans_ext(basename(SummaryStats))
                    )
ResTWAS
}


######### PARSE COMMAND LINE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--LDMatrix"), type="character", default=NULL,
                        help="Path to dose matrix", metavar = "type"),
    optparse::make_option(c("--PhenotypeID"), type="character", default=NULL,
                        help="Name of gene to compute LD on  ", metavar = "type"),
    optparse::make_option(c("--SummaryStats"), type="character", default=NULL,
                        help="comma seperated list of Sumstats to use", metavar = "type"),
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

# convert comma seperated list of summary stat files into a dataframe 
# to loop over and calculate TWAS Z scores 
SummaryStatsList <- strsplit(SummaryStats,',')

##################### LOAD DATA ##########################
# loading finemapping data 
susie_dat <- load_susie_data(FineMappingRes) %>% 
                mutate(variant = str_replace(variant,'chrchr','chr'))

# Loads LD matrix
LD <- readRDS(MatrixLD)

############# COMPUTE TWAS Z SCORE ########################
# computes TWAS Z statistic 
#ResTWAS <-  susie_dat %>%
            #TWAS(SummaryStats,LD,PhenotypeID)

# loop over summary stats and perform TWAS 
# for each phenotype 
ResTWAS <-  SummaryStatsList %>%
                map_dfr(~TWAS(.x, SusieData = susie_dat,LD = LD,PhenotypeID = PhenotypeID))

# write to output
ResTWAS %>% write_tsv(OutFileName)
