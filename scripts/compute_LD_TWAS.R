library(tidyverse)
library(data.table)
library(arrow)
library(bedr)
library(Rfast)
library(optparse)

compute_LD <- function(X) {
  if (is.null(X)) {
    stop("X must be provided.")
  }

  # Mean impute X
  genotype_data_imputed <- apply(X, 2, function(x) {
    pos <- which(is.na(x))
    if (length(pos) != 0) {
      x[pos] <- mean(x, na.rm = TRUE)
    }
    return(x)
  })

  # Check if Rfast package is installed
  if (requireNamespace("Rfast", quietly = TRUE)) {
    # Use Rfast::cora for faster correlation calculation
    R <- Rfast::cora(genotype_data_imputed, large = TRUE)
  } else {
    # Use base R cor function if Rfast is not installed
    R <- cor(genotype_data_imputed)
  }

  colnames(R) <- rownames(R) <- colnames(genotype_data_imputed)
  R

}


########### PARSE COMMAND LINE ARGUMENTS ########
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--DoseMatrix"), type="character", default=NULL,
                        help="Path to dose matrix", metavar = "type"),
  optparse::make_option(c("--PhenotypeID"), type="character", default=NULL,
                        help="Name of gene to compute LD on  ", metavar = "type"),
  optparse::make_option(c("--VariantList"), type="character", default=NULL,
                        help=" ", metavar = "type")

    )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
DosePath <- opt$DoseMatrix
VariantListPath <- opt$VariantList
GeneID <- opt$PhenotypeID
OutFileName <- paste0(opt$PhenotypeID,'.LD.rds')


###### LOAD SUSIE DATA #########
VariantList <- fread(VariantListPath) %>%
        filter(phenotype == GeneID ) %>%  
        select(variant) %>% 
        mutate(variant = str_replace(variant,'chrchr','chr')) %>% 
        pull(variant)


########## BEGIN LD CALCULATION ########

genotype_dat <- fread(DosePath) %>% 
    mutate(variant = paste0(CHROM,'_',POS,'_',REF,'_',ALT)) %>%
    filter(variant %in% VariantList) %>%
    select(-CHROM,-POS,-REF,-ALT) %>% 
    column_to_rownames('variant') %>% 
    t() %>% 
    scale()
LD_matrix <- compute_LD(genotype_dat)
saveRDS(LD_matrix,file = OutFileName)
