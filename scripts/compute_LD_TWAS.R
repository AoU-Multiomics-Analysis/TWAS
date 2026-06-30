library(tidyverse)
library(data.table)
library(arrow)
library(bedr)
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

has_gene_version <- function(gene_id) {
  grepl("[.][0-9]+$", as.character(gene_id))
}

strip_gene_version <- function(gene_id) {
  gene_id <- as.character(gene_id)
  ifelse(has_gene_version(gene_id), sub("[.][0-9]+$", "", gene_id), gene_id)
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
PhenotypeID <- opt$PhenotypeID
if (has_gene_version(PhenotypeID)) {
  message(paste0("Detected gene version suffix in PhenotypeID: ", PhenotypeID))
} else {
  message(paste0("No gene version suffix detected in PhenotypeID: ", PhenotypeID))
}
GeneID <- strip_gene_version(PhenotypeID)
OutFileName <- paste0(opt$PhenotypeID,'.LD.rds')


###### LOAD SUSIE DATA #########
VariantListDat <- fread(VariantListPath)
variant_list_version_count <- sum(has_gene_version(VariantListDat$phenotype), na.rm = TRUE)
if (variant_list_version_count > 0) {
  message(
    paste0(
      "Detected gene version suffixes in ",
      variant_list_version_count,
      " VariantList phenotype rows; stripping for matching"
    )
  )
} else {
  message("No gene version suffixes detected in VariantList phenotype column")
}

VariantList <- VariantListDat %>%
        mutate(phenotype = strip_gene_version(phenotype)) %>%
        filter(phenotype == GeneID ) %>%
        select(variant) %>% 
        mutate(variant = str_replace(variant,'chrchr','chr')) %>% 
        pull(variant)

message(paste0("Matched ", length(VariantList), " variants for PhenotypeID after version stripping: ", GeneID))
if (length(VariantList) == 0) {
  stop(paste0("No variants found for PhenotypeID after version stripping: ", GeneID))
}


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
