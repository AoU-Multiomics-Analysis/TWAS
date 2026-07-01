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

extract_ensembl_gene_id <- function(gene_id) {
  gene_id <- as.character(gene_id)
  ensembl_pattern <- "ENS[A-Z]*G[0-9]+([.][0-9]+)?"
  matches <- stringr::str_extract_all(gene_id, ensembl_pattern)

  extracted <- vapply(seq_along(gene_id), function(i) {
    if (is.na(gene_id[i])) {
      return(NA_character_)
    }
    if (length(matches[[i]]) > 0 && !all(is.na(matches[[i]]))) {
      return(tail(stats::na.omit(matches[[i]]), 1))
    }
    gene_id[i]
  }, character(1))

  extracted
}

strip_gene_version <- function(gene_id) {
  sub("[.][0-9]+$", "", extract_ensembl_gene_id(gene_id))
}

normalize_molecular_trait_id <- function(molecular_trait_id) {
  molecular_trait_id <- as.character(molecular_trait_id)
  stringr::str_replace_all(
    molecular_trait_id,
    "(ENS[A-Z]*G[0-9]+)[.][0-9]+",
    "\\1"
  )
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
GeneID <- normalize_molecular_trait_id(PhenotypeID)
if (!is.na(GeneID) && PhenotypeID != GeneID) {
  message(paste0("Standardized PhenotypeID for matching: ", PhenotypeID, " -> ", GeneID))
}
OutFileName <- paste0(opt$PhenotypeID,'.LD.rds')


###### LOAD SUSIE DATA #########
VariantListDat <- fread(VariantListPath)
variant_list_standardized_count <- sum(
  !is.na(VariantListDat$phenotype) &
    normalize_molecular_trait_id(VariantListDat$phenotype) != as.character(VariantListDat$phenotype),
  na.rm = TRUE
)
if (variant_list_standardized_count > 0) {
  message(
    paste0(
      "Standardized ",
      variant_list_standardized_count,
      " VariantList phenotype rows by stripping embedded Ensembl gene versions for matching"
    )
  )
} else {
  message("No VariantList phenotype rows required molecular trait ID standardization")
}

VariantList <- VariantListDat %>%
        mutate(phenotype = normalize_molecular_trait_id(phenotype)) %>%
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
