library(tidyverse)
library(data.table) 
library(optparse)

########### PARSE COMMAND LINE ARGUMENTS ########
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--BedFile"), type="character", default=NULL,
                        help="Path ot bed file", metavar = "type"),
  optparse::make_option(c("--PhenotypeID"), type="character", default=NULL,
                        help="Name of gene extract  ", metavar = "type")
    )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
BedFilePath<- opt$BedFile
PhenotypeID <- opt$PhenotypeID


############ EXTRACT GENE VALUE ##########
SubsetBed <- fread(BedFilePath) %>% 
            dplyr::filter(gene_id == PhenotypeID)
GeneMeta <- SubsetBed %>%
                dplyr::rename('chrom' = 1) %>% 
                mutate(chrom = str_remove(chrom,'chr')) %>% 
                select(1,2,3)


PhenoDf <- SubsetBed %>% 
            select(-1,-2,-3) %>%
            pivot_longer(!gene_id) %>% 
            dplyr::select(-1) %>% 
            mutate(FID =0,IID = name) %>% 
            select(FID,IID,value)

GeneMeta %>% write_tsv('gene_region.tsv')
PhenoDf %>% write_tsv('pheno.txt')

