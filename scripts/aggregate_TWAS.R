library(tidyverse)
library(data.table)
library(optparse)


########### COMMAND LINE ARGUMENTS ########
option_list <- list(
  optparse::make_option(c("--FilePaths"), type="character", default=NULL,
                        help="Phenotype metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
merged_tsv <- paste0(opt$OutputPrefix,'_TWAS.tsv')


list.files('.')
list.files('localized/')

############### PARSE DATA ################
filepath_df <- fread(opt$FilePaths,header = FALSE) %>% dplyr::rename('path' = 1) %>% pull(path)
number_files <- filepath_df %>% length() 
message(paste0('Number of files found: ',number_files))

TWAS_df <- dplyr::tibble(
        zscore = numeric(),
        pvalue = numeric(),
        stat = numeric(),
        gene = character(),
        GWAS = character()
)


counter <- 0
for (x in filepath_df){
    current_dat <- fread(x)
    TWAS_df <- bind_rows(TWAS_df,current_dat)
    counter <- counter + 1 
    print(counter)
}

TWAS_df %>% fwrite(merged_tsv)
