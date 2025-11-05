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



############### PARSE DATA ################
#filepath_df <- fread(opt$FilePaths,header = FALSE) %>% dplyr::rename('path' = 1) %>% pull(path)
TWAS_files <- list.files('localized/',pattern = "\\.TWAS.txt$")
number_files <- TWAS_files %>% length() 
message(paste0('Number of files found: ',number_files))

TWAS_df <- dplyr::tibble(
        zscore = numeric(),
        pvalue = numeric(),
        stat = numeric(),
        gene = character(),
        GWAS = character()
)


counter <- 0
for (x in TWAS_files){
    current_dat <- fread(paste0('localized/',x))
    TWAS_df <- bind_rows(TWAS_df,current_dat)
    counter <- counter + 1 
    print(counter)
}

TWAS_df %>% fwrite(merged_tsv,sep = '\t')
