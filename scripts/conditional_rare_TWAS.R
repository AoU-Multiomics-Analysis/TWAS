suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
})

get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)

    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }

    getwd()
}

helper_script <- file.path(get_script_dir(), "susie_TWAS.R")
if (!file.exists(helper_script)) {
    helper_script <- "susie_TWAS.R"
}
source(helper_script)

required_option <- function(opt, name) {
    value <- opt[[name]]
    if (is.null(value) || is.na(value) || value == "") {
        stop(paste0("Missing required option --", name), call. = FALSE)
    }
    value
}

option_list <- list(
    optparse::make_option(c("--LDMatrix"), type = "character", default = NULL,
                          help = "Path to RDS LD object. May be a named list of gene LD matrices or a single-gene LD matrix."),
    optparse::make_option(c("--SummaryStats"), type = "character", default = NULL,
                          help = "Path to bgzipped, tabix-indexed GWAS summary statistics."),
    optparse::make_option(c("--SusieRes"), type = "character", default = NULL,
                          help = "Path to SuSiE fine-mapping results with a rare indicator column."),
    optparse::make_option(c("--OutputPrefix"), type = "character", default = NULL,
                          help = "Output prefix."),
    optparse::make_option(c("--RareColumn"), type = "character", default = "rare",
                          help = "Column in SusieRes marking rare variants. Default: rare")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

MatrixLD <- required_option(opt, "LDMatrix")
FineMappingRes <- required_option(opt, "SusieRes")
SummaryStats <- required_option(opt, "SummaryStats")
OutputPrefix <- required_option(opt, "OutputPrefix")
RareColumn <- opt$RareColumn
OutFileName <- paste0(OutputPrefix, ".TWAS.two_predictor.txt")

message("Loading fine mapping")
susie_dat <- fread(FineMappingRes) %>%
    mutate(variant = str_replace(variant, "chrchr", "chr")) %>%
    mutate(chromosome = str_remove_all(chromosome, "chr")) %>%
    mutate(chromosome = paste0("chr", chromosome))

if (!RareColumn %in% colnames(susie_dat)) {
    stop(paste0("Fine-mapping results are missing rare column: ", RareColumn), call. = FALSE)
}

message("Loading LD matrix")
LD <- readRDS(MatrixLD)
gene_ids <- genes_with_ld(LD, susie_dat)

message("Computing two-predictor TWAS statistics")
ResTWAS <- susie_dat %>%
    filter(molecular_trait_id %in% gene_ids) %>%
    mutate(geneID = molecular_trait_id) %>%
    group_by(molecular_trait_id) %>%
    group_modify(~TWAS_two_predictor(SummaryStats, ., LD, rare_col = RareColumn)) %>%
    ungroup()

ResTWAS %>% write_tsv(OutFileName)
