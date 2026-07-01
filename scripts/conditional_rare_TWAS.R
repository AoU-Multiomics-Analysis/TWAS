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
                          help = "Path to SuSiE fine-mapping results with a rare indicator column or gvs_max_af."),
    optparse::make_option(c("--OutputPrefix"), type = "character", default = NULL,
                          help = "Output prefix."),
    optparse::make_option(c("--RareColumn"), type = "character", default = "rare",
                          help = "Column in SusieRes marking rare variants. If missing, infer from gvs_max_af < 0.01. Default: rare"),
    optparse::make_option(c("--GeneFilter"), type = "character", default = NULL,
                          help = "Optional path to gene filter TSV. Always filters Coloc == TRUE when provided. Embedded Ensembl gene IDs are extracted for gene filtering; molecular trait IDs are preserved for LD matching."),
    optparse::make_option(c("--ChosenLabel"), type = "character", default = NULL,
                          help = "Optional chosen_label value to filter GeneFilter."),
    optparse::make_option(c("--QTLType"), type = "character", default = NULL,
                          help = "Optional type value to filter GeneFilter.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

MatrixLD <- required_option(opt, "LDMatrix")
FineMappingRes <- required_option(opt, "SusieRes")
SummaryStats <- required_option(opt, "SummaryStats")
OutputPrefix <- required_option(opt, "OutputPrefix")
RareColumn <- opt$RareColumn
GeneFilter <- opt$GeneFilter
ChosenLabel <- opt$ChosenLabel
QTLType <- opt$QTLType
OutFileName <- paste0(OutputPrefix, ".TWAS.two_predictor.txt")

message("Loading fine mapping")
susie_dat <- fread(FineMappingRes) %>%
    mutate(variant = str_replace(variant, "chrchr", "chr")) %>%
    mutate(chromosome = str_remove_all(chromosome, "chr")) %>%
    mutate(chromosome = paste0("chr", chromosome)) %>%
    standardize_susie_gene_ids()

susie_dat <- ensure_rare_column(susie_dat, rare_col = RareColumn)
susie_dat <- filter_susie_genes(
    susie_dat,
    gene_filter_path = GeneFilter,
    chosen_label = ChosenLabel,
    qtl_type = QTLType
)

message("Loading LD matrix")
LD <- readRDS(MatrixLD) %>%
    standardize_ld_gene_ids()
gene_ids <- genes_with_ld(LD, susie_dat)

message("Computing two-predictor TWAS statistics")
ResTWAS <- susie_dat %>%
    filter(molecular_trait_id %in% gene_ids) %>%
    mutate(geneID = molecular_trait_id) %>%
    group_by(molecular_trait_id) %>%
    group_modify(~TWAS_two_predictor(SummaryStats, ., LD, rare_col = RareColumn)) %>%
    ungroup()

output_twas_genes <- dplyr::n_distinct(ResTWAS$gene)
output_twas_rows <- nrow(ResTWAS)

ResTWAS %>% write_tsv(OutFileName)

message(
    paste0(
        "Wrote output TWAS file with ",
        output_twas_genes,
        " genes and ",
        output_twas_rows,
        " rows: ",
        OutFileName
    )
)
