suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(optparse)
})

option_list <- list(
    optparse::make_option(c("--TWASFiles"), type = "character", default = NULL,
                          help = "Text file listing localized TWAS result files."),
    optparse::make_option(c("--MetadataFiles"), type = "character", default = NULL,
                          help = "Text file listing localized metadata files, in the same order as TWASFiles."),
    optparse::make_option(c("--OutputPrefix"), type = "character", default = "TWAS",
                          help = "Prefix for aggregate outputs.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

required_option <- function(x, name) {
    if (is.null(x) || is.na(x) || trimws(as.character(x)) == "") {
        stop(paste0("Missing required option --", name), call. = FALSE)
    }
    x
}

safe_label <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(is.na(x) | x == "", "missing", x)
}

twas_file_list <- required_option(opt$TWASFiles, "TWASFiles")
metadata_file_list <- required_option(opt$MetadataFiles, "MetadataFiles")
output_prefix <- required_option(opt$OutputPrefix, "OutputPrefix")

twas_files <- readLines(twas_file_list, warn = FALSE)
metadata_files <- readLines(metadata_file_list, warn = FALSE)
twas_files <- twas_files[twas_files != ""]
metadata_files <- metadata_files[metadata_files != ""]

if (length(twas_files) == 0) {
    stop("No TWAS files provided for aggregation.", call. = FALSE)
}

if (length(twas_files) != length(metadata_files)) {
    stop("TWASFiles and MetadataFiles must contain the same number of entries.", call. = FALSE)
}

all_results <- vector("list", length(twas_files))
all_metadata <- vector("list", length(twas_files))

for (i in seq_along(twas_files)) {
    current_twas <- fread(twas_files[[i]])
    current_metadata <- fread(metadata_files[[i]])

    if (nrow(current_metadata) != 1) {
        stop(paste0("Metadata file must contain exactly one data row: ", metadata_files[[i]]), call. = FALSE)
    }

    metadata_cols <- setdiff(names(current_metadata), names(current_twas))
    current_twas <- bind_cols(
        current_twas,
        current_metadata[rep(1, nrow(current_twas)), metadata_cols, with = FALSE]
    )

    all_results[[i]] <- current_twas
    all_metadata[[i]] <- current_metadata
}

merged_results <- bind_rows(all_results)
merged_metadata <- bind_rows(all_metadata)

all_output <- paste0(output_prefix, ".all_TWAS.tsv")
metadata_output <- paste0(output_prefix, ".manifest_run_metadata.tsv")
fwrite(merged_results, all_output, sep = "\t")
fwrite(merged_metadata, metadata_output, sep = "\t")

dir.create("by_qtl_type", showWarnings = FALSE)
qtl_types <- sort(unique(merged_results$qtl_type))

for (qtl_type in qtl_types) {
    qtl_type_safe <- unique(merged_results$qtl_type_safe[merged_results$qtl_type == qtl_type])
    if (length(qtl_type_safe) != 1) {
        qtl_type_safe <- safe_label(qtl_type)[[1]]
    }

    output_path <- file.path("by_qtl_type", paste0(output_prefix, ".", qtl_type_safe, "_TWAS.tsv"))
    qtl_results <- merged_results %>%
        filter(.data$qtl_type == .env$qtl_type)
    fwrite(qtl_results, output_path, sep = "\t")
}

message(
    paste0(
        "Aggregated ",
        length(twas_files),
        " TWAS files into ",
        all_output,
        " and ",
        length(qtl_types),
        " QTL-type files."
    )
)
