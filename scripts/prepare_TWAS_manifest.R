suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(optparse)
    library(stringr)
})

option_list <- list(
    optparse::make_option(c("--QTLManifest"), type = "character", default = NULL,
                          help = "TSV with QTL inputs. Required columns: susie_path, ld_matrix_path, qtl_type."),
    optparse::make_option(c("--GWASManifest"), type = "character", default = NULL,
                          help = "TSV with GWAS inputs. Required columns: summary_stats_path, summary_stats_index_path, chosen_label."),
    optparse::make_option(c("--OutputPrefix"), type = "character", default = "TWAS",
                          help = "Prefix used for generated pair output names.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

required_option <- function(x, name) {
    if (is.null(x) || is.na(x) || trimws(as.character(x)) == "") {
        stop(paste0("Missing required option --", name), call. = FALSE)
    }
    x
}

normalize_colnames <- function(x) {
    tolower(gsub("[^a-zA-Z0-9]+", "_", trimws(x)))
}

find_column <- function(dat, aliases, manifest_name) {
    normalized_names <- normalize_colnames(names(dat))
    normalized_aliases <- normalize_colnames(aliases)
    idx <- match(normalized_aliases, normalized_names, nomatch = 0)
    idx <- idx[idx > 0]

    if (length(idx) == 0) {
        stop(
            paste0(
                manifest_name,
                " is missing required column. Accepted names: ",
                paste(aliases, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    names(dat)[idx[1]]
}

clean_value <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[x == ""] <- NA_character_
    x
}

safe_label <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(is.na(x) | x == "", "missing", x)
}

write_column <- function(dat, column, path) {
    writeLines(as.character(dat[[column]]), path, useBytes = TRUE)
}

qtl_manifest <- required_option(opt$QTLManifest, "QTLManifest")
gwas_manifest <- required_option(opt$GWASManifest, "GWASManifest")
output_prefix <- required_option(opt$OutputPrefix, "OutputPrefix")

qtl_dat <- fread(qtl_manifest)
gwas_dat <- fread(gwas_manifest)

qtl_cols <- list(
    susie_path = find_column(qtl_dat, c("susie_path", "susie_file", "susie", "finemapping", "fine_mapping", "fine_mapping_path", "susie_res"), "QTLManifest"),
    ld_matrix_path = find_column(qtl_dat, c("ld_matrix_path", "ld_matrix", "ld_path", "ld_file", "ldmatrix"), "QTLManifest"),
    qtl_type = find_column(qtl_dat, c("qtl_type", "qtltype", "type", "label"), "QTLManifest")
)

gwas_cols <- list(
    summary_stats_path = find_column(gwas_dat, c("summary_stats_path", "summary_stats", "sumstats", "sum_stats", "sumstats_path", "summary_stats_file"), "GWASManifest"),
    summary_stats_index_path = find_column(gwas_dat, c("summary_stats_index_path", "summary_stats_index", "sumstats_index", "sum_stats_index", "sumstats_index_path", "index_path", "tbi", "summary_stats_tbi"), "GWASManifest"),
    chosen_label = find_column(gwas_dat, c("chosen_label", "chosenlabel", "trait", "trait_label", "label"), "GWASManifest")
)

qtl_clean <- qtl_dat %>%
    transmute(
        qtl_row = dplyr::row_number(),
        susie_path = clean_value(.data[[qtl_cols$susie_path]]),
        ld_matrix_path = clean_value(.data[[qtl_cols$ld_matrix_path]]),
        qtl_type = clean_value(.data[[qtl_cols$qtl_type]])
    ) %>%
    filter(!is.na(susie_path), !is.na(ld_matrix_path), !is.na(qtl_type))

gwas_clean <- gwas_dat %>%
    transmute(
        gwas_row = dplyr::row_number(),
        summary_stats_path = clean_value(.data[[gwas_cols$summary_stats_path]]),
        summary_stats_index_path = clean_value(.data[[gwas_cols$summary_stats_index_path]]),
        chosen_label = clean_value(.data[[gwas_cols$chosen_label]])
    ) %>%
    filter(!is.na(summary_stats_path), !is.na(summary_stats_index_path), !is.na(chosen_label))

if (nrow(qtl_clean) == 0) {
    stop("No valid QTL manifest rows remain after filtering missing required values.", call. = FALSE)
}

if (nrow(gwas_clean) == 0) {
    stop("No valid GWAS manifest rows remain after filtering missing required values.", call. = FALSE)
}

pairs <- tidyr::crossing(qtl_clean, gwas_clean) %>%
    mutate(
        qtl_type_safe = safe_label(qtl_type),
        chosen_label_safe = safe_label(chosen_label),
        pair_id = paste0("qtl", qtl_row, "_gwas", gwas_row),
        output_prefix = paste(output_prefix, pair_id, qtl_type_safe, chosen_label_safe, sep = ".")
    ) %>%
    select(
        pair_id,
        qtl_row,
        gwas_row,
        susie_path,
        ld_matrix_path,
        qtl_type,
        qtl_type_safe,
        summary_stats_path,
        summary_stats_index_path,
        chosen_label,
        chosen_label_safe,
        output_prefix
    )

fwrite(pairs, "twas_manifest_pairs.tsv", sep = "\t")

write_column(pairs, "pair_id", "pair_id.txt")
write_column(pairs, "susie_path", "susie_path.txt")
write_column(pairs, "ld_matrix_path", "ld_matrix_path.txt")
write_column(pairs, "qtl_type", "qtl_type.txt")
write_column(pairs, "qtl_type_safe", "qtl_type_safe.txt")
write_column(pairs, "summary_stats_path", "summary_stats_path.txt")
write_column(pairs, "summary_stats_index_path", "summary_stats_index_path.txt")
write_column(pairs, "chosen_label", "chosen_label.txt")
write_column(pairs, "chosen_label_safe", "chosen_label_safe.txt")
write_column(pairs, "output_prefix", "output_prefix.txt")

message(
    paste0(
        "Prepared ",
        nrow(pairs),
        " TWAS manifest pairs from ",
        nrow(qtl_clean),
        " QTL rows and ",
        nrow(gwas_clean),
        " GWAS rows."
    )
)
