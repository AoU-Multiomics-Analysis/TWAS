suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(arrow)
    library(bedr)
    })

######### FUNCTIONS ######################

load_susie_data <- function(path) {
susie_data <- arrow::read_parquet(path) %>% 
        select(variant,posterior_mean,chromosome,position,ref,alt,pip,cs_id)
susie_data
    
}




extract_interval <- function(finemapping_data) {
    start <- min(as.numeric(finemapping_data$position))
    end <- max(as.numeric(finemapping_data$position))
    chromosome <- unique(finemapping_data$chromosome)[1]
    interval <- paste0(chromosome,':',start,'-',end)
    interval
}
# take variants from finemapping results and query GWAS file
# returns merged dataframe with fine-mapped variants and 
# gwas summary statistics
extract_gwas_data <- function(input_susie,input_gwas) {
    require(bedr)
    header <- strsplit(readLines(input_gwas, n = 1),"\t")[[1]]
    gwas_dat <- tabix(
                    extract_interval(input_susie),
                    input_gwas,
                    check.chr = FALSE
                )
    colnames(gwas_dat) <- header
    cleaned_gwas_dat <- gwas_dat %>% 
                mutate(
                        POS = as.numeric(POS),
                        Pvalue = as.numeric(Pvalue),
                        BETA = as.numeric(BETA),
                        SE = as.numeric(SE)
                    )  
    output <- input_susie %>% 
                mutate(position = as.numeric(position)) %>% 
                left_join(cleaned_gwas_dat,by = c('position' = 'POS','chromosome' ='CHR')) %>% 
                mutate(allele_match = case_when(ref == REF & alt == ALT ~ TRUE,
                                              ref == REF & alt == ALT ~ FALSE))  
    output
}

calculate_TWAS_Z <- function(variant_df,LD_matrix) {
    components <- calculate_TWAS_components(variant_df, LD_matrix)
    components %>% select(zscore, pvalue, stat)
}

as_numeric_scalar <- function(x) {
    as.numeric(x)[1]
}

calculate_TWAS_components <- function(variant_df, LD_matrix) {
    if (nrow(variant_df) == 0) {
        stop("variant_df must contain at least one variant.")
    }

    missing_variants <- union(
        setdiff(variant_df$variant, rownames(LD_matrix)),
        setdiff(variant_df$variant, colnames(LD_matrix))
    )
    if (length(missing_variants) > 0) {
        stop(
            paste0(
                "Variants missing from LD matrix: ",
                paste(utils::head(missing_variants, 10), collapse = ", ")
            )
        )
    }

    subset_LD <- LD_matrix[variant_df$variant,variant_df$variant,drop = FALSE]
    denom <- t(variant_df$posterior_mean) %*% data.matrix(subset_LD) %*%  variant_df$posterior_mean
    stat <- t(variant_df$posterior_mean) %*%  variant_df$z

    zscore <- stat/sqrt(denom)
    pvalue <- pchisq(zscore * zscore,1,lower.tail = FALSE)
    output <- data.frame(
        zscore = as_numeric_scalar(zscore),
        pvalue = as_numeric_scalar(pvalue),
        stat = as_numeric_scalar(stat),
        denom = as_numeric_scalar(denom)
    )
    output
}

parse_rare_indicator <- function(x) {
    if (any(is.na(x))) {
        stop("Rare indicator contains missing values.")
    }

    if (is.logical(x)) {
        return(x)
    }

    if (is.numeric(x) || is.integer(x)) {
        return(x != 0)
    }

    x_clean <- tolower(trimws(as.character(x)))
    true_values <- c("true", "t", "1", "yes", "y", "rare")
    false_values <- c("false", "f", "0", "no", "n", "common")
    allowed_values <- c(true_values, false_values, NA_character_)

    unknown_values <- setdiff(unique(x_clean), allowed_values)
    if (length(unknown_values) > 0) {
        stop(
            paste0(
                "Could not parse rare indicator values: ",
                paste(unknown_values, collapse = ", ")
            )
        )
    }

    x_clean %in% true_values
}

ensure_rare_column <- function(variant_df, rare_col = "rare", af_col = "gvs_max_af", af_threshold = 0.01) {
    if (rare_col %in% colnames(variant_df)) {
        return(variant_df)
    }

    if (!af_col %in% colnames(variant_df)) {
        stop(
            paste0(
                "Missing rare indicator column: ", rare_col,
                ". Could not infer it because ", af_col, " is also missing."
            )
        )
    }

    max_af <- suppressWarnings(as.numeric(variant_df[[af_col]]))
    missing_af <- is.na(max_af)
    if (any(missing_af)) {
        message(
            paste0(
                "Filtering ",
                sum(missing_af),
                " variants with missing or non-numeric ",
                af_col,
                " before rare annotation"
            )
        )
        variant_df <- variant_df[!missing_af, ]
        max_af <- max_af[!missing_af]
    }

    if (nrow(variant_df) == 0) {
        stop(paste0("No variants remain after filtering missing or non-numeric ", af_col, " values."))
    }

    message(paste0("Rare column ", rare_col, " not found; defining it as ", af_col, " < ", af_threshold))
    variant_df[[rare_col]] <- max_af < af_threshold
    variant_df
}

is_empty_option <- function(x) {
    is.null(x) || length(x) == 0 || is.na(x) || trimws(as.character(x)) == ""
}

parse_filter_logical <- function(x, column_name) {
    if (is.logical(x)) {
        return(replace_na(x, FALSE))
    }

    if (is.numeric(x) || is.integer(x)) {
        return(!is.na(x) & x != 0)
    }

    x_clean <- tolower(trimws(as.character(x)))
    true_values <- c("true", "t", "1", "yes", "y")
    false_values <- c("false", "f", "0", "no", "n", "", "na", "nan")
    allowed_values <- c(true_values, false_values)

    unknown_values <- setdiff(unique(x_clean[!is.na(x_clean)]), allowed_values)
    if (length(unknown_values) > 0) {
        stop(
            paste0(
                "Could not parse logical values in ", column_name, ": ",
                paste(unknown_values, collapse = ", ")
            )
        )
    }

    !is.na(x_clean) & x_clean %in% true_values
}

extract_ensembl_gene_id <- function(gene_id) {
    gene_id <- as.character(gene_id)
    ensembl_pattern <- "ENS[A-Z]*G[0-9]+(?:\\.[0-9]+)?"
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
    sub("\\.[0-9]+$", "", extract_ensembl_gene_id(gene_id))
}

standardize_susie_gene_ids <- function(susie_dat, id_col = "molecular_trait_id") {
    if (!id_col %in% colnames(susie_dat)) {
        stop(paste0("Missing required SuSiE gene ID column: ", id_col))
    }

    original_ids <- as.character(susie_dat[[id_col]])
    cleaned_ids <- strip_gene_version(original_ids)
    changed <- !is.na(original_ids) & !is.na(cleaned_ids) & original_ids != cleaned_ids

    if (any(changed)) {
        example_pairs <- unique(paste0(original_ids[changed], " -> ", cleaned_ids[changed]))
        message(
            paste0(
                "Standardized ",
                dplyr::n_distinct(original_ids[changed]),
                " SuSiE molecular trait IDs to Ensembl gene IDs. Examples: ",
                paste(utils::head(example_pairs, 3), collapse = "; ")
            )
        )
    }

    susie_dat[[id_col]] <- cleaned_ids
    susie_dat
}

standardize_ld_gene_ids <- function(LD) {
    if (is_single_ld_matrix(LD) || is.null(names(LD))) {
        return(LD)
    }

    original_names <- names(LD)
    cleaned_names <- strip_gene_version(original_names)
    changed <- !is.na(original_names) & !is.na(cleaned_names) & original_names != cleaned_names

    if (anyDuplicated(cleaned_names)) {
        duplicated_names <- unique(cleaned_names[duplicated(cleaned_names)])
        stop(
            paste0(
                "Cleaning LD matrix names created duplicate gene IDs: ",
                paste(utils::head(duplicated_names, 10), collapse = ", "),
                ". Provide a single gene-level LD matrix per cleaned gene ID."
            )
        )
    }

    if (any(changed)) {
        example_pairs <- unique(paste0(original_names[changed], " -> ", cleaned_names[changed]))
        message(
            paste0(
                "Standardized ",
                sum(changed),
                " LD matrix names to Ensembl gene IDs. Examples: ",
                paste(utils::head(example_pairs, 3), collapse = "; ")
            )
        )
    }

    names(LD) <- cleaned_names
    LD
}

filter_susie_genes <- function(susie_dat, gene_filter_path = NULL, chosen_label = NULL, qtl_type = NULL) {
    if (is_empty_option(gene_filter_path)) {
        message(
            paste0(
                "No gene filter provided; keeping ",
                dplyr::n_distinct(susie_dat$molecular_trait_id),
                " genes in SuSiE data"
            )
        )
        return(susie_dat)
    }

    message(paste0("Loading gene filter: ", gene_filter_path))
    gene_filter <- fread(gene_filter_path)
    required_columns <- c("Gene", "Coloc")
    missing_columns <- setdiff(required_columns, colnames(gene_filter))
    if (length(missing_columns) > 0) {
        stop(
            paste0(
                "Gene filter is missing required columns: ",
                paste(missing_columns, collapse = ", ")
            )
        )
    }

    input_filter_genes <- gene_filter %>%
        mutate(Gene = as.character(Gene)) %>%
        filter(!is.na(Gene), Gene != "") %>%
        mutate(.gene_match_id = strip_gene_version(Gene)) %>%
        filter(!is.na(.gene_match_id), .gene_match_id != "") %>%
        distinct(.gene_match_id) %>%
        nrow()

    gene_filter <- gene_filter %>%
        mutate(.coloc_keep = parse_filter_logical(Coloc, "Coloc")) %>%
        filter(.coloc_keep)

    message(
        paste0(
            "After Coloc filter: ",
            nrow(gene_filter),
            " rows and ",
            dplyr::n_distinct(strip_gene_version(gene_filter$Gene)),
            " unique genes"
        )
    )

    if ("chosen_label" %in% colnames(gene_filter)) {
        message(
            paste0(
                "Gene filter includes ",
                dplyr::n_distinct(gene_filter$chosen_label),
                " chosen_label values after Coloc filtering"
            )
        )
    }

    if (!is_empty_option(chosen_label)) {
        if (!"chosen_label" %in% colnames(gene_filter)) {
            stop("Gene filter is missing required column for --ChosenLabel: chosen_label")
        }
        chosen_label_clean <- trimws(as.character(chosen_label))
        message(paste0("Received ChosenLabel option: '", chosen_label_clean, "'"))
        gene_filter <- gene_filter %>%
            filter(trimws(as.character(.data$chosen_label)) == .env$chosen_label_clean)
        message(
            paste0(
                "After ChosenLabel filter: ",
                nrow(gene_filter),
                " rows and ",
                dplyr::n_distinct(strip_gene_version(gene_filter$Gene)),
                " unique genes"
            )
        )
    } else if ("chosen_label" %in% colnames(gene_filter)) {
        message("No --ChosenLabel provided; keeping all chosen_label values")
    }

    if (!is_empty_option(qtl_type)) {
        if (!"type" %in% colnames(gene_filter)) {
            stop("Gene filter is missing required column for --QTLType: type")
        }
        qtl_type_clean <- trimws(as.character(qtl_type))
        message(paste0("Received QTLType option: '", qtl_type_clean, "'"))
        gene_filter <- gene_filter %>%
            filter(trimws(as.character(.data$type)) == .env$qtl_type_clean)
        message(
            paste0(
                "After QTLType filter: ",
                nrow(gene_filter),
                " rows and ",
                dplyr::n_distinct(strip_gene_version(gene_filter$Gene)),
                " unique genes"
            )
        )
    }

    keep_genes <- gene_filter %>%
        mutate(Gene = as.character(Gene)) %>%
        filter(!is.na(Gene), Gene != "") %>%
        distinct(Gene) %>%
        pull(Gene)
    keep_gene_match_ids <- strip_gene_version(keep_genes)
    keep_gene_match_ids <- unique(keep_gene_match_ids[!is.na(keep_gene_match_ids) & keep_gene_match_ids != ""])

    if (length(keep_gene_match_ids) == 0) {
        stop("No genes remain after applying gene filter.")
    }

    message(
        paste0(
            "Gene filter criteria retained ",
            length(keep_gene_match_ids),
            " of ",
            input_filter_genes,
            " unique genes in the filter file"
        )
    )

    output <- susie_dat %>%
        mutate(.gene_match_id = strip_gene_version(molecular_trait_id)) %>%
        filter(.gene_match_id %in% keep_gene_match_ids) %>%
        select(-.gene_match_id)

    if (nrow(output) == 0) {
        stop("No SuSiE variants remain after applying gene filter.")
    }

    message(
        paste0(
            "Gene filter matched ",
            dplyr::n_distinct(output$molecular_trait_id),
            " genes and ",
            nrow(output),
            " variants in SuSiE data"
        )
    )
    output
}

is_single_ld_matrix <- function(LD) {
    is.matrix(LD) || inherits(LD, "Matrix")
}

genes_with_ld <- function(LD, susie_dat) {
    if (!is_single_ld_matrix(LD)) {
        return(names(LD))
    }

    gene_ids <- unique(susie_dat$molecular_trait_id)
    if (length(gene_ids) != 1) {
        stop("LDMatrix is a single matrix, but SusieRes contains multiple genes.")
    }

    gene_ids
}

get_gene_ld <- function(LD, gene_id) {
    if (is_single_ld_matrix(LD)) {
        return(LD)
    }

    LD[[gene_id]]
}

empty_twas_row <- function(gene_id, gwas_name = NA_character_) {
  tibble(
    gene = gene_id,
    GWAS = gwas_name,
    zscore = NA_real_,
    pvalue = NA_real_
  )
}

empty_two_predictor_twas_row <- function(gene_id, gwas_name = NA_character_, status = NA_character_) {
  tibble(
    gene = gene_id,
    GWAS = gwas_name,
    n_variants = NA_integer_,
    n_common = NA_integer_,
    n_rare = NA_integer_,
    full_zscore = NA_real_,
    full_pvalue = NA_real_,
    full_stat = NA_real_,
    full_denom = NA_real_,
    common_zscore = NA_real_,
    common_pvalue = NA_real_,
    common_stat = NA_real_,
    common_denom = NA_real_,
    rare_zscore = NA_real_,
    rare_pvalue = NA_real_,
    rare_stat = NA_real_,
    rare_denom = NA_real_,
    Vcc = NA_real_,
    Vrr = NA_real_,
    Vcr = NA_real_,
    Vcor = NA_real_,
    chisq_joint = NA_real_,
    p_joint = NA_real_,
    rare_cond_stat = NA_real_,
    rare_cond_var = NA_real_,
    z_rare_cond = NA_real_,
    p_rare_cond = NA_real_,
    status = status
  )
}

calculate_two_predictor_TWAS <- function(variant_df, LD_matrix, rare_col = "rare") {
    variant_df <- ensure_rare_column(variant_df, rare_col = rare_col)

    variant_df <- variant_df %>%
        mutate(
            posterior_mean = as.numeric(posterior_mean),
            z = as.numeric(z),
            rare_indicator = parse_rare_indicator(.data[[rare_col]])
        )

    full_components <- calculate_TWAS_components(variant_df, LD_matrix)
    common_df <- variant_df %>% filter(!rare_indicator)
    rare_df <- variant_df %>% filter(rare_indicator)

    base_result <- tibble(
        n_variants = nrow(variant_df),
        n_common = nrow(common_df),
        n_rare = nrow(rare_df),
        full_zscore = full_components$zscore,
        full_pvalue = full_components$pvalue,
        full_stat = full_components$stat,
        full_denom = full_components$denom
    )

    if (nrow(common_df) == 0 || nrow(rare_df) == 0) {
        return(
            bind_cols(
                base_result,
                tibble(
                    common_zscore = NA_real_,
                    common_pvalue = NA_real_,
                    common_stat = NA_real_,
                    common_denom = NA_real_,
                    rare_zscore = NA_real_,
                    rare_pvalue = NA_real_,
                    rare_stat = NA_real_,
                    rare_denom = NA_real_,
                    Vcc = NA_real_,
                    Vrr = NA_real_,
                    Vcr = NA_real_,
                    Vcor = NA_real_,
                    chisq_joint = NA_real_,
                    p_joint = NA_real_,
                    rare_cond_stat = NA_real_,
                    rare_cond_var = NA_real_,
                    z_rare_cond = NA_real_,
                    p_rare_cond = NA_real_,
                    status = "need both common and rare variants"
                )
            )
        )
    }

    common_components <- calculate_TWAS_components(common_df, LD_matrix)
    rare_components <- calculate_TWAS_components(rare_df, LD_matrix)

    common_rare_LD <- LD_matrix[common_df$variant, rare_df$variant, drop = FALSE]
    Vcc <- common_components$denom
    Vrr <- rare_components$denom
    Vcr <- as_numeric_scalar(
        t(common_df$posterior_mean) %*%
            data.matrix(common_rare_LD) %*%
            rare_df$posterior_mean
    )

    if (!is.finite(Vcc) || !is.finite(Vrr) || Vcc <= 0 || Vrr <= 0) {
        return(
            bind_cols(
                base_result,
                tibble(
                    common_zscore = common_components$zscore,
                    common_pvalue = common_components$pvalue,
                    common_stat = common_components$stat,
                    common_denom = common_components$denom,
                    rare_zscore = rare_components$zscore,
                    rare_pvalue = rare_components$pvalue,
                    rare_stat = rare_components$stat,
                    rare_denom = rare_components$denom,
                    Vcc = Vcc,
                    Vrr = Vrr,
                    Vcr = Vcr,
                    Vcor = NA_real_,
                    chisq_joint = NA_real_,
                    p_joint = NA_real_,
                    rare_cond_stat = NA_real_,
                    rare_cond_var = NA_real_,
                    z_rare_cond = NA_real_,
                    p_rare_cond = NA_real_,
                    status = "non-positive predictor variance"
                )
            )
        )
    }

    V <- matrix(c(Vcc, Vcr, Vcr, Vrr), nrow = 2, byrow = TRUE)
    s <- c(common_components$stat, rare_components$stat)
    chisq_joint <- tryCatch(
        as_numeric_scalar(t(s) %*% solve(V) %*% s),
        error = function(e) NA_real_
    )
    p_joint <- ifelse(
        is.na(chisq_joint),
        NA_real_,
        pchisq(chisq_joint, df = 2, lower.tail = FALSE)
    )

    rare_cond_stat <- rare_components$stat - Vcr / Vcc * common_components$stat
    rare_cond_var <- Vrr - Vcr^2 / Vcc
    z_rare_cond <- ifelse(
        is.finite(rare_cond_var) && rare_cond_var > 0,
        rare_cond_stat / sqrt(rare_cond_var),
        NA_real_
    )
    p_rare_cond <- ifelse(
        is.na(z_rare_cond),
        NA_real_,
        pchisq(z_rare_cond^2, df = 1, lower.tail = FALSE)
    )

    bind_cols(
        base_result,
        tibble(
            common_zscore = common_components$zscore,
            common_pvalue = common_components$pvalue,
            common_stat = common_components$stat,
            common_denom = common_components$denom,
            rare_zscore = rare_components$zscore,
            rare_pvalue = rare_components$pvalue,
            rare_stat = rare_components$stat,
            rare_denom = rare_components$denom,
            Vcc = Vcc,
            Vrr = Vrr,
            Vcr = Vcr,
            Vcor = Vcr / sqrt(Vcc * Vrr),
            chisq_joint = chisq_joint,
            p_joint = p_joint,
            rare_cond_stat = rare_cond_stat,
            rare_cond_var = rare_cond_var,
            z_rare_cond = z_rare_cond,
            p_rare_cond = p_rare_cond,
            status = ifelse(is.na(chisq_joint) || is.na(z_rare_cond), "singular covariance", "ok")
        )
    )
}

TWAS <- function(GWAS_path,SusieData,LD) {
message('Extracting GWAS summary statistics for:')
message(GWAS_path)
PhenotypeID <- SusieData %>% distinct(geneID) %>% pull(geneID)
message(PhenotypeID)
  GWAS <- tryCatch(
    extract_gwas_data(SusieData, GWAS_path),
    error = function(e) {
      message("GWAS extraction failed for ", PhenotypeID, " : ", e$message)
      return(NULL)
    }
  )

  if (is.null(GWAS) || nrow(GWAS) == 0) {
    return(
      empty_twas_row(
        gene_id = PhenotypeID,
        gwas_name = tools::file_path_sans_ext(basename(GWAS_path))
      )
    )
  }
GeneLD <- get_gene_ld(LD, PhenotypeID)
  if (is.null(GeneLD)) {
    return(
      empty_twas_row(
        gene_id = PhenotypeID,
        gwas_name = tools::file_path_sans_ext(basename(GWAS_path))
      )
    )
  }
    
# filters gwas data based on number of measurements 
# for each variant that are present in gwas data. 
# Sometimes variants have multiple measurements in a gwas 
# but im not qutie sure why
FilteredGWAS <- GWAS %>%
    filter(allele_match == TRUE) %>% 
    mutate(z = BETA/SE) %>% 
    group_by(variant) %>% 
    filter(dplyr::n() == 1) %>% 
    ungroup()

  if (nrow(FilteredGWAS) == 0) {
    return(
      empty_twas_row(
        gene_id = PhenotypeID,
        gwas_name = tools::file_path_sans_ext(basename(GWAS_path))
      )
    )
  }
#message('Computing TWAS Z')
ResTWAS <- FilteredGWAS %>% 
            calculate_TWAS_Z(GeneLD) %>% 
            mutate(gene = PhenotypeID,
                   GWAS = tools::file_path_sans_ext(basename(GWAS_path))
                    )

ResTWAS
}

TWAS_two_predictor <- function(GWAS_path, SusieData, LD, rare_col = "rare") {
message('Extracting GWAS summary statistics for:')
message(GWAS_path)
PhenotypeID <- SusieData %>% distinct(geneID) %>% pull(geneID)
message(PhenotypeID)
  GWAS <- tryCatch(
    extract_gwas_data(SusieData, GWAS_path),
    error = function(e) {
      message("GWAS extraction failed for ", PhenotypeID, " : ", e$message)
      return(NULL)
    }
  )

  gwas_name <- tools::file_path_sans_ext(basename(GWAS_path))
  if (is.null(GWAS) || nrow(GWAS) == 0) {
    return(
      empty_two_predictor_twas_row(
        gene_id = PhenotypeID,
        gwas_name = gwas_name,
        status = "no GWAS variants"
      )
    )
  }
GeneLD <- get_gene_ld(LD, PhenotypeID)
  if (is.null(GeneLD)) {
    return(
      empty_two_predictor_twas_row(
        gene_id = PhenotypeID,
        gwas_name = gwas_name,
        status = "no LD matrix for gene"
      )
    )
  }
    
FilteredGWAS <- GWAS %>%
    filter(allele_match == TRUE) %>% 
    mutate(z = BETA/SE) %>% 
    group_by(variant) %>% 
    filter(dplyr::n() == 1) %>% 
    ungroup()

  if (nrow(FilteredGWAS) == 0) {
    return(
      empty_two_predictor_twas_row(
        gene_id = PhenotypeID,
        gwas_name = gwas_name,
        status = "no uniquely matched GWAS variants"
      )
    )
  }

ResTWAS <- tryCatch(
    calculate_two_predictor_TWAS(FilteredGWAS, GeneLD, rare_col = rare_col) %>% 
        mutate(
            gene = PhenotypeID,
            GWAS = gwas_name,
            .before = 1
        ),
    error = function(e) {
        message("Two-predictor TWAS failed for ", PhenotypeID, " : ", e$message)
        empty_two_predictor_twas_row(
            gene_id = PhenotypeID,
            gwas_name = gwas_name,
            status = e$message
        )
    }
)

ResTWAS
}

run_susie_TWAS <- function() {
    ######### PARSE COMMAND LINE ARGUMENTS #########
    option_list <- list(
        optparse::make_option(c("--LDMatrix"), type="character", default=NULL,
                            help="Path to RDS object that is a list and contains a matrix for each gene in analysis", metavar = "type"),
        optparse::make_option(c("--SummaryStats"), type="character", default=NULL,
                            help="path to indexed summary stats", metavar = "type"),
        optparse::make_option(c("--SusieRes"), type="character", default=NULL,
                            help="Path to finemapping data for a gene", metavar = "type"),
        optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                            help="Path to finemapping data for a gene", metavar = "type"),
        optparse::make_option(c("--GeneFilter"), type="character", default=NULL,
                            help="Optional path to gene filter TSV. Always filters Coloc == TRUE when provided. Embedded Ensembl gene IDs are extracted and trailing GENCODE version suffixes are ignored for gene matching.", metavar = "type"),
        optparse::make_option(c("--ChosenLabel"), type="character", default=NULL,
                            help="Optional chosen_label value to filter GeneFilter.", metavar = "type"),
        optparse::make_option(c("--QTLType"), type="character", default=NULL,
                            help="Optional type value to filter GeneFilter.", metavar = "type")

        )

    opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
    MatrixLD <- opt$LDMatrix
    FineMappingRes <- opt$SusieRes
    SummaryStats <- opt$SummaryStats
    GeneFilter <- opt$GeneFilter
    ChosenLabel <- opt$ChosenLabel
    QTLType <- opt$QTLType
    OutFileName <- paste0(opt$OutputPrefix,'.TWAS.txt')


    ##################### LOAD DATA ##########################
    # loading finemapping data 
    message('Loading fine mapping')
    susie_dat <- fread(FineMappingRes) %>% 
                    mutate(variant = str_replace(variant,'chrchr','chr')) %>% 
                    mutate(chromosome = str_remove_all(chromosome,'chr')) %>% 
                    mutate(chromosome = paste0('chr',chromosome)) %>%
                    standardize_susie_gene_ids()
    susie_dat <- filter_susie_genes(
        susie_dat,
        gene_filter_path = GeneFilter,
        chosen_label = ChosenLabel,
        qtl_type = QTLType
    )



    # Loads LD matrix
    message('Loading LD matrix')
    LD <- readRDS(MatrixLD) %>%
        standardize_ld_gene_ids()

    ############# COMPUTE TWAS Z SCORE ########################
    # computes TWAS Z statistic 
    #ResTWAS <-  susie_dat %>%
                #TWAS(SummaryStats,LD,PhenotypeID)

    # loop over summary stats and perform TWAS 
    # for each phenotype 
    #ResTWAS <-  SummaryStatsList %>%
                    #map_dfr(~TWAS(.x, SusieData = susie_dat,LD = LD,PhenotypeID = PhenotypeID))

    ResTWAS <- susie_dat %>% 
        filter(molecular_trait_id %in% genes_with_ld(LD, susie_dat)) %>%    
        mutate(geneID = molecular_trait_id) %>% 
        group_by(molecular_trait_id) %>% 
        group_modify(~TWAS(SummaryStats,.,LD))
    ResTWAS %>% write_tsv(OutFileName)
}

if (sys.nframe() == 0) {
    run_susie_TWAS()
}
