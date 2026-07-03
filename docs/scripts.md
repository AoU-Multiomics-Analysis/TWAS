# Script Reference

This document describes the R scripts in `scripts/`. The WDL workflows call these scripts inside the pipeline Docker image, but they can also be run directly when the required dependencies and input files are available.

## `CreatePhenoFile.R`

Extracts expression data for a single gene from a BED-format expression matrix and prepares files for heritability estimation.

### Inputs

- `--BedFile`: BED-format gene expression file with columns `chrom`, `start`, `end`, `gene_id`, followed by per-sample expression values.
- `--PhenotypeID`: Gene ID to extract. This must match a `gene_id` value in the BED file.

### Outputs

- `gene_region.tsv`: Three-column file with `chrom`, `start`, and `end`, defining the +/- 1 Mb cis-window around the gene.
- `pheno.txt`: Phenotype file with `FID`, `IID`, and expression value, formatted for PLINK/GCTA.

## `compute_LD_TWAS.R`

Computes an LD correlation matrix for the fine-mapped variants associated with one gene.

### Inputs

- `--DoseMatrix`: Tab-separated genotype dosage matrix. Required columns are `CHROM`, `POS`, `REF`, and `ALT`, followed by per-sample dosages.
- `--PhenotypeID`: Gene or molecular trait ID whose fine-mapped variants should be used.
- `--VariantList`: File with at least `phenotype` and `variant` columns.

### Molecular Trait ID Matching

`--PhenotypeID` and `VariantList$phenotype` are normalized before matching. The script preserves the molecular trait ID but strips embedded Ensembl gene versions.

Examples:

```text
ENSG00000111647.13                                      -> ENSG00000111647
chr10:100262063:100267571:clu_12500_-:ENSG00000095485.18 -> chr10:100262063:100267571:clu_12500_-:ENSG00000095485
A0JNW5_ENSG00000111647.13                              -> A0JNW5_ENSG00000111647
```

### Output

- `<PhenotypeID>.LD.rds`: RDS file containing the LD matrix for the gene's fine-mapped variants.

## `susie_TWAS.R`

Runs the standard SuSiE-weighted TWAS analysis.

For each gene, the script extracts overlapping GWAS variants, aligns alleles, computes variant Z-scores as `BETA / SE`, and then calculates:

```text
s = t(w) %*% z
var(s) = t(w) %*% LD %*% w
z_TWAS = s / sqrt(var(s))
```

where `w` is the SuSiE posterior mean vector.

### Inputs

- `--LDMatrix`: RDS object containing either a named list of LD matrices or a single-gene LD matrix.
- `--SusieRes`: Tab-separated SuSiE fine-mapping results. Required columns include `molecular_trait_id`, `variant`, `posterior_mean`, `chromosome`, `position`, `ref`, `alt`, `pip`, and `cs_id`.
- `--SummaryStats`: Bgzipped, tabix-indexed GWAS summary statistics with `CHR`, `POS`, `REF`, `ALT`, `BETA`, `SE`, and `Pvalue`.
- `--OutputPrefix`: Prefix for the output file.
- `--GeneFilter`: Optional TSV of genes to run. When provided, the script always keeps `Coloc == TRUE`, optionally filters `chosen_label` and `type`, and restricts TWAS to genes in the `Gene` column.
- `--ChosenLabel`: Optional value to match against `GeneFilter$chosen_label`.
- `--QTLType`: Optional value to match against `GeneFilter$type`.

### Gene ID Matching

The script uses two match IDs:

- For `--GeneFilter`, it extracts the embedded Ensembl gene ID from `molecular_trait_id` and `Gene`, then strips GENCODE versions.
- For LD lookup and TWAS grouping, it preserves the molecular trait ID and strips only embedded Ensembl gene versions.

This supports eQTL IDs, splicing IDs, and pQTL IDs without collapsing multiple protein or splice traits from the same gene.

### Output

- `<OutputPrefix>.TWAS.txt`: Tab-separated TWAS output with `zscore`, `pvalue`, `stat`, `gene`, and `GWAS`.

## `conditional_rare_TWAS.R`

Runs the two-predictor rare/common TWAS follow-up. It partitions each gene's SuSiE weights into common and rare components, then reports full, common-only, rare-only, joint, and rare-conditional statistics.

The statistical details are described in [Conditional Rare TWAS](conditional-rare-twas.md).

### Inputs

- `--LDMatrix`: RDS object containing either a named list of LD matrices or a single-gene LD matrix.
- `--SusieRes`: Tab-separated SuSiE fine-mapping results. It should include the columns required by `susie_TWAS.R`, plus either a rare indicator column or `gvs_max_af`.
- `--SummaryStats`: Bgzipped, tabix-indexed GWAS summary statistics.
- `--OutputPrefix`: Prefix for the output file.
- `--RareColumn`: Rare indicator column. Defaults to `rare`.
- `--GeneFilter`: Optional TSV of genes to run. When provided, the script always keeps `Coloc == TRUE`, optionally filters `chosen_label` and `type`, and restricts TWAS to genes in the `Gene` column.
- `--ChosenLabel`: Optional value to match against `GeneFilter$chosen_label`.
- `--QTLType`: Optional value to match against `GeneFilter$type`.

### Rare Annotation

If `--RareColumn` is present, it is parsed as a rare/common indicator. Supported true values include `TRUE`, `1`, `yes`, and `rare`; supported false values include `FALSE`, `0`, `no`, and `common`.

If `--RareColumn` is absent, the script falls back to:

```text
rare = gvs_max_af < 0.01
```

Before doing this fallback, variants with missing or non-numeric `gvs_max_af` values are removed.

### Output

- `<OutputPrefix>.TWAS.two_predictor.txt`: Tab-separated output containing full, common, rare, joint, and rare-conditional TWAS statistics.

Key output columns include:

- `full_zscore`, `full_pvalue`: TWAS using all variants.
- `common_zscore`, `common_pvalue`: TWAS using common variants only.
- `rare_zscore`, `rare_pvalue`: TWAS using rare variants only.
- `chisq_joint`, `p_joint`: 2-df joint test of common and rare predictors.
- `z_rare_cond`, `p_rare_cond`: 1-df rare component test conditional on the common component.
- `Vcc`, `Vrr`, `Vcr`, `Vcor`: covariance terms for the common and rare scores.
- `n_variants`, `n_common`, `n_rare`: variant counts used for each gene.
- `status`: calculation status for the gene.

## `aggregate_TWAS.R`

Aggregates multiple per-gene or per-trait TWAS result files into one merged TSV.

### Inputs

- `--FilePaths`: Plain-text file listing individual TWAS result files, one per line.
- `--OutputPrefix`: Prefix for the merged output file.

### Output

- `<OutputPrefix>_TWAS.tsv`: Merged TWAS result file.

## `prepare_TWAS_manifest.R`

Expands QTL and GWAS manifests into a QTL x GWAS pair table for the manifest-driven WDL workflow.

### Inputs

- `--QTLManifest`: TSV with QTL inputs.
- `--GWASManifest`: TSV with GWAS inputs.
- `--OutputPrefix`: Prefix used to generate unique per-pair output names.

The QTL manifest must contain columns equivalent to:

```text
susie_path	ld_matrix_path	qtl_type
```

QTL manifest columns:

| Canonical column | Accepted aliases | Description |
|------|------|-------------|
| `susie_path` | `susie_file`, `susie`, `finemapping`, `fine_mapping`, `fine_mapping_path`, `susie_res` | Path to a SuSiE fine-mapping result file |
| `ld_matrix_path` | `ld_matrix`, `ld_path`, `ld_file`, `ldmatrix` | Path to the matching LD matrix RDS |
| `qtl_type` | `qtltype`, `type`, `label` | QTL type label, such as `eQTL`, `sQTL`, or `pQTL` |

The SuSiE file should contain the columns required by `susie_TWAS.R`. The LD
matrix should use names that match the SuSiE `molecular_trait_id` values after
the molecular-trait ID normalization described above.

The GWAS manifest must contain columns equivalent to:

```text
summary_stats_path	summary_stats_index_path	chosen_label
```

GWAS manifest columns:

| Canonical column | Accepted aliases | Description |
|------|------|-------------|
| `summary_stats_path` | `summary_stats`, `sumstats`, `sum_stats`, `sumstats_path`, `summary_stats_file` | Path to a bgzipped, tabix-indexed GWAS summary statistics file |
| `summary_stats_index_path` | `summary_stats_index`, `sumstats_index`, `sum_stats_index`, `sumstats_index_path`, `index_path`, `tbi`, `summary_stats_tbi` | Path to the tabix index for the summary statistics |
| `chosen_label` | `chosenlabel`, `trait`, `trait_label`, `label` | Trait label for the GWAS |

When a gene filter is supplied to the WDL workflow, each pair's `qtl_type` is
passed to `--QTLType` and each pair's `chosen_label` is passed to
`--ChosenLabel`. The gene filter is therefore applied separately for each QTL x
GWAS pair.

Blank rows or rows with missing required values are removed before creating the
cross-product.

### Outputs

- `twas_manifest_pairs.tsv`: Cross-product of valid QTL and GWAS manifest rows. It includes `pair_id`, row indices, raw paths, labels, safe labels, and per-pair output prefixes.
- One text file per pair-table column, used by WDL `read_lines()` to scatter over pairs.

Example pair output:

```text
pair_id	qtl_row	gwas_row	qtl_type	chosen_label	output_prefix
qtl1_gwas1	1	1	eQTL	monocyte count	TWAS.qtl1_gwas1.eqtl.monocyte_count
qtl2_gwas1	2	1	pQTL	monocyte count	TWAS.qtl2_gwas1.pqtl.monocyte_count
```

## `aggregate_manifest_TWAS.R`

Aggregates manifest-driven TWAS outputs. It combines each pair result with its metadata, then writes one global aggregate file and one aggregate file per QTL type.

### Inputs

- `--TWASFiles`: Text file listing localized TWAS result files.
- `--MetadataFiles`: Text file listing localized per-run metadata files in the same order as `--TWASFiles`.
- `--OutputPrefix`: Prefix for aggregate output files.

### Outputs

- `<OutputPrefix>.all_TWAS.tsv`: All TWAS results from all QTL x GWAS pairs.
- `<OutputPrefix>.manifest_run_metadata.tsv`: One metadata row per QTL x GWAS run.
- `by_qtl_type/<OutputPrefix>.<qtl_type>_TWAS.tsv`: One aggregate file per QTL type.
