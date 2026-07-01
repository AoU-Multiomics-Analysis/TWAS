# Workflow Reference

This document describes the WDL workflows in `workflows/`.

## `ComputeHeritability.wdl`

Workflow: `EstimateHeritability`

Estimates cis-SNP heritability for one gene using PLINK2 and GCTA REML. The workflow extracts expression for the target gene, defines a cis-window, subsets genotype data, builds a genetic relationship matrix, and estimates heritability.

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `pvar` | File | PLINK2 `.pvar` variant information file |
| `pgen` | File | PLINK2 `.pgen` genotype file |
| `psam` | File | PLINK2 `.psam` sample information file |
| `ExpressionBed` | File | BED-format gene expression matrix |
| `Covars` | File | Quantitative covariates file for GCTA REML |
| `PhenotypeID` | String | Gene ID to compute heritability for |
| `Memory` | Int | Memory in GB |
| `NumPrempt` | Int | Number of preemptible retries |

### Outputs

| Name | Type | Description |
|------|------|-------------|
| `HeritabilityEstimate` | File | GCTA `.hsq` file containing the heritability estimate |

## `compute_LD_TWAS.wdl`

Workflow: `PreprocessLD`

Computes the LD matrix for the fine-mapped variants associated with one gene. This is a preprocessing step for standard and conditional TWAS.

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `DoseMatrix` | File | Tab-separated genotype dosage matrix |
| `VariantList` | File | File with `phenotype` and `variant` columns; embedded Ensembl gene versions are stripped while preserving molecular trait IDs |
| `PhenotypeID` | String | Gene or molecular trait ID; embedded Ensembl gene versions are stripped while preserving molecular trait IDs |
| `Memory` | Int | Memory in GB |
| `NumPrempt` | Int | Number of preemptible retries |

### Outputs

| Name | Type | Description |
|------|------|-------------|
| `MatrixLD` | File | RDS file containing the LD matrix |

## `susie_TWAS.wdl`

Workflow: `TWAS`

Runs standard SuSiE-weighted TWAS for one GWAS trait using the LD matrix, SuSiE fine-mapping results, and GWAS summary statistics.

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `LDMatrix` | File | RDS LD matrix file |
| `SumStats` | File | Bgzipped GWAS summary statistics |
| `SumStatsIndex` | File | Tabix index for `SumStats` |
| `FineMapping` | File | SuSiE fine-mapping result file |
| `NameGWAS` | String | Output prefix for the GWAS |
| `GeneFilter` | File? | Optional gene filter TSV |
| `ChosenLabel` | String | Optional `chosen_label` filter; blank means no filter |
| `QTLType` | String | Optional `type` filter; blank means no filter |
| `Memory` | Int | Memory in GB |
| `NumPrempt` | Int | Number of preemptible retries |

### Outputs

| Name | Type | Description |
|------|------|-------------|
| `ResTWAS` | File | Standard TWAS result file, `<NameGWAS>.TWAS.txt` |

## `conditional_rare_TWAS.wdl`

Workflow: `ConditionalRareTWASAnalysis`

Runs the two-predictor rare/common TWAS analysis. The workflow uses the same inputs as standard TWAS, with an additional `RareColumn` input to define rare variants when the fine-mapping file already contains a rare annotation.

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `LDMatrix` | File | RDS LD matrix file or named list of gene LD matrices |
| `SumStats` | File | Bgzipped GWAS summary statistics |
| `SumStatsIndex` | File | Tabix index for `SumStats` |
| `FineMapping` | File | SuSiE fine-mapping result file |
| `NameGWAS` | String | Output prefix for the GWAS |
| `RareColumn` | String | Rare indicator column in `FineMapping`; defaults to `rare` |
| `GeneFilter` | File? | Optional gene filter TSV |
| `ChosenLabel` | String | Optional `chosen_label` filter; blank means no filter |
| `QTLType` | String | Optional `type` filter; blank means no filter |
| `Memory` | Int | Memory in GB |
| `NumPrempt` | Int | Number of preemptible retries |

### Outputs

| Name | Type | Description |
|------|------|-------------|
| `ResTWAS` | File | Conditional rare TWAS result file, `<NameGWAS>.TWAS.two_predictor.txt` |

## `aggregate_TWAS.wdl`

Workflow: `AggregateTWASWorkflow`

Aggregates multiple TWAS result files into one TSV. The workflow downloads result files from Google Cloud Storage with `gsutil`, then calls `aggregate_TWAS.R`.

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `TWASResFOFN` | File | File of file names listing GCS paths to TWAS result files |
| `OutputPrefix` | String | Prefix for the merged output |
| `Memory` | Int | Memory in GB |
| `NumThreads` | Int | Number of CPU threads |

### Outputs

| Name | Type | Description |
|------|------|-------------|
| `AggregatedTWAS` | File | Merged TWAS result file, `<OutputPrefix>_TWAS.tsv` |
| `filelist` | File | Localized file list used during aggregation |

## Runtime Image

The workflows use:

```text
ghcr.io/aou-multiomics-analysis/twas:main
```

The Docker image is built from `envs/Dockerfile` and published by the GitHub Actions workflow in `.github/workflows/docker-image.yml`.
