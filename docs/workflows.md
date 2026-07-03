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

## `manifest_TWAS.wdl`

Workflow: `ManifestTWAS`

Runs TWAS over the cross-product of a QTL manifest and a GWAS manifest. This is the preferred workflow for running multiple QTL data types against multiple GWAS traits in one submission.

The workflow has three stages:

1. Expand the QTL and GWAS manifests into one row per QTL x GWAS pair.
2. Localize each pair's SuSiE file, LD matrix, GWAS summary statistics, and summary-stat index with `gsutil cp` when paths begin with `gs://`.
3. Aggregate all pair outputs globally and separately by QTL type.

### QTL Manifest

The QTL manifest is a tab-separated file with a header. Each row describes one
QTL dataset to run.

| Column | Description |
|------|-------------|
| `susie_path` | Path to SuSiE fine-mapping results. Supports `gs://` or task-visible local paths |
| `ld_matrix_path` | Path to the matching LD matrix RDS. Supports `gs://` or task-visible local paths |
| `qtl_type` | QTL type label, such as `eQTL`, `sQTL`, or `pQTL`; passed to `--QTLType` when `GeneFilter` is supplied |

Example:

```text
susie_path	ld_matrix_path	qtl_type
gs://bucket/qtl/eqtl.susie.tsv.gz	gs://bucket/qtl/eqtl.ld.rds	eQTL
gs://bucket/qtl/pqtl.susie.tsv.gz	gs://bucket/qtl/pqtl.ld.rds	pQTL
```

### GWAS Manifest

The GWAS manifest is a tab-separated file with a header. Each row describes one
GWAS trait to run.

| Column | Description |
|------|-------------|
| `summary_stats_path` | Path to bgzipped, tabix-indexed GWAS summary statistics. Supports `gs://` or task-visible local paths |
| `summary_stats_index_path` | Path to the tabix index for `summary_stats_path`; localized next to the summary statistics as `<summary_stats_basename>.tbi` |
| `chosen_label` | Trait label passed to `--ChosenLabel` when `GeneFilter` is supplied |

Example:

```text
summary_stats_path	summary_stats_index_path	chosen_label
gs://bucket/gwas/monocyte.tsv.gz	gs://bucket/gwas/monocyte.tsv.gz.tbi	monocyte count
gs://bucket/gwas/leukocyte.tsv.gz	gs://bucket/gwas/leukocyte.tsv.gz.tbi	leukocyte quantity
```

If the QTL manifest has 2 rows and the GWAS manifest has 3 rows, the workflow
creates 6 pair-level TWAS tasks. Each pair receives the QTL row's `qtl_type` and
the GWAS row's `chosen_label`.

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `QTLManifest` | File | QTL manifest TSV |
| `GWASManifest` | File | GWAS manifest TSV |
| `OutputPrefix` | String | Prefix for pair names and aggregate outputs |
| `AnalysisType` | String | `conditional_rare` or `standard`; defaults to `conditional_rare` |
| `RareColumn` | String | Rare indicator column for conditional rare TWAS; defaults to `rare` |
| `GeneFilter` | File? | Optional gene filter TSV |
| `NumPrempt` | Int | Number of preemptible retries for pair-level TWAS tasks |
| `TWASMemory` | Int | Memory in GB for pair-level TWAS tasks |
| `TWASDiskGB` | Int | Disk in GB for pair-level TWAS tasks |
| `AggregateMemory` | Int | Memory in GB for aggregation |
| `AggregateThreads` | Int | CPU threads for aggregation |
| `AggregateDiskGB` | Int | Disk in GB for aggregation |

### Outputs

| Name | Type | Description |
|------|------|-------------|
| `PairManifest` | File | Expanded QTL x GWAS pair table |
| `PairResults` | Array[File] | One TWAS result per QTL x GWAS pair |
| `PairMetadata` | Array[File] | One metadata file per QTL x GWAS pair |
| `AllResults` | File | Aggregate across all QTL types and GWAS traits |
| `ByQTLTypeResults` | Array[File] | Aggregate files split by QTL type |
| `RunMetadata` | File | One metadata row per pair-level run |

## Runtime Image

The workflows use:

```text
ghcr.io/aou-multiomics-analysis/twas:main
```

The Docker image is built from `envs/Dockerfile` and published by the GitHub Actions workflow in `.github/workflows/docker-image.yml`.
