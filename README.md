# TWAS Pipeline

This repository contains scripts and WDL workflows for running a **Transcriptome-Wide Association Study (TWAS)** using SuSiE fine-mapping results and GWAS summary statistics. The pipeline estimates SNP heritability, computes LD matrices, runs the TWAS analysis, and aggregates results across genes.

---

## Repository Structure

```
TWAS/
├── scripts/          # R scripts used by the WDL workflows
├── workflows/        # WDL workflow definitions
├── envs/             # Dockerfile for the runtime environment
└── .github/
    └── workflows/    # GitHub Actions CI configuration
```

---

## Scripts

### `CreatePhenoFile.R`

Extracts expression data for a single gene from a BED-format expression matrix and prepares the input files required for heritability estimation.

**Inputs:**
- `--BedFile` – Path to a BED-format gene expression file (tab-separated, with columns `chrom`, `start`, `end`, `gene_id`, and per-sample expression values).
- `--PhenotypeID` – Gene ID to extract (must match a `gene_id` value in the BED file).

**Outputs:**
- `gene_region.tsv` – A three-column file (`chrom`, `start`, `end`) defining the ±1 Mb cis-window around the gene.
- `pheno.txt` – A phenotype file (`FID`, `IID`, expression value) compatible with PLINK/GCTA.

---

### `compute_LD_TWAS.R`

Computes the linkage disequilibrium (LD) matrix for the variants associated with a specific gene, using a dosage matrix. The LD matrix is used downstream in the TWAS Z-score calculation.

**Inputs:**
- `--DoseMatrix` – Path to a tab-separated dosage matrix (columns: `CHROM`, `POS`, `REF`, `ALT`, followed by per-sample dosages).
- `--PhenotypeID` – Gene ID whose fine-mapped variants should be used.
- `--VariantList` – Path to a file listing variants with at least a `phenotype` column and a `variant` column (output from fine-mapping).

**Outputs:**
- `<PhenotypeID>.LD.rds` – An RDS file containing the LD correlation matrix for the gene's cis-variants.

---

### `susie_TWAS.R`

Performs the core TWAS analysis for a gene. For each provided set of GWAS summary statistics, it integrates SuSiE fine-mapping posterior means with GWAS Z-scores and an LD matrix to compute a TWAS Z-score and p-value.

**Inputs:**
- `--LDMatrix` – Path to an RDS file (output of `compute_LD_TWAS.R`) containing a named list of LD matrices, one per gene.
- `--SusieRes` – Path to a tab-separated file of SuSiE fine-mapping results (must include columns `molecular_trait_id`, `variant`, `posterior_mean`, `chromosome`, `position`, `ref`, `alt`, `pip`, `cs_id`).
- `--SummaryStats` – Path to a bgzipped, tabix-indexed GWAS summary statistics file (columns: `CHR`, `POS`, `REF`, `ALT`, `BETA`, `SE`, `Pvalue`).
- `--OutputPrefix` – Prefix for the output file name (typically the GWAS name).

**Outputs:**
- `<OutputPrefix>.TWAS.txt` – A tab-separated file with columns `zscore`, `pvalue`, `stat`, `gene`, and `GWAS` for each gene.

---

### `aggregate_TWAS.R`

Aggregates multiple per-gene TWAS result files (output of `susie_TWAS.R`) into a single merged TSV file.

**Inputs:**
- `--FilePaths` – Path to a plain-text file listing the paths to individual `.TWAS.txt` result files (one per line).
- `--OutputPrefix` – Prefix for the output merged file.

**Outputs:**
- `<OutputPrefix>_TWAS.tsv` – A single merged TSV file containing all TWAS results.

---

## Workflows

### `ComputeHeritability.wdl`

**Workflow:** `EstimateHeritability`

Estimates the cis-SNP heritability (h²) of a gene's expression using PLINK2 and GCTA REML. For a given gene, it: (1) extracts expression and defines the cis-window using `CreatePhenoFile.R`, (2) subsets genotype data to the cis-window using PLINK2, (3) constructs a genetic relationship matrix (GRM), and (4) runs restricted maximum likelihood (REML) to estimate h².

**Inputs:**

| Name | Type | Description |
|------|------|-------------|
| `pvar` | File | PLINK2 `.pvar` variant information file |
| `pgen` | File | PLINK2 `.pgen` genotype file |
| `psam` | File | PLINK2 `.psam` sample information file |
| `ExpressionBed` | File | BED-format gene expression matrix |
| `Covars` | File | Quantitative covariates file for GCTA REML |
| `PhenotypeID` | String | Gene ID to compute heritability for |
| `Memory` | Int | Memory in GB for the runtime |
| `NumPrempt` | Int | Number of preemptible retries |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `HeritabilityEstimate` | File | GCTA `.hsq` file containing the heritability estimate |

---

### `compute_LD_TWAS.wdl`

**Workflow:** `PreprocessLD`

Computes the LD matrix for a gene's cis-variants by running `compute_LD_TWAS.R`. This is a preprocessing step required before running the TWAS analysis.

**Inputs:**

| Name | Type | Description |
|------|------|-------------|
| `DoseMatrix` | File | Tab-separated genotype dosage matrix |
| `VariantList` | File | File listing variants to include, with `phenotype` and `variant` columns |
| `PhenotypeID` | String | Gene ID for which to compute the LD matrix |
| `Memory` | Int | Memory in GB for the runtime |
| `NumPrempt` | Int | Number of preemptible retries |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `MatrixLD` | File | RDS file containing the LD correlation matrix for the gene |

---

### `susie_TWAS.wdl`

**Workflow:** `TWAS`

Runs the TWAS analysis for a gene against a single set of GWAS summary statistics by executing `susie_TWAS.R`. Requires the LD matrix (from `compute_LD_TWAS.wdl`) and SuSiE fine-mapping results as input.

**Inputs:**

| Name | Type | Description |
|------|------|-------------|
| `LDMatrix` | File | RDS LD matrix file (output of `compute_LD_TWAS.wdl`) |
| `SumStats` | File | Bgzipped, tabix-indexed GWAS summary statistics file |
| `SumStatsIndex` | File | Tabix index (`.tbi`) for the summary statistics file |
| `FineMapping` | File | SuSiE fine-mapping results file |
| `NameGWAS` | String | Name/prefix for the GWAS (used in output file naming) |
| `Memory` | Int | Memory in GB for the runtime |
| `NumPrempt` | Int | Number of preemptible retries |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `ResTWAS` | File | Tab-separated TWAS results file (`<NameGWAS>.TWAS.txt`) |

---

### `aggregate_TWAS.wdl`

**Workflow:** `AggregateTWASWorkflow`

Aggregates per-gene TWAS result files into a single merged TSV. It downloads result files from Google Cloud Storage using `gsutil` and then runs `aggregate_TWAS.R` to combine them.

**Inputs:**

| Name | Type | Description |
|------|------|-------------|
| `TWASResFOFN` | File | File of file names (FOFN): a plain-text file listing GCS paths to individual `.TWAS.txt` result files |
| `OutputPrefix` | String | Prefix for the merged output file |
| `Memory` | Int | Memory in GB for the runtime |
| `NumThreads` | Int | Number of CPU threads |

**Outputs:**

| Name | Type | Description |
|------|------|-------------|
| `AggregatedTWAS` | File | Merged TSV file containing all TWAS results (`<OutputPrefix>_TWAS.tsv`) |
| `filelist` | File | List of localized file paths used in the aggregation |

---

## GitHub Actions

### `docker-image.yml` — Docker Image CI

Builds and pushes the Docker image to the GitHub Container Registry (`ghcr.io`) whenever code is pushed to or a pull request is opened against the `main` branch.

The image is built from `envs/Dockerfile` and published as `ghcr.io/aou-multiomics-analysis/twas:<tag>`. This image is used as the runtime environment in all WDL workflows.
