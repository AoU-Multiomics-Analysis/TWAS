# TWAS Pipeline

This repository contains R scripts and WDL workflows for running transcriptome-wide association studies (TWAS) from SuSiE fine-mapping results and GWAS summary statistics.

The main analysis combines fine-mapped molecular QTL weights with GWAS Z-scores to test whether genetically predicted molecular trait levels are associated with a GWAS trait. The repository also includes a conditional rare-variant TWAS follow-up that separates each gene model into common and rare components and tests whether the rare component contributes signal after accounting for the common component.

## Analysis Overview

The pipeline is organized around five steps.

1. Prepare per-gene phenotype inputs for optional cis-heritability estimation.
2. Compute an LD matrix for the fine-mapped variants used in each gene model.
3. Run standard SuSiE-weighted TWAS using all variants in the gene model.
4. Run conditional rare/common TWAS by partitioning the model into common and rare variants.
5. Aggregate per-trait or per-gene TWAS output files for downstream analysis.

At a high level, the standard TWAS statistic is a weighted sum of GWAS Z-scores:

```text
s = t(w) %*% z
var(s) = t(w) %*% LD %*% w
z_TWAS = s / sqrt(var(s))
```

where `w` is the SuSiE posterior mean effect vector, `z` is the GWAS Z-score vector, and `LD` is the variant correlation matrix.

The conditional rare TWAS extends this idea by using two predictors:

```text
w_common
w_rare
```

It reports the full TWAS, common-only TWAS, rare-only TWAS, a 2-df joint common+rare test, and a 1-df rare-conditional test. See [Conditional Rare TWAS](docs/conditional-rare-twas.md) for the statistical details and interpretation.

## Documentation

- [Script Reference](docs/scripts.md): command-line inputs and outputs for each R script.
- [Workflow Reference](docs/workflows.md): WDL workflow inputs, outputs, and runtime behavior.
- [Conditional Rare TWAS](docs/conditional-rare-twas.md): model definition, covariance calculations, conditional rare test, output columns, and interpretation.

## Repository Structure

```text
TWAS/
+-- scripts/          # R scripts used directly or by WDL tasks
+-- workflows/        # WDL workflow definitions
+-- docs/             # Detailed pipeline and method documentation
+-- envs/             # Docker runtime definition
+-- .github/
    +-- workflows/    # GitHub Actions CI configuration
```

## Core Inputs

The main TWAS scripts expect:

- SuSiE fine-mapping results with variant IDs, molecular trait IDs, posterior mean weights, genomic positions, and alleles.
- GWAS summary statistics with `CHR`, `POS`, `REF`, `ALT`, `BETA`, `SE`, and `Pvalue`.
- LD matrices for the variants used in each gene model.

For conditional rare TWAS, fine-mapping results should include either a rare-variant indicator column, usually `rare`, or a `gvs_max_af` column. If the rare indicator is absent, the script defines rare variants as `gvs_max_af < 0.01` after filtering variants with missing or non-numeric `gvs_max_af` values.

## Gene And Molecular Trait ID Handling

The TWAS and LD scripts use two related ID normalizations:

- Gene filtering extracts the embedded Ensembl gene ID and strips the GENCODE version.
- TWAS grouping and LD matching preserve the molecular trait ID, but strip embedded Ensembl gene versions.

This distinction matters for pQTL and splicing inputs, where multiple molecular traits can map to the same gene.

Supported examples:

```text
Input ID                                                   Gene filter ID      LD/TWAS ID
ENSG00000111647.13                                      -> ENSG00000111647    ENSG00000111647
chr10:100262063:100267571:clu_12500_-:ENSG00000095485.18 -> ENSG00000095485    chr10:100262063:100267571:clu_12500_-:ENSG00000095485
A0JNW5_ENSG00000111647.13                              -> ENSG00000111647    A0JNW5_ENSG00000111647
```

This lets a gene-level filter retain all matching molecular traits while avoiding LD-name collisions for multiple proteins or splice events from the same gene.

## Runtime

The WDL workflows use the Docker image built from [envs/Dockerfile](envs/Dockerfile). GitHub Actions builds and publishes this image to the GitHub Container Registry as:

```text
ghcr.io/aou-multiomics-analysis/twas:<tag>
```

The workflows currently reference:

```text
ghcr.io/aou-multiomics-analysis/twas:main
```
