# Conditional Rare TWAS

The conditional rare TWAS analysis asks whether rare-variant expression-prediction signal contributes to the GWAS association after accounting for the common-variant expression-prediction signal for the same gene.

This is useful when a full TWAS model mixes common and rare fine-mapped variants. The full model can be significant because of common variants, rare variants, or both. The conditional rare test separates those components and asks whether the rare component carries residual association beyond the common component.

## Inputs to the Test

For one gene, let:

```text
z          = vector of GWAS variant Z-scores
LD         = variant LD/correlation matrix
w_common   = SuSiE posterior mean weights for common variants, zero elsewhere
w_rare     = SuSiE posterior mean weights for rare variants, zero elsewhere
```

Rare variants are defined from the configured rare indicator column, usually `rare`. If that column is absent, variants with missing or non-numeric `gvs_max_af` are removed and the fallback definition is:

```text
rare = gvs_max_af < 0.01
```

## Standard TWAS Score

For any weight vector `w`, the TWAS score is:

```text
s = t(w) %*% z
```

Under the null hypothesis of no GWAS association, the approximate variance of this score is:

```text
Var(s) = t(w) %*% LD %*% w
```

The standard TWAS Z-score is:

```text
z_TWAS = s / sqrt(Var(s))
```

The conditional rare script calculates this statistic three ways:

- `full_zscore`: all variants in the gene model.
- `common_zscore`: common variants only.
- `rare_zscore`: rare variants only.

## Two-Predictor Score Vector

The conditional analysis builds two scores:

```text
s_common = t(w_common) %*% z
s_rare   = t(w_rare) %*% z
```

Then it collects them into:

```text
s = c(s_common, s_rare)
```

Because both scores are linear functions of the same GWAS Z-score vector, they can be correlated through LD. Their covariance matrix is:

```text
V = rbind(
  c(t(w_common) %*% LD %*% w_common,
    t(w_common) %*% LD %*% w_rare),
  c(t(w_rare) %*% LD %*% w_common,
    t(w_rare) %*% LD %*% w_rare)
)
```

The output columns use these names:

```text
Vcc = V[1, 1] = Var(s_common)
Vrr = V[2, 2] = Var(s_rare)
Vcr = V[1, 2] = Cov(s_common, s_rare)
Vcor = Vcr / sqrt(Vcc * Vrr)
```

`Vcor` is not a phenotype correlation. It is the LD- and weight-induced correlation between the common and rare TWAS scores.

## Joint Common+Rare Test

The joint test asks whether the two-dimensional score vector has evidence of association:

```text
chisq_joint = t(s) %*% solve(V) %*% s
p_joint = pchisq(chisq_joint, df = 2, lower.tail = FALSE)
```

This is a 2-df omnibus test. It can detect signal in the common component, rare component, or both, but it does not by itself say which component is driving the association.

## Conditional Rare Test

The conditional rare test removes from the rare score the part that is linearly predictable from the common score:

```text
s_rare_cond = s_rare - Vcr / Vcc * s_common
```

This is the residual rare score after projecting `s_rare` onto `s_common`.

The coefficient is:

```text
beta = Vcr / Vcc
```

because the best linear projection of `s_rare` on `s_common` under the null is:

```text
E[s_rare | s_common] = Cov(s_rare, s_common) / Var(s_common) * s_common
```

After subtracting this projection, the residual is uncorrelated with the common score:

```text
Cov(s_rare - beta * s_common, s_common) = 0
```

Solving for `beta` gives:

```text
Cov(s_rare, s_common) - beta * Var(s_common) = 0
beta = Vcr / Vcc
```

## Conditional Rare Variance

The conditional rare variance is the variance of the residual rare score:

```text
var_rare_cond = Var(s_rare - Vcr / Vcc * s_common)
```

Using the variance rule:

```text
Var(A - bB) = Var(A) + b^2 Var(B) - 2b Cov(A, B)
```

with:

```text
A = s_rare
B = s_common
b = Vcr / Vcc
```

we get:

```text
var_rare_cond = Vrr + (Vcr / Vcc)^2 * Vcc - 2 * (Vcr / Vcc) * Vcr
              = Vrr - Vcr^2 / Vcc
```

This is the amount of rare-score variance left after accounting for its covariance with the common score.

## Conditional Rare Z-Score

The final conditional rare Z-score is:

```text
z_rare_cond = s_rare_cond / sqrt(var_rare_cond)
```

and the p-value is:

```text
p_rare_cond = pchisq(z_rare_cond^2, df = 1, lower.tail = FALSE)
```

This is a 1-df test of the rare component after accounting for the common component.

## Interpreting Patterns

Some useful patterns:

- `rare_zscore` and `z_rare_cond` are similar: rare signal is not strongly explained by the common predictor.
- `rare_zscore` is strong but `z_rare_cond` is weak: rare marginal signal overlaps with or is explained by the common predictor.
- `rare_zscore` is weak but `z_rare_cond` is strong: conditioning revealed rare signal that was masked by covariance with the common predictor.
- `common_zscore` and `z_rare_cond` have opposite signs: common and rare components may point in opposite association directions.
- `p_joint` is significant but `p_rare_cond` is not: the gene-level association may be driven mainly by the common component.
- `p_rare_cond` is significant but `p_joint` is not: the targeted 1-df rare test is stronger than the 2-df omnibus test; interpret carefully and inspect `n_rare`, `Vcor`, and variant-level evidence.

## One-Variant Rare Components

When `n_rare == 1`, the rare-only TWAS component is essentially the GWAS Z-score of that one rare variant, oriented by the sign of its prediction weight. The conditional rare statistic then asks whether that single variant's oriented association remains after accounting for covariance with the common predictor.

These hits can be biologically interesting, especially when the variant is fine-mapped with high PIP, but they should be interpreted as single-variant rare signals rather than polygenic rare burden-like signals.

## Output Status Values

The `status` column reports whether the two-predictor calculation succeeded for a gene.

Common values:

- `ok`: common and rare variants were both present, covariance terms were valid, and the statistics were computed.
- `need both common and rare variants`: the gene had only common or only rare variants after filtering.
- `non-positive predictor variance`: one of the predictor variances was zero or invalid.
- `singular covariance`: the joint covariance matrix could not be inverted or the conditional variance was invalid.
- `no GWAS variants`: no GWAS variants were extracted for the gene interval.
- `no LD matrix for gene`: no LD matrix was available for the cleaned gene ID.
- `no uniquely matched GWAS variants`: variants were found, but none passed allele matching and duplicate filtering.

## Practical Filters

For downstream rare-component analyses, a typical first-pass filter is:

```text
status == "ok"
n_rare > 0
is.finite(z_rare_cond)
is.finite(p_rare_cond)
```

Then choose a significance threshold appropriate for the analysis, for example FDR across all tested gene-trait pairs:

```text
p.adjust(p_rare_cond, method = "BH") < 0.05
```

It is often useful to stratify or annotate by:

- `n_rare == 1` versus `n_rare > 1`.
- `abs(Vcor)`, because high common-rare score correlation can make conditioning more influential.
- Whether `common_zscore` and `z_rare_cond` have the same or opposite sign.
- Fine-mapping support for the rare variants, such as PIP and credible set membership.
