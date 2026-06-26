# Code Review — tglow-r v0.1.20
**Date:** 2026-06-25
**Scope:** `calculate_eigenfeatures`, `find_eigenmarkers`, `plot_eigenmarkers` (last commit: 653bb53)

---

## Finding 1 — `"ward.D"` crashes in `calculate_feature_clustering` inner dispatch
**File:** `R/clustering.r` ~line 259 | **Severity: High**

The outer guard at line ~239 accepts `method %in% c("complete", "ward.D", "ward.D2")`, runs `hclust`, and stores the result. But the inner dispatch inside the resolution loop checks `c("complete", "ward", "ward.D2")` — using `"ward"` instead of `"ward.D"`. Passing `method="ward.D"` falls through to `stop("No valid cluster method")` every time, despite being documented as valid. `find_eigenmarkers` calls this function, so this bug surfaces there too.

---

## Finding 2 — NA cluster values propagate silently through `calculate_eigenfeatures`
**File:** `R/clustering.r` ~line 303 | **Severity: High**

```r
clusters <- unique(dataset@assays[[assay]]@features[, cluster.col])
```

No `na.omit()` or `!is.na()` guard. If any feature lacks a cluster assignment, `NA` enters `clusters`, becomes a column name in the output matrix, and causes `matrix[, NA]` / `meta[NA, ...]` assignments to silently corrupt results or produce a malformed assay with no error.

---

## Finding 3 — PCA sign indeterminacy in eigenfeatures is not anchored
**File:** `R/clustering.r` ~line 328–331 | **Severity: High**

```r
pca <- irlba::prcomp_irlba(data[,selector,drop=F], n=1, center=F, scale=F)
matrix[, cluster] <- pca$x[, 1]
```

SVD-based PCA has sign indeterminacy: PC1 direction is arbitrary. No sign-fixing (e.g. `* sign(mean(pca$x[,1]))`) is applied. Two runs on identical data — or across platforms / `irlba` versions — can produce opposite-signed eigenfeatures. Because `plot_eigenmarkers` uses a diverging red/blue colour scale, the same cluster can appear as a "positive" or "negative" marker purely from numerical noise.

---

## Finding 4 — `center=F, scale=F` likely captures mean rather than variance structure
**File:** `R/clustering.r` ~line 328 | **Severity: Medium-High**

The commented-out fallback used `center=TRUE, scale.=TRUE`. The live call uses neither. For HCI intensity data (non-zero means), PC1 without centering primarily captures the mean offset of the feature block rather than the dominant axis of covariation between features. Two phenotypically distinct cell populations with similar mean intensity levels could produce a non-significant eigenmarker not because there is no signal, but because the eigenfeature is poorly constructed.

---

## Finding 5 — `groups.x` unconditionally overwritten, breaking `...`-passed grouping
**File:** `R/plotting.r` ~line 1485 | **Severity: Medium**

```r
df.plot <- plot_markers(markers.eigen, return.data = TRUE, ...)$data
df.plot$groups.x <- df.plot$feature   # unconditional overwrite
```

Any `grouping.x` passed in `...` is processed by `plot_markers` and stored in `df.plot$groups.x`, then immediately discarded. The `tapply`-based clustering and the axis ordering then both run on `feature` regardless of what the caller requested, making `grouping.x` silently ignored in `plot_eigenmarkers`.

---

## Finding 6 — Silent NA injection when `markers$feature` not in `eigenmarkers$features`
**File:** `R/plotting.r` ~line 1471 | **Severity: Medium**

```r
markers$clust <- eigenmarkers$features[markers$feature, eigenmarkers$clustering.col]
```

`markers$feature` is used as a row-index into `eigenmarkers$features`. If any feature name is absent (e.g. `markers` comes from a different assay than the one used to build `eigenmarkers`), R silently returns `NA`. This `NA` propagates into `aggregate()` as a spurious cluster group and ultimately into the heatmap with no warning.

---

## Finding 7 — Hardcoded positional column index `[, 3]`
**File:** `R/plotting.r` ~line 1481 | **Severity: Medium**

```r
markers.eigen$estimate <- marker.means[rownames(markers.eigen), 3]
```

`aggregate()` with two grouping variables returns columns `feature`, `class`, `x` — column 3 is `x` today, but this is a fragile positional contract. If a grouping variable is added or the call is refactored, column 3 silently shifts to a grouping column and `estimate` is populated with cluster IDs instead of numeric effect sizes. Should be `marker.means[rownames(markers.eigen), "x"]`.

---

## Finding 8 — `cluster.col` existence not validated before use
**File:** `R/clustering.r` ~line 303 | **Severity: Medium**

`calculate_eigenfeatures` does not check whether `cluster.col` exists as a column in the features slot before indexing into it. The peer function `calculate_feature_clustering` (line ~180) validates `feature.group` with an explicit `%in% colnames(...)` guard and a descriptive `stop()`. Passing a nonexistent `cluster.col` produces a cryptic `subscript out of bounds` with no hint about which argument was wrong.

---

## Finding 9 — Division by zero when `pca$totalvar == 0`
**File:** `R/clustering.r` ~line 335 | **Severity: Low-Medium**

```r
meta[cluster, "var_exp"] <- pca$sdev[1]^2 / pca$totalvar
```

No guard against `totalvar == 0`. If a cluster's features are constant across all cells, `totalvar` is 0 and the result is `NaN`, which silently enters the features slot of the new assay. Any downstream filtering or display logic using `var_exp` will operate on `NaN` without warning.

---

## Finding 10 — `matrix` variable shadows `base::matrix()`
**File:** `R/clustering.r` ~line 306 | **Severity: Low**

```r
matrix <- matrix(NA, nrow = ..., ncol = nclust)
```

The local variable `matrix` shadows the base R constructor for the rest of the function. The current code doesn't call `matrix()` again inside the function, so there is no runtime failure today — but any future edit that tries to construct a matrix inside this function will silently operate on the pre-computed data object instead. Rename to e.g. `eigen.mat` or `result.mat`.

---

## Additional notes

- `df.plot$size <- -log10(df.plot$padj)` in `plot_eigenmarkers` is dead code — `size` is never mapped in the `aes()`. The docstring's claim that significance drives point sizing is inaccurate for the tile plot.
- There is a commented-out `prcomp` line left in `calculate_eigenfeatures` (line ~327) that should be removed.
- `return.object=F` uses bare `F` rather than `FALSE` — the linter is configured to ignore this, but it is a latent global-variable-masking risk.
