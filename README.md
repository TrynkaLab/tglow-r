# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package structure is inspired by Seurat.

Please see the [wiki](https://github.com/TrynkaLab/tglow-r/wiki/) for a detailled walkthrough of the analsysis steps and capabilities of the package, some basic information is outlined below.

*Very important note:* We make no claims on the statistical validity of applying some of the approaches on any given dataset and this package is "use at your own risk". As the HCI feature space is so diverse and to maintain flexibility you can in principle run any data through the package but this also means you can easily end up violating statistical assumptions. If in doubt, reach out to your local friendly statistician for advice if any given method is valid.


# Installation

Please see the wiki for detailled [install instructions](https://github.com/TrynkaLab/tglow-r/wiki/Installation)

For the basic install (does not cover optional dependencies) and given blas, lapack, nlopt and libxml2 are available to the system (see install instructions otherwise):

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_git("https://github.com/TrynkaLab/tglow-r.git")
```

# Getting started

Please see the wiki [wiki](https://github.com/TrynkaLab/tglow-r/wiki/) for a detailled walkthrough of the analsysis steps and capabilities of the package.

A very minimal workflow on a qc'ed object might look like this:

```
data(tglow_example)

# Subset to the first 50 features on the raw assay
tglow@assays[["raw"]] <- tglow@assays$raw[,1:50]

# Select 2000 random cells
tglow <- tglow[sample(nrow(tglow), 2000),]

# Scale data
tglow@assays[["raw"]] <- scale_assay(tglow@assays$raw)

# Run PCA
tglow <- calculate_pca(tglow, assay="raw")

# Cluster cells
tglow <- calculate_clustering(tglow, "PCA.raw")

# UMAP
tglow <- calculate_umap(tglow, reduction="PCA.raw")

# Plot
tglow_dimplot(tglow, "UMAP.PCA.raw", ident="clusters_res_0.1")
tglow_dimplot(tglow, "UMAP.PCA.raw", ident="drug")
tglow_dimplot(tglow, "UMAP.PCA.raw", ident="donor")

# Find markers
markers <- find_markers(tglow, "clusters_res_0.1", assay="raw", slot="scale.data")

# Calculate drug effects with a linear mixed model
drug.effect <- calculate_lmm(tglow,
                             assay="raw",
                             slot="data",
                             covariates=c("drug", "Metadata_well", "donor"),
                             formula= ~ drug + (1|Metadata_well) + (1|donor))
```


# List of functions

> NOTE on naming scheme: Functions with the prefix `tglow_` either retlate to one of the Tglow s4 objects or are specific to the Tglow image processing pipeline. The functions without a spcific prefix should be more generically applicable. 

A full index of functions can be found in [Function List](vingettes/function_list.md)


# Authors
- Olivier Bakker
- Julie Matte
- Madeline Ohl