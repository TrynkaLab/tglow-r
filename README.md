# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package strucutre is heavily inspired by Seurat.

# Installation & dependencies
Currently, just source the scripts. Will update a list of dependencies later.

# Data structure
Data is organized into a TglowDataset object, which stores image / well level metadata alongside the features. Features are stored in a slot called assays, which have the class TglowAssay. These are very similar to Seurat Assays. TglowAssay objects store the cell-feature level data and make a distinction between numeric data used for analysis which is stored as a matrix, and cell level metadata such as object IDs stored in a dataframe.

# Workflow

1. Loading data into a tglow object, making sure the metadata and features are properly assigned.
2. QC at the image level to identify outlier images
3. QC features, to remove lowly varying ones and those with a lot of NA's
4. QC at the cell level, to identify outlier cells 
   1. Based on marker features and expection of what cells should look like (size, shape, intensity)
   2. Based on PCA outliers within a QC group which reflects a biological condition.
5. Normalization and scaling. (uses BoxCox transform)
6. Covariate / batch regression (currently only linear models)
7. Dimensionality reduction and UMAP
8. Clustering (Leiden based on ANNOY knn grapph)
