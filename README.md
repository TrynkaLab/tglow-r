# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package strucutre is heavily inspired by Seurat.

# Installation & dependencies
Note: Will update a list of dependencies later and fix the dependencies during install.

Clone the repo


## On Sanger farm22

If using R from the headnode or jammy64 you should be able to install directly through the cloned repo. If using Rstudio server, the install can be a bit funky, and the following workarround works by copying the repo to a local temp, and building from there. The Rlibs user must be set when calling `rstudio start` for this to work.

``` 
module load HGI/softpack/groups/cell_activation_tc/tglow-r/5
```

Then launch R.
The latest dev version is on `/software/teamtrynka/installs/tglow-r-core`
```
# Install tglow R package 
path <- tempdir()
path <- paste0(path, "/tglowr")

# Copy the package to a wirtable tempdir
dir.create(paste0(path))
file.copy(from = list.files("</path/to/repo>", full.names = TRUE), 
          to = path, recursive = TRUE)

# Build and document
devtools::document(pkg=path)
res  <- devtools::build(pkg=path, path=path)

# Find the user library to install in
lib  <- strsplit(Sys.getenv("R_LIBS_USER"), ":")[[1]][2]

# Remove the old version, install the new version, restart R
remove.packages("tglowr", lib=lib)
install.packages(res, lib=lib)
.rs.restartR()

# Load the package
library(tglowr, lib=lib)
```



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
