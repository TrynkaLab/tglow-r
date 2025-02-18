# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package strucutre is inspired by Seurat.

Please also see the [wiki](https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core/-/wikis/home) for a detailled walkthrough of the analsysis steps and capabilities of the package.

*Very important note:* We make no claims on the statistical validity of applying some of the approaches on any given dataset and this package is "use at your own risk". As the HCI feature space is so diverse and to maintain flexibility you can in principle run any data through the pacakge but this also means you can easily end up violating statistical assumptions. If in doubt, reach out to your local friendly statistician for advice if any given method is valid.

# Installation from source

> NOTE: For Sanger farm22 install, see below

This will install the latest development version, we don't yet have a release, but for stability you can checkout a specific commit using the `ref` argument in `remotes::install_git()`

If you need to build some dependencies from source, make sure there is a BLAS/LAPACK, nlopt (nlopt), libxml2 (igraph) lib available if it isn't already, otherwise dependencies likely will not install. However this will depend heaviliy on your setup. Below is a minimal example using conda.

```
conda install -c conda-forge blas lapack nlopt libxml2
conda install R
```

Then launch R

> NOTE: For now repo is private, make sure you are on VPN when calling this.
```
install.packages("remotes")
install.packages("git2r")
install.packages("getPass")

# For now it is in a private git, so will need to provide credentials
remotes::install_git("https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core.git",
 credentials=git2r::cred_user_pass(readline(prompt="Username: "), getPass::getPass()))
```
This unlocks all of the core functionality

#### Installing suggested packages

##### ggrastr
To enable rasterization of plots with many points (this is fully optional, and plotting will work without it)
```
# Optional if installing ggrastr to enable rasterization of plots with many points
conda install -c conda-forge r-ragg 
```

Then launch R
```
install.packages("ggrastr")
```

##### EBImage and RBioFormats
These packages are required if you want to be able to retrieve example images of objects from .ome.tiff files
organized by /plate/row/col/field.ome.tiff (this is fully optional and the rest of the package will work without it)

```
# Optional if installing EBImage and RBioFormats to fetch example images
conda install -c conda-forge r-rcurl r-rjava fftw  
```

Then launch R
```
install.packages("BiocManager")

BiocManager::install("EBImage")
BiocManager::install("RBioFormats")
```

## On Sanger farm22 - latest dev version - reccomended
A version compatible with the tglow-r softpack module comes pre-installed in `/software/teamtrynka/installs/tglow-rlibs` and can be loaded as such. 


#### Using pre-installed version
> NOTE: The version here changes often at the moment, so might not be the most stable.

``` 
module load HGI/softpack/groups/cell_activation_tc/tglow-r/6
```

Then launch R.
```
library(tglowr, lib="/software/teamtrynka/installs/tglow-rlibs")
```

#### Installing into your personal library from git
Alternatively you can install if using R from the headnode or jammy64 directly through gitlab. This module should have all dependencies pre-installed.
```
module load HGI/softpack/groups/cell_activation_tc/tglow-r/6
```

Then launch R.
```
install.packages("remotes")
install.packages("git2r")
install.packages("getPass")

# For now it is in a private git, so will need to provide credentials
remotes::install_git("https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core.git",
 credentials=git2r::cred_user_pass(readline(prompt="Username: "), getPass::getPass()))
```

#### Installing from the sources directly
Alternatively you can also use the latest code from `/software/teamtrynka/installs/tglow-r-core` in combination with `devtools::document()` and `devtools::build()`.
This ensures you have the aboslute latest version of the code (which is both good and bad), the pre-installed version or the git version might lag behind a bit depending on the commit / build frequency.

```
module load HGI/softpack/groups/cell_activation_tc/tglow-r/6
```

Then launch R.
```
# Install the latest tglow R package in a library of your choiche
tmpdir <- "/path/to/install/dir"
path <- paste0(tmpdir, "/build-dir")
libdir <- paste0(tmpdir)

dir.create(path)
dir.create(libdir)

# Copy the package to a wirtable tempdir
file.copy(
    from = list.files("/software/teamtrynka/installs/tglow-r-core", full.names = TRUE),
    to = path, recursive = TRUE
)

# Build and document
devtools::document(pkg = path)
res <- devtools::build(pkg = path, path = path)

# Remove the old version, install the new version, restart R
remove.packages("tglowr", lib = libdir)
install.packages(res, lib = libdir)

system(paste0("chmod -R 770 ", libdir))

```

##  On Sanger farm22 - Installing on Rstudio server - not reccomended (stop gap solution)
If using Rstudio server, the install can be a bit funky, and the following workarround works by copying the repo to a local temp, and building from there into a local tmp library inside the container. Alternatively you can ofcourse install into your personal library as well prior to launching Rstudio Sever and provide that in RStudio start, but this is a bit annoying for development, as you would need to re-start when re-installing.

> NOTE: Long term this is not reccomended as it produces overheads and can fill up /tmp in case of crashes. You can also consider installing into a /lustre folder as a very hacky solution that is more persistent

``` 
module load HGI/softpack/groups/cell_activation_tc/tglow-r/6
```

Then launch R.
The latest dev version is in `/software/teamtrynka/installs/tglow-r-core`
```
# Install the latest tglow R package in a tmpdir
tmpdir <- tempdir()
path   <- paste0(tmpdir, "/tglowr")
libdir <- paste0(tmpdir, "/rlibs-tglow")

dir.create(path)
dir.create(libdir)

# Copy the package to a wirtable tempdir
file.copy(from = list.files("/software/teamtrynka/installs/tglow-r-core", full.names = TRUE), 
          to = path, recursive = TRUE)

# Build and document
devtools::document(pkg=path)
res  <- devtools::build(pkg=path, path=path)

# Remove the old version, install the new version, restart R
remove.packages("tglowr", lib=libdir)
install.packages(res, lib=libdir)
```

You can then load the library using `library("tglowr", lib=paste0(tempdir(), "/rlibs-tglow"))` for as long as the session persists.

# List of functions

Note on naming scheme: Functions with the prefix `tglow_` either retlate to one of the Tglow s4 objects or are specific to the Tglow image processing pipeline.
The functions without a spcific prefix should be more generically applicable. 

#### S4 Objects
- TglowDataset
- TglowAssay
- TglowMatrix
- TglowReduction
- TglowFilter

#### Generic Methods

- colnames
- ncol
- nrow
- rownames

#### Data Structures and Manipulation

- objectIds
- objectIds<-
- isAvailable
- isValid
- getDataByObject
- getImageData
- getImageDataByObject

### Constructors
- TglowMatrix
- TglowAssayFromList
- TglowAssayFromMatrix
- TglowDatasetFromList
- TglowDatasetFromMatrices

#### Data Import and Export

- read_cellprofiler_dir
- read_cellprofiler_fileset_a
- read_cellprofiler_fileset_b
- tglow_read_binmat
- tglow_read_imgs
- tglow_read_imgs_aics


#### Data Aggregation and Merging

- add_features_to_assay
- add_global_ids
- aggregate_assay
- aggregate_by_imagecol
- aggregate_metadata
- aggregate_tglow_matrix
- merge_filesets

#### Data Transformation and Scaling

- apply_boxcox
- boxcox_transform
- fast_colscale
- mod_zscore
- scale_assay

#### Filtering and Data Cleaning

- apply_feature_filters
- apply_image_filters
- calculate_feature_filters
- calculate_object_filters
- filter_agg_coef_var
- filter_agg_coef_var_multicol
- filter_agg_inf
- filter_agg_inf_median
- filter_agg_inf_median_sum
- filter_agg_inf_mutlicol
- filter_vec_max
- filter_vev_max_sum
- filter_vec_min
- filter_vec_min_sum
- filter_vec_mod_z
- filter_vec_mod_z_sum
- filter_vec_mod_z_perc
- filter_agg_na
- filter_agg_na_multicol
- filter_agg_unique_val
- filter_agg_unique_val_multicol
- filter_agg_zero_var
- filter_agg_zero_var_multicol
- filter_agg_blacklist

#### Statistical Analysis

- correct_lm
- correct_lmm
- correct_lm_per_featuregroup
- calculate_lm
- calculate_lmm
- calculate_pca
- calculate_umap
- calculate_clustering
- find_markers
- find_outliers_pca
- lm_matrix
- lmm_matrix

#### Image Processing

- img_composite
- img_max_per_channel
- img_max_project
- img_norm

#### Plotting and Visualization

- apply_color
- hex_to_rgb
- plot_boxline
- plot_hex
- plot_hist_dens_grouped
- plot_img
- plot_img_set
- plot_simple_hm
- plot_xy
- tglow_plot_execution_time
- tglow_plot_location_hex
- tglow_dimplot
- theme_plain

#### Utility Functions

- nearest_index
- tglow_assay_from_list
- tglow_build_img_index
- tglow_dataset_from_list
- tglow_filters_from_table
- get_feature_meta_from_names
- fetch_representative_object
- fetch_representative_object_quantiles
- list_has_overlap
    

# Authors
- Olivier Bakker
- Julie Matte
- Madeline Ohl