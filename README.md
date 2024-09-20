# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package strucutre is heavily inspired by Seurat.

# Installation & dependencies
Note: Will update a list of dependencies later and automate the dependencies during install.  

Note: For now repo is private, make sure you are on VPN when calling this.
```
library(remotes)

remotes::install_git("https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core.git")
```

## On Sanger farm22 - latest dev version

If using R from the headnode or jammy64 you should be able to install directly through gitlab, using the cloned repo or using the latest code from `/software/teamtrynka/installs/tglow-r-core`. If using Rstudio server, the install can be a bit funky, and the following workarround works by copying the repo to a local temp, and building from there into a local tmp library inside the container.

``` 
module load HGI/softpack/groups/cell_activation_tc/tglow-r/5
```

Then launch R.
The latest dev version is in `/software/teamtrynka/installs/tglow-r-core`
```
# Install the latest tglow R package in a tmpdir
tmpdir <- tempdir()
path   <- paste0(tmpdir, "/tglowr")
libdir <- paste0("/tmp/rlibs-tglow")

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

# Overview of workflow

1. Loading data into a TglowDataset, making sure the metadata and features are properly assigned.
2. QC at the image level to identify outlier images
3. QC features, to remove lowly varying ones and those with a lot of NA's
4. QC at the cell level, to identify outlier cells 
   1. Based on marker features and expection of what cells should look like (size, shape, intensity)
   2. Based on PCA outliers within a QC group which reflects a biological condition.
5. Normalization and scaling. (BoxCox transform, z-score / modified z)
6. Optional scaling to control samples
7. Optional aggregation (mean, median, sum)
8. Covariate / batch regression (currently only linear models)
9. PCA and UMAP
10. Clustering (Louvain / Leiden based on ANNOY knn grapph)
11. Finding cluster markers using t-test
12. Finding associations using linear models

# Using TglowDataset
For more detaills also see the function definitions

## Data structure
Data is organized into a TglowDataset object, which stores image / well level metadata alongside the features. Features are stored in a slot called assays, which have the class TglowAssay. These are very similar to Seurat Assays. TglowAssay objects store the cell-feature level data and make a distinction between numeric data used for analysis which is stored as a matrix, and cell level metadata such as object IDs stored in a dataframe.

## Operations on TglowDataset

##### Show
Show will show some usefull info on the object
```
> data(tglow_example)
> tglow

TglowData with:  5000  objects (cells),  360  images and  1  assays 
- Assays:
	$raw:	TglowAssay with: 5000 objects and 989 features. Size: 0.04 Gb 
- Active assay: raw 
- Reductions: 
- Object size: 0.05 Gb
```

##### Slicing

TglowDatasets and assays can be sliced by row
```
# Select the 50th cell
tglow[50,]
```

They can also be sliced by row, but this is not reccomended unless using column names, as assays can have different number of columns. A warning is raised when you try to slice columns with a non-character.
```
> tglow[50,34]
TglowData with:  1  objects (cells),  1  images and  1  assays 
- Assays:
	$raw:	TglowAssay with: 1 objects and 1 features. Size: 0 Gb 
- Active assay: raw 
- Reductions: 
- Object size: 0 Gb

Warning message:
In .local(x, i, j, ..., drop) :
  Assuming all assays have the same column order
```


TglowAssays can also be slices.
```
# Select the 50th cell
tglow$raw[50,1]
```

##### Accessing assays

You can acess TglowAssays from the `@assays` slot by using `$` or by `[[]]`
```
tglow@assays[["raw"]]
tglow$raw
tglow[["raw"]]
```

##### Accessing feature data
You can access the data using slicing `[]` or by using `$` on a TglowAssay or a TglowMatrix. If you use slicing on a TglowAssay, a new assay is returned.
If you use `$` a list with items `data` and `scale.data` is returned.
```
# Returns a slice as a TglowAssay
tglow$raw[50,1]

# Returns a slice as a TglowMatrix
tglow$raw@data[50,1]

# Returns a list with data and scale data items
tglow$raw$cell_AreaShape_Area

# Returns a numeric with the value
tglow$raw@data$cell_AreaShape_Area
```
##### Setting assays

You can set TglowAssays using the `@assays` slot. The assays slot is just a list, so you can put anything in it but if you want it to work properly
`new.assay` must be a TglowAssay
```
# Create a new assay from raw, with just the first 10 features
new.assay <- new("TglowAssay",
      data=tglow$raw@data[,1:10],
      scale.data=NULL,
      features=tglow$raw@features[1:10,]
)
tglow@assays[["new.assay"]] <- new.assay
```

##### Accessing metadata
Metadata is stored in two main places. @meta for cell level metadata, and @image.meta with image level metadata (not unique per object).
To get metadata per cell object there is a convience function `getDataByObject()` which grabs data from any assay, slot or metadata item
as long as its colnames are uniquely findable. There is also a slot @image.ids whose values should match rownames(image.meta) which indicates
which objects belong to which images. 

```
# Get drug, dose and cell size from metadata per object
data <- getDataByObject(tglow, c("drug", "dose", "cell_AreaShape_Area"), assay="raw", slot="data")

# Metadata can also be directly accessed
dim(tglow@meta)  
dim(tglow@image.meta)  

# Image data can also be accessed using
?getImageData()
?getImageDataByObject()
``` 


# List of functions


#### S4 Objects
- TglowDataset
- TglowAssay
- TglowMatrix
- TglowReduction
- TglowFilter


#### IO
- read_cellprofiler_dir
- read_cellprofiler_fileset_a
- read_cellprofiler_fileset_b


#### Constructors
- TglowMatrix
- tglow_dataset_from_list
- tglow_assay_from_list
- tglow_filters_from_table


#### Utils
- add_global_ids
- aggregate_assay
- aggregate_by_imagecol
- aggregate_metadata
- aggregate_tglow_matrix
- get_feature_meta_from_names
- merge_filesets
- getDataByObject
- getImageData
- getImageDataByObject


#### QC & Filtering functions
- apply_feature_filters
- apply_image_filters
- calculate_feature_filters
- calculate_object_filters
- find_outliers_pca
- find_outliers_pca_fixed


#### Transform functions
- fast_colscale
- boxcox_transform
- apply_boxcox
- mod_zscore
- scale_assay


#### Clustering & reductions
- apply_clustering
- calculate_pca
- calculate_umap

#### Regression functions
- lm_matrix
- calculate_lm
- apply_correction_lm
- find_markers

#### Filters
- filter_coef_var
- filter_coef_var_sum
- filter_inf
- filter_inf_median
- filter_inf_median_sum
- filter_inf_mutlicol
- filter_max
- filter_max_sum
- filter_min
- filter_min_sum
- filter_mod_z
- filter_mod_z_perc
- filter_na
- filter_na_multicol
- filter_near_zero_var
- filter_near_zero_var_sum
- filter_sum
- filter_unique_val
- filter_unique_val_sum
- filter_zero_var
- filter_zero_var_sum