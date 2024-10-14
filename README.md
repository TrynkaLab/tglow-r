# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package strucutre is heavily inspired by Seurat.

*Very important note:* We make no claims on the statistical validity of applying some of the approaches on any given dataset and this package is "use at your own risk". As the HCI feature space is so diverse and to maintain flexibility you can in principle run any data through the pacakge but this also means you can easily end up violating statistical assumptions. If in doubt, reach out to your local friendly statistician for advice if any given method is valid.

# Installation & dependencies
> NOTE: While dependencies should be installed automatically, but this is currently untested so your milage may vary.

> NOTE: For now repo is private, make sure you are on VPN when calling this.

This will install the latest development version, we don't yet have a release, but for stability you can checkout a specific commit using the `ref` argument in `remotes::install_git()`
```
library(remotes)

remotes::install_git("https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core.git")
```

## On Sanger farm22 - latest dev version - reccomended
A version compatible with the tglow-r softpack module comes pre-installed in `/software/teamtrynka/installs/tglow-rlibs` and can be loaded as such. 

> NOTE: The version here changes often at the moment, so might not be the most stable.

``` 
module load HGI/softpack/groups/cell_activation_tc/tglow-r/6
```

Then launch R.
```
library(tglowr, lib="/software/teamtrynka/installs/tglow-rlibs")
```

Alternatively you can install if using R from the headnode or jammy64 directly through gitlab as above or using the cloned repo or using the latest code from `/software/teamtrynka/installs/tglow-r-core` in combination with `devtools::document()` and `devtools::build()`.


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

# Overview of workflow

1. Loading data into a TglowDataset, making sure the metadata and features are properly assigned.
2. QC at the image level to identify outlier images
3. QC features, to remove lowly varying ones and those with a lot of NA's
4. QC at the cell level, to identify outlier cells 
   1. Based on marker features and expection of what cells should look like (size, shape, intensity)
   2. Based on PCA outliers within a QC group which reflects a biological condition.
5. Normalization and scaling. (BoxCox transform, z-score / modified z)
6. Optional scaling to control samples
7. Optional aggregation to a grouping variable (mean, median, sum)
8. Covariate / batch regression (currently only linear or linear mixed models)
9. PCA and UMAP
10. Clustering (Louvain / Leiden based on ANNOY knn grapph)
11. Finding cluster markers using t-test
12. Finding associations using linear or linear mixed models

# Loading data into a TglowDataset
Currently built arround the output of the tglow-pipeline, will write some more general constructors soon.
Have a look at the help text for `read_cellprofiler_dir`, `tglow_dataset_from_list`, `read_cellprofiler_fileset_a` and `read_cellprofiler_fileset_b` 
as there might be some options & patterns to set depending how you export the data from cellprofiler.

>NOTE: Currently I'd reccomend type 'B' as it is the best tested, and the least work to setup in the cellprofiler pipeline and keeps relationships with children

```
path    <- "../../pipeline_disulfram/results/cellprofiler_v1"

# Read the data in the new format, merging strategy takes applies the function
# to the child objects, na.rm controls if NA's should be removed when calculating
# this.
output  <- read_cellprofiler_dir(path, pattern=".zip", type="B", 
                                 merging.strategy="mean", na.rm=T)
                                 
# Convert to tglow object
tglow <- tglow_dataset_from_list(output, assay="cells")

# Check if the dataset is valid
isValid(tglow)
```

The package comes with a bundled tglow object for testing which can be loaded with `data(tglow_example)`

## Checking validity of a TglowDataset
There are a couple of assumptions made downstream to enable functionality, that might not be met if the TglowDataset is improperly constructed.
Given construction is based on matching patterns, this can give issues sometimes if the default patterns are not appropriate for your data.

There is a utility method `isValid()` that you can use to check the key assumptions are met. The method is implemented for 
TglowDataset,TglowAssay,TglowMatrix and TglowReduction. When something is invalid, a warning is raised with more detailled info on what is wrong.
After fixing the issue, make sure to run `isValid()` again, as it returns FALSE after encountering its first problem in the hirearchy.

>NOTE: Invalidity of an assay might not mean that anything is terribly broken depending on what it is. It might just mean that some operations like 
slicing a dataset do not work as expected if rownames of a reduction are not set for example. In principle, anything generated through the proper 
functions should yield a valid Tglow class object. If you find this is not the case, please raise an issue. 


# Using TglowDataset
For more detaills also see the function definitions

## Data structure
Data is organized into a TglowDataset object, which stores image / well level metadata alongside the features. Features are stored in a slot called assays, which have the class TglowAssay. These are very similar to Seurat Assays. TglowAssay objects store the cell-feature level data and make a distinction between numeric data used for analysis which is stored as a matrix, and cell level metadata such as object IDs stored in a dataframe.

## Operations on TglowDataset

#### Show
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

#### Slicing

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

#### Accessing assays

You can acess TglowAssays from the `@assays` slot by using `$` or by `[[]]`
```
tglow@assays[["raw"]]
tglow$raw
tglow[["raw"]]
```

#### Accessing feature data
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
#### Setting assays

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

#### Accessing metadata
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

#### Getting and setting object IDs

Object ID's can be viewed in two ways
```
# Accessing the slot directly
tglow@object.ids

# Or using objectIds
objectIds(tglow)
```

The method `objectIds` is implemented for `TglowDataset`, `TglowAssay` and `TglowReduction`


To set object ID's you can use the `<-` operator. 

```
# Set object ID's for all assays, reductions and metadata
objectIds(tglow) <- paste0("O", 1:nrow(tglow))
```

> NOTE: You can also do it manually through the slots, but this is not reccomended, as it can lead to issues when not all slots are set properly, as it assumed all slots have the rownames set to enable easy slicing by object ID. Similarly, you could call `objectIds(tglow@assays[[1]]) <- 1:nrow(tglow)` but this will likely break downstream functionality.


#### Matching TglowDatasets together based on matching metadata

You can align two datasets on a matching ID. Filesets might not always be read in the same order, so the default `ObjectNumber_Global` id's are not guaranteed to match when using cellprofiler results from different runs. Given we use the same cellpose masks and if you have configured cellprofiler to NOT relabel cells  you can use the ObjectNumber to match between datasets

```
# Fetch the data you need to match in set 1
rn                   <- getDataByObject(tglow.new, c("plate_id", "well", "Metadata_field", "cell_ObjectNumber_Global"))
rn$ObjectNumber      <- gsub("FS\\d+_I\\d+_O(\\d+)", "\\1", rn$cell_ObjectNumber_Global)

# Assign the new ID's to set 1
objectIds(tglow.new) <- paste0(rn$plate_id, "_", rn$well, "_", rn$Metadata_field, "_", rn$ObjectNumber)

# Fetch the data you need to match in set 2
rn               <- getDataByObject(tglow, c("plate_id", "well", "Metadata_field", "cell_ObjectNumber_Global"))
rn$ObjectNumber  <- gsub("FS\\d+_I\\d+_O(\\d+)", "\\1", rn$cell_ObjectNumber_Global)

# Assign the new ID's to set 2
objectIds(tglow) <- paste0(rn$plate_id, "_", rn$well, "_", rn$Metadata_field, "_", rn$ObjectNumber)

# Overlap the ID's 
ol              <- intersect(tglow@object.ids, tglow.new@object.ids)

tglow.new <- tglow.new[ol,]
tglow     <- tglow[ol,]

```


### Setting up filters
Filters can be easily configured based on a filter table, making it easy to template sets of operations. Filters are NOT applied in order, but run independently. If you do want to run filters in order, you will have to run successive iterations, but this is easy enough to do. An easy way to maintain filters and edit them is to store them in a google sheet and load them into R. Then using the function `tglow_filters_from_table` to create the filter objects. The filter table should have the following columns, and one sheet for feature level filters, and one for object level filters. Exact layouts are customizable, see the help of `tglow_filters_from_table`


| Keyword         | Description                                                                                                                                                                    |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| name            | Filter name                                                                                                                                                                     |
| column_pattern  | What features to apply the filters to for feature filters, or what features to use to calculate the filters for object filters. The pattern 'all' is a special case that applies to all features |
| metadata_group  | Optional - Calculate filters within in a group of objects (not all filters respect this)                                                                                        |
| type            | The filter function above                                                                                                                                                       |
| value           | Threshold value passed ot filter                                                                                                                                                |
| transpose       | If transpose is true, data are first transposed, so the columns become the objects, not the features                                                                            |
| note            | Place to store extra info                                                                                                                                                       |
| active          | Should the filter be applied at runtime                                                                                                                                         |

#### Available filters

| Filter Types                  | Description                                                              | Respects Grouping | Note                                                                                       |
|-------------------------------|--------------------------------------------------------------------------|-------------------|--------------------------------------------------------------------------------------------|
| filter_coef_var               | Coefficient of variation                                                | FALSE             |                                                                                            |
| filter_coef_var_sum           | Coefficient of variation sum. All columns must pass                    | FALSE             |                                                                                            |
| filter_inf                    | Infinite values                                                         | FALSE             |                                                                                            |
| filter_inf_median             | Infinite median value                                                   | FALSE             |                                                                                            |
| filter_inf_median_sum         | Infinite median value sum. All columns must pass                       | FALSE             |                                                                                            |
| filter_inf_mutlicol           | Infinite values  - multiple columns                                      | FALSE             |                                                                                            |
| filter_max                    | Maximum value                                                           | FALSE             |                                                                                            |
| filter_max_sum                | Maximum value sum. All columns must pass                                | FALSE             |                                                                                            |
| filter_min                    | Minimum value                                                           | FALSE             |                                                                                            |
| filter_min_sum                | Minimum value sum. All columns must pass                                | FALSE             |                                                                                            |
| filter_mod_z                  | Absolute modified z-score < thresh                                      | TRUE              |                                                                                            |
| filter_mod_z_sum              | Absolute modified z-score < thresh sum. All columns must pass          | TRUE              |                                                                                            |
| filter_mod_z_perc             | Absolute modified z-score < thresh sum. Percentage of columns must pass  | TRUE              |                                                                                            |
| filter_na                     | NA filter                                                               | FALSE             |                                                                                            |
| filter_na_multicol            | NA filter - multiple columns                                            | FALSE             |                                                                                            |
| filter_near_zero_var          | Near zero variance from caret                                           | FALSE             | These are quite slow and intensive to compute, recommend filter_coef_var instead          |
| filter_near_zero_var_sum      | Near zero variance from caret sum. All columns must pass               | FALSE             | These are quite slow and intensive to compute, recommend filter_coef_var instead          |
| filter_unique_val             | Minimal number of unique values                                          | FALSE             |                                                                                            |
| filter_unique_val_sum         | Minimal number of unique values sum. All columns must pass              | FALSE             |                                                                                            |
| filter_zero_var               | Exactly 0 variance                                                      | FALSE             | Should also be covered by filter_coef_var                                                  |
| filter_zero_var_sum           | Exactly 0 variance sum. All columns must pass                           | FALSE             | Should also be covered by filter_coef_var                                                  |
| filter_blacklist              | Always FALSE                                                            | FALSE             |                                                                                            |






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
  
#### Data Import and Export

- read_cellprofiler_dir
- read_cellprofiler_fileset_a
- read_cellprofiler_fileset_b
- tglow_read_binmat
- tglow_read_imgs


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
- filter_blacklist
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
- filter_mod_z_sum
- filter_na
- filter_na_multicol
- filter_near_zero_var
- filter_near_zero_var_sum
- filter_sum
- filter_unique_val
- filter_unique_val_sum
- filter_zero_var
- filter_zero_var_sum

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
- find_outliers_pca_fixed
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