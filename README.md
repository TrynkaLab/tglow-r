# Tglow: R package for analyzing HCI features

This repo contains an R package for analyzing (single cell) HCI imaging data. The package strucutre is heavily inspired by Seurat.

*Very important note:* We make no claims on the statistical validity of applying some of the approaches on any given dataset and this package is "use at your own risk". As the HCI feature space is so diverse and to maintain flexibility you can in principle run any data through the pacakge but this also means you can easily end up violating statistical assumptions. If in doubt, reach out to your local friendly statistician for advice if any given method is valid.

# Installation & dependencies
> NOTE: For now repo is private, make sure you are on VPN when calling this.

This will install the latest development version, we don't yet have a release, but for stability you can checkout a specific commit using the `ref` argument in `remotes::install_git()`

If you need to build some dependencies from source, make sure there is a BLAS/LAPACK, nlopt (nlopt), libxml2 (igraph) lib available if it isn't already, otherwise dependencies likely will not install. However this will depend heaviliy on your setup. Below is a minimal example using conda.

```
conda install -c conda-forge blas lapack nlopt 
conda install R
```

Then launch R
```
install.packages("remotes")

remotes::install_git("https://gitlab.internal.sanger.ac.uk/TrynkaLab/tglow-r-core.git")
```
This unlocks the core functionality

#### Installing suggested packages
To enable the suggested packages, manually install ggrastr, EBImage and RBioFormats
```
# Optional if installing ggrastr to enable rasterization of plots with many points
conda install -c conda-forge r-ragg
```

Then launch R
```
BiocManager::install("EBImage")
BiocManager::install("RBioFormats")
install.packages("ggrastr")
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
3. QC features using a dynamic, configurable filter system, to remove lowly varying ones and those with a lot of NA's
4. QC at the cell level using a dynamic, configurable filter system, to identify outlier cells 
   1. Based on customizable marker features and expection of what cells should look like (size, shape, intensity)
   2. Based on PCA outliers within a QC group which reflects a biological condition.
5. Feature transformation (BoxCox transform)
6. Normalization and scaling (z-score / modified z)
7. Optional scaling to control samples
8. Optional aggregation to a grouping variable (mean, median, sum)
9. Covariate / batch regression (currently only linear or linear mixed models)
   1.  In linear and linear mixed modes, regressions can be run on subgroups of objects
   2.  In linear mode, can correct groups of features for specific features (e.g. correct all nucleus features for nucleus_intensity and all mito feature for mito_intensity)
10. PCA and UMAP reductions
11. Clustering (Louvain / Leiden based on rcpp ANNOY knn or exact knn grapph)
12. Finding cluster markers using t-test
13. Finding associations using linear or linear mixed models

### On the whishlist
1. Integrate support for Milo
2. Add additional automated outlier filtering approaches beyond PCA
3. Improve and standardize plotting "ecosystem"
4. Add plate overview plots
5. Create configuration object for fetching standard columns like object x/y/z position, plate/well/row/col/field (currently configured on a per function basis)

# Loading data into a TglowDataset

> NOTE: The package comes with a bundled tglow object for testing which can be loaded with `data(tglow_example)` if you just want to play arround

The package is designed around manipulating image level and single object level data at the same time. This has as an advantage that storage heavy string based metadata don't get replicated unnecessarily in the matrix which starts to matter when processing hundreds of thousands to millions of objects in a typical HCI analysis. It also makes it easier to perform certain steps on the image level (such as QC etc) out of the box without aggregating. If you just have single cell level data, you could always add a dummy image assay, and the downstream functionaly should still work.

Data can easily be retrieved at both the image (data and metadata) and at the object (data and metadata) level using `getDataByObject()` which returns a data.frame with one row per object. More detaills on that below.

## Loading generic paired image + object data
For this example, I will asumme the folowing information is available. By default, checks on validtiy are done on all these objects, so at minimum a warning is raised if one of these does not meet assumptions.

- `objects`: A numeric matrix (not data.frame) with object level data, rows are objects, columns are imaging features. Rownames have unique object ids, colnames have unique feature ids
- `images`: A numeric matrix (not data.frame) with image level data, rows are images, columns are imaging features. Rownames have unique image ids, colnames have unique feature ids
- `image.ids`: A character vector of length `nrow(objects)` describing how rows in `images` connect to rows in `objects`. The values in this vector must match the rownames of the images.
- `object.meta` (optional): This is a optional data frame with any non-numeric/numeric metadata for each object. Must be in the same order as `objects` and assumes rownames are set to the same as `objects`, and columns are unique metadata items or imaging features, colnames must be set
- `image.meta` (optional): as `object.meta` but then paired to the `images` matrix

The minimal example would then look as follows:

```
tglow <- TglowDatasetFromMatrices(objects, images, image.ids)
```

## Loading from the tglow nextflow pipeline / cellprofiler results
To load cellprofiler features produced by the tglow-pipeline or any other cellprofiler pipeline, you can follow the steps below.
First have a look at the help text for `read_cellprofiler_dir`, `TglowDatasetFromList`, , `read_cellprofiler_fileset_a` and `read_cellprofiler_fileset_b` as there might be some options & patterns to set depending how you export the data from cellprofiler.

The function `read_cellprofiler_dir` scans the directory tree for a given pattern to find "filesets" which are individual exports from a cellprofiler instance. When using the tglow-pipeline there will be one of these filesets per imaging well. But this depends on the way your cellprofiler pipeline is setup. How exactly this is setup really doesn't matter for downstream analysis. 

By default, these readers matche certain patterns in feature names to decide what to put as metadata, and what to use as features. Any non-numeric variables are put as metadata, as a `TglowAssay` can only store numeric data.

>NOTE: Currently I'd reccomend type 'B' in combination with tglow-pipeline results, and the least work to setup in the cellprofiler pipeline and keeps relationships with children.

### Type A
This assumes that all the object level information is in one text file (controlled by `pat.cells`), and does not read and match child files that might be exported. It also assumes the image data is in a seperate file (controlled by `pat.img`), and uses the cellprofiler '_Image.txt' file to identify unique filesets.

```
# Read the dataset into a list
output  <- read_cellprofiler_dir("/path/to/results",
                                 pattern="_Image.txt",
                                 type="A", 
                                 pat.img = "_Image.txt", 
                                 pat.cells = "_cell.txt")

# Convert to tglow object
tglow <- TglowDatasetFromList(output, assay="cells")

# Check if the dataset is valid
isValid(tglow)
```
 
### Type B
Read the data in the type 'B' format, which is one .zip archive per fileset with a .tsv file for each object measured. One file is assumed to be the parent object (e.g. cells) and each child object (e.g. nuclei, cytoplasm etc) is in a seperate file and has a ID column linking back to the parent object. The merging strategy takes applies the function to the child objects in cases where one parent has multiple children. If the merging strategy is set, one matrix is returned with where each row is a parent object, and the child object values represent the mean. The parameter na.rm controls if NA's should be removed when calculating the child mean values.

```
output  <- read_cellprofiler_dir("/path/to/results",
                                 pattern=".zip",
                                 type="B", 
                                 merging.strategy="mean",
                                 na.rm=T)
                                 
# Convert to tglow object
tglow <- TglowDatasetFromList(output, assay="cells")

# Check if the dataset is valid
isValid(tglow)
```


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
Data is organized into a TglowDataset object, which stores image level metadata seperate from the features. Features are stored in a slot called assays, which have the class TglowAssay. These are structurally similar to Seurat Assays. TglowAssay objects store the numeric cell-feature level data. Any other cell level metadata not relevant for describing biology, such as object IDs or absolute locations of a cell bounding box, should be stored on the @meta slot. 

A single TglowDataset will always have the same objects across its assays, but assays can differ in features, making it easy to subset qc and manipulate featuresets on the same object. To subset on the image or object level, TglowDatasets can be sliced and subsetted in various ways. More details on that below.

The object structure is as follows

*TglowDataset*
- @assays: list of TglowAssay's
  - TglowAssay: stores feature level information, rows are objects, columns are features
    - @data: Matrix with object level features
    - @scale.data: Scaled version of @data (usually mean 0 variance 1, but other options are available)
    - @features: data frame storing feature metadata
- @meta: data.frame storing cell level metadata (id's, clusterings etc.)
- @image.meta: data.frame storing image level metadata (conditions, drugs, donors, well ids etc.)
- @image.data: TglowAssay for storing raw image features
- @image.data.trans: TglowAssay for storing BoxCox transformed image features
- @image.data.norm: TglowAssay for storing normalized image features
- @object.ids: Character vector with the id's of the objects
- @image.ids: Charachter vector with the images each object comes from
- @reduction: List to store reductions
  - TglowReduction: Stores PCA/UMAP in a semi standardized format
    - @x: Matrix with reduction coordinates, rows are objects, columns are dimensions
    - @sdev: SD of the components
    - @sdev_total: Total standard deviations
    - @object: Flexible slot for storing PCA/UMAP output objects should that be needed
- @graph: Slot for storing the kNN graph, currently not formalized
- @active.assay: Name of an active assay, currently not in use

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

They can also be sliced by column, but this is not reccomended unless using column names, as assays can have different number of columns. A warning is raised when you try to slice columns with a non-character.
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

Individual TglowAssays can also be sliced, in this case it is safe to use integers to select features.
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

#### Setting and maniuplating assays

You can set TglowAssays using the `@assays` slot. The assays slot is just a list, so you can put anything in it but if you want it to work properly
`new.assay` must be a TglowAssay. Note at the moment you can only add assays using `tglow@assays[["new.assay"]] <- new.assay` and not the other operators.
```
# Create a new assay from raw, with just the first 10 features
new.assay <- new("TglowAssay",
      data=tglow$raw@data[,1:10],
      scale.data=NULL,
      features=tglow$raw@features[1:10,]
)
tglow@assays[["new.assay"]] <- new.assay
```

You can manipulate specific slots in existing assays as well. For example, replacing the data slot with the scale.data slot.
```
tglow@assays[["new.assay"]]@data <- tglow$new.assay@scale.data
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
Filters can be easily configured based on a filter table, making it easy to template sets of operations. Filters are NOT applied seqeuntially, but run independently. If you do want to run filters in seqeuntially, you will have to run successive iterations, but this is easy enough to do. An easy way to maintain filters and edit them is to store them in a google sheet and load them into R. Then using the function `tglow_filters_from_table` to create the filter objects. The filter table should have the following columns, and one sheet for feature level filters, and one for object level filters. Exact layouts are customizable, see the help of `tglow_filters_from_table`

There are two flavors of filters:
- filter_vec_x: Accepts a vector and returns a logical vector of the same length (i.e. 'which objects for this feature are > 0')
- filter_agg_x: Accepts a vector and aggregates on a statistic and returns a single logical (i.e 'is the variance of this feature > 0')

Then there are the filter modifiers
- filter_vec_x_sum: Applies the filter to multiple columns, returning a logical of nrow(input), where T only if all columns for that row are T, otherwise F
- filter_agg_x_mutlticol: Applies a filter to data with multiple columns and returns a logical vector of ncol(input). If you want to apply these at the object level (i.e. 'filter objects with >x% of NA features'), make sure to set `transpose=T` in the filter definition, if you want to filter features (i.e. 'filter features with >x% of NA objects') leave `transpose=F`.

##### Example

I want to filter objects where _mito features have more then 50% NA's and overall features objects have no more then 10% NA's. Another example can be found in /vingettes/example.r
```
data("tglow_example")

filters <- list()
# Filter cells which have >50% NA in mitochondria features
filters[["mito.na"]]    <- new("TglowFilter",
                               name="mito.na",
                               column_pattern="_mito",
                               func="filter_agg_na_multicol",
                               threshold=0.5,
                               transpose=T)

# Filter cells which have >10% in any features
filters[["general.na"]] <- new("TglowFilter",
                               name="general.na",
                               column_pattern="all",
                               func="filter_agg_na_multicol",
                               threshold=0.1,
                               transpose=T)

res <- calculate_object_filters(tglow, filters, "raw")
```
#### Defining custom filters
You can also define custom filters at runtime by loading a new function into the global environment. Just make sure it has the following signature `function(vec, thresh, grouping)`

```
# Create a new filter function
my_filter <- function(vec, thresh, grouping=NULL) {
  return(vec == thresh)
}

# Add it as a filter object
filters[["my.filter"]] <- new("TglowFilter",
                               name="my.filter",
                               column_pattern="all",
                               func="my_filter",
                               threshold=10,
                               transpose=F)

```


#### Available filters
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


| filter types | description | respects grouping | note |
|--------------|-------------|-------------------|------|
| filter_agg_coef_var | Coefficient of variation | FALSE | |
| filter_agg_coef_var_multicol | Coefficient of variation - multiple columns | FALSE | |
| filter_agg_inf | Infinite values | FALSE | |
| filter_agg_inf_median | Infinite median value | FALSE | |
| filter_agg_inf_median_sum | Infinite median value sum. All columns must pass | FALSE | |
| filter_agg_inf_mutlicol | Infinite values  - multiple columns | FALSE | |
| filter_vec_max | Maximum value | FALSE | |
| filter_vev_max_sum | Maximum value sum. All columns must pass | FALSE | |
| filter_vec_min | Minimum value | FALSE | |
| filter_vec_min_sum | Minimum value sum. All columns must pass | FALSE | |
| filter_vec_mod_z | Absolute modified z-score < thresh | TRUE | |
| filter_vec_mod_z_sum | Absolute modified z-score < thresh sum. All columns must pass | TRUE | |
| filter_vec_mod_z_perc | Absolute modified z-score < thresh sum. Pecentage of columns must pass | TRUE | |
| filter_agg_na | NA filter | FALSE | |
| filter_agg_na_multicol | NA filter - multiple columns | FALSE | |
| filter_agg_unique_val | Minimal number of unique values | FALSE | |
| filter_agg_unique_val_multicol | Minimal number of unique values sum. - multiple columns | FALSE | |
| filter_agg_zero_var | Exactly 0 variance | FALSE | Should also be covered by filter_coef_var |
| filter_agg_zero_var_multicol | Exactly 0 variance sum.- multiple columns | FALSE | Should also be covered by filter_coef_var |
| filter_agg_blacklist | Always FALSE | FALSE | |

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