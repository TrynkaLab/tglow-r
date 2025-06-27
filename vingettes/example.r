# Load the package
library(tglowr)
library(data.table)
library(ggplot2)
library(cowplot)

# Used for the function as.color()
library(network)

#-------------------------------------------------------------------------------
# This example data is included with the package
data("tglow_example")

#-------------------------------------------------------------------------------
# Filter images
# Subset data to remove fields with less than 10 cells
tglow <- apply_image_filters(tglow, tglow@image.meta$Count_cell >= 10)

# BoxCox and scale image data
tglow <- apply_boxcox(tglow, "image.data")

# Note: I have not yet added the sample metadata to the object, so this is just blank for now
tglow@image.meta$qc_group <- paste0(
    tglow@image.meta$Metadata_plate,
    "_", tglow@image.meta$donor,
    "_", tglow@image.meta$drug
)

# Find outlier images
# Set the thresholds here to something low just to illustrate the point
# return.pcs can be used to make the per qc group plots if need be
res <- find_outliers_pca(tglow,
    assay = "image.data.trans",
    qc.group = tglow@image.meta$qc_group,
    pc.n = 5,
    return.pcs = T
)

# Example plot for QC group 1
plot(res$pcs[[1]]$pcs[, 1],
    res$pcs[[1]]$pcs[, 2],
    pch = 20,
    col = ifelse(res$pcs[[1]]$outliers, "orange", "blue")
)

# Apply the outlier filter
tglow <- apply_image_filters(tglow, !res$outliers)

tglow
#-------------------------------------------------------------------------------
# Filter Features
# Load filters (can be templated as CSV, or loaded from google sheets)
filter.table <- data.frame(
    name = c("na_all", "inf_median", "min_uniuqe_val", "coef_var", "zero_var"),
    column_pattern = c("all", "all", "all", "all", "all"),
    type = c("filter_agg_na", "filter_agg_inf_median", "filter_agg_unique_val", "filter_agg_coef_var", "filter_agg_zero_var"),
    value = c(0.100, NA, 2.000, 0.001, NA),
    note = c("10% max NA cells in a featue", "Filters features with infinite median", "minimal number of unique values", "Filter for a coefficent of variation", "Filter for zero variance"),
    active = c(TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
)

# Convert to tglow filter objects and calculate filters
filters <- tglow_filters_from_table(filter.table, active.col = 6)
res <- calculate_feature_filters(tglow, filters, assay = "raw", slot = "data")

# Check how many have been removed
rm.per.filter <- sort(nrow(res) - colSums(res), decreasing = T)
rm.per.filter

# Apply feature filters to all assays
tglow <- apply_feature_filters(tglow, res)

#-------------------------------------------------------------------------------
# Filter cells
# Load filters (can be templated as CSV, or loaded from google sheets)
filter.table <- data.frame(
    name = c(
        "NA filter", "No infitnte feature filter", "Cell size minor min", "Cell size minor max",
        "Cell size major min", "Cell size major max", "Neighbour touchign max",
        "Cell Equivalent Diameter", "Nucleus count min", "Nucleus count max"
    ),
    column_pattern = c(
        "all", "all", "cell_AreaShape_MinorAxisLength",
        "cell_AreaShape_MinorAxisLength", "cell_AreaShape_MajorAxisLength",
        "cell_AreaShape_MajorAxisLength",
        "cell_Neighbors_NumberOfNeighbors_Adjacent",
        "cell_AreaShape_EquivalentDiameter",
        "cell_Children_nucl_Count",
        "cell_Children_nucl_Count"
    ),
    metadata_group = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
    type = c(
        "filter_agg_na_multicol", "filter_agg_inf_mutlicol", "filter_vec_min_sum",
        "filter_vec_max_sum", "filter_vec_min_sum", "filter_vec_max_sum",
        "filter_vec_max_sum", "filter_vec_min_sum",
        "filter_vec_min_sum", "filter_vec_max_sum"
    ),
    value = c(0.1, 0.0, 10.0, 150.0, 10.0, 250.0, 4.0, 30.0, 1.0, 1.0),
    transpose = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    note = c(
        "Percentage of NA features in a cell",
        "No infinite values",
        NA, NA, NA, NA, NA, NA, NA, NA
    ),
    active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
)

# Convert to tglow filter objects
filters <- tglow_filters_from_table(filter.table, func.col = 4, thresh.col = 5, trans.col = 6, active.col = 8)
res <- calculate_object_filters(tglow, filters, assay = "raw", slot = "data")

# Check how many have been removed
rm.per.filter <- sort(nrow(res) - colSums(res), decreasing = T)
rm.per.filter

# Apply the cell folters
tglow <- tglow[rowSums(res) == ncol(res), ]

#-------------------------------------------------------------------------------
# Find outlier cellls
res <- find_outliers_pca(tglow,
    assay = "raw",
    qc.group = tglow@image.meta[tglow@image.ids, "qc_group"],
    return.pcs = T
)

tglow <- tglow[!res$outliers, ]

#-------------------------------------------------------------------------------
# Boxcox transform
tglow <- apply_boxcox(tglow, assay = "raw", assay.out = "trans")

#-------------------------------------------------------------------------------
# PCA - UMAP
tglow <- calculate_pca(tglow, assay = "trans", pc.n = 10)
tglow <- calculate_umap(tglow, reduction = "PCA.trans", pc.n = 10, downsample = NULL)

# Just a very quick and dirty example plot
um <- tglow@reduction$UMAP.PCA.trans@x
plot(um[, 1],
    um[, 2],
    col = as.color(getDataByObject(tglow, "cell_AreaShape_EquivalentDiameter", assay = "trans", drop = T)),
    pch = 20
)

#-------------------------------------------------------------------------------
# Clustering
tglow <- calculate_clustering(tglow, reduction = "PCA.trans", resolution = 0.2, k = 10, method = "louvain")

table(tglow@meta$clusters)

plot(um[, 1], um[, 2], col = tglow@meta$clusters, pch = 20)

#-------------------------------------------------------------------------------
# Find marker features for clusters
markers <- find_markers(tglow, "clusters", assay = "trans", slot = "scale.data")

boxplot(tglow@assays$trans@scale.data$cell_Intensity_MassDisplacement_ER ~ tglow@meta$clusters)

#-------------------------------------------------------------------------------
# Convert clusters to character label
tglow@meta$clusters <- as.character(tglow@meta$clusters)

# Regress out effects of covariates
covariates <- c("clusters")
tglow <- apply_correction_lm(tglow, "trans", slot = "scale.data", slot.covar = "scale.data", covariates = covariates)

# PCA after correction
tglow <- calculate_pca(tglow, assay = "trans.lm.corrected", pc.n = 10)
tglow <- calculate_umap(tglow, reduction = "PCA.trans.lm.corrected", pc.n = 10, downsample = NULL)

# Just a very quick and dirty example plot
um <- tglow@reduction$UMAP.PCA.trans.lm.corrected@x
plot(um[, 1], um[, 2], col = tglow@meta$clusters, pch = 20)

# Plot the boxplots of mass displacement
boxplot(tglow@assays$trans.lm.corrected@data$cell_Intensity_MassDisplacement_ER ~ tglow@meta$clusters)
boxplot(tglow@assays$trans@data$cell_Intensity_MassDisplacement_ER ~ tglow@meta$clusters)

# Check the means of mass displacement ER
aggregate(tglow@assays$trans.lm.corrected@data$cell_Intensity_MassDisplacement_ER, by = list(tglow@meta$clusters), FUN = mean)
aggregate(tglow@assays$trans@data$cell_Intensity_MassDisplacement_ER, by = list(tglow@meta$clusters), FUN = mean)

#-------------------------------------------------------------------------------
# Scale an assay relative to a control sample
tglow@assays[["trans.rel"]] <- scale_assay(tglow$trans,
    grouping = getDataByObject(tglow, "drug"),
    reference.group = tglow@meta$clusters == 3,
    scale.method = "median"
)

# After
df.plot <- getDataByObject(tglow, c(
    "drug",
    "clusters",
    "cell_Intensity_MassDisplacement_ER"
),
assay = "trans.rel",
slot = "scale.data"
)
p1 <- ggplot(data = df.plot, mapping = aes(x = clusters, y = cell_Intensity_MassDisplacement_ER)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_boxplot() +
    facet_wrap(~drug)

# Before
df.plot <- getDataByObject(tglow, c(
    "drug",
    "clusters",
    "cell_Intensity_MassDisplacement_ER"
),
assay = "trans",
slot = "scale.data"
)

p2 <- ggplot(data = df.plot, mapping = aes(x = clusters, y = cell_Intensity_MassDisplacement_ER)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_boxplot() +
    facet_wrap(~drug)

plot_grid(plotlist = list(p2, p1), nrow = 2)

#-------------------------------------------------------------------------------
# Find associations between metadata items and features
covariates <- c("drug", "dose")
res <- calculate_lm(tglow,
    assay = "trans",
    slot = "scale.data",
    slot.covar = "scale.data",
    covariates = covariates
)
View(res$model.stats)
#-------------------------------------------------------------------------------
# Aggregate whole tglow object on time_plate_well
tglow.agg <- aggregate_by_imagecol(tglow, "time_plate_well", method = "mean")
tglow.agg <- calculate_pca(tglow.agg, "trans")
tglow.agg <- calculate_umap(tglow.agg, reduction = "PCA.trans")
tglow.agg <- apply_clustering(tglow.agg, reduction = "PCA.trans", resolution = 1, k = 10, method = "louvain")

um <- tglow.agg@reduction$UMAP.PCA.trans@x

plot(um[, 1],
    um[, 2],
    col = tglow.agg@meta$clusters,
    pch = 20, cex = 2
)

find_markers(tglow.agg, ident = "clusters", assay = "trans", slot = "scale.data", return.top = 3)
#-------------------------------------------------------------------------------
