#-------------------------------------------------------------------------------
#' Create a new TglowMatrix
#' @param matrix The data matrix
#' @export
TglowMatrix <- function(data) {
  new("TglowMatrix", data = data)
}

#-------------------------------------------------------------------------------
#' Convert a tglow list to a tglow assay
#'
#' @param output Tglow list object obtained from tglow.read.dir
#' @param assay Which assay to convert
#' @param meta.cols Which object level collumns in the features to consider as metadata
#' @param col.object Which collumn to use as the identifier. Must be in meta.cols
#'
#' @export
tglow.assay.from.list <- function(output, assay, meta.cols, col.object) {
  if (!col.object %in% meta.cols) {
    stop(paste0("Object ID column col.object must be in meta.cols but ", col.object, " is not"))
  }

  row.meta <- output[[assay]][, meta.cols]
  data <- as.matrix(output[[assay]][, !colnames(output[[assay]]) %in% meta.cols])
  rownames(data) <- row.meta[, col.object]

  features <- output$features[!rownames(output$features) %in% meta.cols, ]
  features$analyze <- T
  features <- features[colnames(data), ]
  rownames(features) <- colnames(data)

  assay <- new("TglowAssay",
    data = TglowMatrix(data),
    features = features
  )

  return(assay)
}

#-------------------------------------------------------------------------------
#' Tglow dataset from list
#'
#' @param output Tglow list object obtained from tglow.read.dir
#' @param assay Which assay to convert
#' @param meta.patterns Grep patterns in collumn names asspciated with object level metadata items
#' @param img.feature.patterns Grep patterns in image column names assocated with image level features to analyze
#' @param col.object The collumn name in the features which contains the per object object identifier
#' @param col.img.id The collumn name in the features which contains the per object image identifier
#' @param col.meta.img.id The collumn name in the image level data which contains the image id
#'
#' @returns A populated TglowDataset
#'
#' @export
tglow.dataset.from.list <- function(
    output,
    assay,
    meta.patterns = c("ImageNumber", "ObjectNumber", "Object_Number", "Parent"),
    img.feature.patterns = c("^AreaOccupied", "^ImageQuality", "^Granularity"),
    col.object = "cell_ObjectNumber_Global",
    col.img.id = "cell_ImageNumber_Global",
    col.meta.img.id = "ImageNumber_Global") {
  # Perform some sanity checks
  if (sum(c(assay, "meta") %in% names(output)) != 2) {
    stop("Assay and meta items not found in list")
  }

  if (!col.img.id %in% colnames(output[[assay]])) {
    stop("Image ID column not found in assay")
  }

  if (!col.meta.img.id %in% colnames(output[["meta"]])) {
    stop("Image ID column not found in metadata")
  }

  # Columns in output data considered metadata and not features
  meta.cols <- unlist(sapply(meta.patterns, grep, colnames(output[[assay]]), value = T))

  # Set main assay
  main.assay <- tglow.assay.from.list(output, assay = assay, meta.cols = meta.cols, col.object = col.object)
  dataset <- new("TglowDataset")
  dataset@assays[["raw"]] <- main.assay
  dataset@active.assay <- "raw"
  dataset@object.ids <- output[[assay]][, col.object]

  # Filter meta
  dataset@meta <- output[[assay]][, meta.cols]
  rownames(dataset@meta) <- output[[assay]][, col.object]

  # Image level data
  image.cols <- unlist(sapply(img.feature.patterns, grep, colnames(output$meta), value = T))
  dataset@image.data <- new("TglowAssay",
    data = TglowMatrix(as.matrix(output$meta[, image.cols])),
    features = get.feature.meta.from.names(image.cols)
  )
  rownames(dataset@image.data@data) <- output$meta[, col.meta.img.id]
  dataset@image.meta <- output$meta[, !colnames(output$meta) %in% image.cols]
  rownames(dataset@image.meta) <- output$meta[, col.meta.img.id]

  dataset@image.ids <- dataset@meta[, col.img.id]
  names(dataset@image.ids) <- dataset@meta[, col.object]
  return(dataset)
}
