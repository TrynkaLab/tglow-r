#' @include utils.r
NULL

#-------------------------------------------------------------------------------
#' Create a new TglowMatrix
#' @param matrix The data matrix
#' @export
TglowMatrix <- function(matrix) {
  if (!is.matrix(matrix)) {
    stop("Data must be matrix")
  }

  if (!is.numeric(matrix)) {
    stop("Data must me numeric")
  }

  new("TglowMatrix", data = matrix)
}

#-------------------------------------------------------------------------------
#' Create a new TglowDataset
#'
#' @description
#' Create a new TglowDataset object from a matrix with image level data and a matrix of object level data
#'
#' @param objects The object level data matrix, must be a matrix, and be numeric
#' @param image The image level data matrix, must be a matrix, and be numeric
#' @param image.ids The of nrow(objects) de which rownames of image corresponds to which object
#' @param object.meta A data.frame with metadata, one row per object. Defaults to data.frame(id=rownames(objects))
#' @param image.meta A data.frame with metadata, one row per image. Defaults to data.frame(id=rownames(images))
#' @param assay.out The assay name to store objects under
#'
#' @details
#' Only the data slot on assay is populated with the values in `objects` if this is not what you want you can
#' always re-assign data to the scale.data slot afterwards by `dataset$raw@scale.data <- dataset$raw@data`.
#' Just to note, the assumption is made that the scale.data slot is actually scaled to mean=0 variance=1 when
#' running PCA's etc, so overriding with non-scaled data may change your results.
#'
#' @returns A populated TglowDataset
#' @export
TglowDatasetFromMatrices <- function(objects, images, image.ids, object.meta = NULL, image.meta = NULL, assay.out = "raw") {
  stop("Not yet implemented")

  if (!check_dimnames(objects)) {
    stop("row and colnames on objects cannot be NULL")
  }

  if (!check_dimnames(images)) {
    stop("row and colnames on image cannot be NULL")
  }

  dataset <- new("TglowDataset")
  dataset@assays[[assay.out]] <- TglowAssayFromMatrix(objects)

  dataset@image.data <- TglowAssayFromMatrix(images)

  dataset@object.ids <- rownames(objects)
  names(dataset@object.ids) <- dataset@object.ids

  dataset@image.ids <- image.ids
  names(dataset@image.ids) <- dataset@object.ids

  if (is.null(object.meta)) {
    dataset@meta <- data.frame(id = rownames(objects), row.names = rownames(objects))
  } else {
    datset@meta <- object.meta
  }

  if (is.null(image.meta)) {
    dataset@image.meta <- data.frame(id = rownames(images), row.names = rownames(images))
  } else {
    datset@image.meta <- image.meta
  }

  if (!tglor::isValid(dataset)) {
    warning("Dataset did not pass checks in isValid, returning anyway")
  }

  return(dataset)
}

#-------------------------------------------------------------------------------
#' Create a new TglowDataset
#' @param objects The object level matrix
#' @param objects The scaled object level matrix (optional)
#' @param features The data frame with features
#' 
#' @returns A TglowAssay
#' @export
TglowAssayFromMatrix <- function(objects, scaled.objects = NULL, features = NULL) {
  if (!check_dimnames(objects)) {
    stop("row and colnames on objects cannot be NULL")
  }

  if (is.null(features)) {
    features <- data.frame(id = colnames(objects))
  }

  if (!is.null(scaled.objects)) {
    if (!check_dimnames(scaled.objects)) {
      stop("row and colnames on scaled.objects cannot be NULL")
    }
    scaled.objects <- TglowMatrix(scaled.objects)
  }

  assay <- new("TglowAssay",
    data = TglowMatrix(objects),
    scale.data = scaled.objects,
    features = features
  )

  if (!tglowr::isValid(assay)) {
    warning("Assay did not pass validity check")
  }

  return(assay)
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
TglowAssayFromList <- function(output, assay, meta.cols, col.object) {
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
#' @param output Tglow list object obtained from [read_cellprofiler_dir()]
#' @param assay Which assay to convert (name of the list item)
#' @param meta.patterns Grep patterns in collumn names asspciated with object level metadata items
#' @param img.feature.patterns Grep patterns in image column names assocated with image level features to analyze
#' @param col.object The collumn name in the features which contains the per object object identifier
#' @param col.img.id The collumn name in the features which contains the per object image identifier
#' @param col.meta.img.id The collumn name in the image level data which contains the image id
#'
#' @returns A populated TglowDataset
#'
#' @export
TglowDatasetFromList <- function(
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
  main.assay <- TglowAssayFromList(output, assay = assay, meta.cols = meta.cols, col.object = col.object)
  dataset <- new("TglowDataset")
  dataset@assays[["raw"]] <- main.assay
  dataset@active.assay <- "raw"
  dataset@object.ids <- output[[assay]][, col.object]
  names(dataset@object.ids) <- dataset@object.ids

  # Filter meta
  dataset@meta <- output[[assay]][, meta.cols]
  rownames(dataset@meta) <- output[[assay]][, col.object]

  # Image level data
  image.cols <- unlist(sapply(img.feature.patterns, grep, colnames(output$meta), value = T))
  dataset@image.data <- new("TglowAssay",
    data = TglowMatrix(as.matrix(output$meta[, image.cols])),
    features = tglowr::get_feature_meta_from_names(image.cols)
  )
  rownames(dataset@image.data@data) <- output$meta[, col.meta.img.id]
  dataset@image.meta <- output$meta[, !colnames(output$meta) %in% image.cols]
  rownames(dataset@image.meta) <- output$meta[, col.meta.img.id]

  dataset@image.ids <- dataset@meta[, col.img.id]
  names(dataset@image.ids) <- dataset@object.ids
  return(dataset)
}

#-------------------------------------------------------------------------------
#' Check if dimnames are set
check_dimnames <- function(object) {
  if (is.null(colnames(object))) {
    warning("colnames of object cannot be NULL")
    return(FALSE)
  }

  if (is.null(rownames(object))) {
    warning("rownames of object cannot be NULL")
    return(FALSE)
  }

  return(TRUE)
}
