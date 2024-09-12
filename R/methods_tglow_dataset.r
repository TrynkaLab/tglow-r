#' @include generics.r
#-------------------------------------------------------------------------------
setMethod(
  "show", signature("TglowDataset"),
  function(object) {
    cat(
      "TglowData with: ",
      nrow(object@meta), " cells, ",
      nrow(object@image.meta), " images and ",
      length(object@assays), " assays \n"
    )

    cat("Assays: \n")
    print(object@assays)

    cat("\nActive assay: ", object@active.assay, "\n")
    cat("\nReductions: ", names(object@reduction), "\n")
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "nrow", signature("TglowDataset"),
  function(x) {
    length(x@object.ids)
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "[",
  "TglowDataset",
  function(x, i, j, drop = F) {
    object <- x

    # Select rows
    if (!missing(i)) {
      # Filter assays
      for (assay in seq_along(object@assays)) {
        object@assays[[assay]] <- object@assays[[assay]][i, , drop = F]
      }

      # Filter meta
      object@meta <- object@meta[i, , drop = F]
      object@object.ids <- object@object.ids[i, drop = F]

      # Select images
      object@image.ids <- object@image.ids[i, drop = F]
      object@image.data <- object@image.data[unique(object@image.ids), , drop = F]
      object@image.meta <- object@image.meta[unique(object@image.ids), , drop = F]

      # Filter dimension reductions
      if (length(object@reduction) >= 1) {
        for (k in seq_along(object@reduction)) {
          object@reduction[[k]] <- object@reduction[[k]][i, , drop = F]
        }
      }
    }

    # Select columns
    if (!missing(j)) {
      if (class(j) != "character") {
        warning("Assuming all assays have the same column order")
      }

      # Filter assays
      for (assay in seq_along(object@assays)) {
        object@assays[[assay]] <- object@assays[[assay]][, j, drop = F]
      }
    }

    object
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "getImageDataAndFeatures", signature("TglowDataset"),
  function(object, j, assay, slot, drop) {
    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) > 0) {
      stop("j can not have NA")
    }

    is.meta <- (j %in% colnames(object@image.meta)) | (j %in% colnames(object@image.data@data))
    j.meta <- j[is.meta]
    j.feature <- j[!is.meta]

    data <- data.frame(slot(object@assays[[assay]], slot))[, j.feature, drop = F]

    if (sum(is.meta) > 0) {
      meta <- data.frame(getImageDataByCell(object, j.meta, slot = slot, drop = F))
      data <- cbind(data, meta)
    }

    return(data[, j, drop = drop])
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "getImageDataByCell", signature("TglowDataset"),
  function(object, j, slot, drop) {
    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) > 0) {
      stop("j can not have NA")
    }

    data <- getImageData(object, j, slot, drop = F)[object@image.ids, , drop = F]
    rownames(data) <- object@object.ids
    return(data[, j, drop = drop])
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "getImageData", signature("TglowDataset"),
  function(object, j, slot, drop) {
    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) > 0) {
      stop("j can not have NA")
    }

    is.meta <- j %in% colnames(object@image.meta)
    is.data <- j %in% colnames(object@image.data@data)

    meta <- data.frame(object@image.meta[, j[is.meta], drop = F])
    data <- data.frame(slot(object@image.data, slot)[, j[is.data], drop = F])

    output <- cbind(meta, data)
    return(output[, j, drop = drop])
  }
)
