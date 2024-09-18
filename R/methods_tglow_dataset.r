#' @include generics.r
#'
#'

#-------------------------------------------------------------------------------
setMethod(
  "$", signature("TglowDataset"),
  function(x, name) {
    return(x@assays[[name]])
  }
)

#' @export
.DollarNames.TglowDataset <- function(x, pattern = "") {
  # x is the .CompletionEnv
  return(names(x@assays))
}

#-------------------------------------------------------------------------------
setMethod(
  "[[", signature("TglowDataset"),
  function(x, i) {
    return(x@assays[[i]])
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "show", signature("TglowDataset"),
  function(object) {
    cat(
      "TglowData with: ",
      nrow(object@meta), " objects (cells), ",
      nrow(object@image.meta), " images and ",
      length(object@assays), " assays \n"
    )

    cat("- Assays:\n")
    for (assay in names(object@assays)) {
      cat("\t$", assay, ":", sep = "")
      cat("\t", capture.output(show(object@assays[[assay]])), "\n", sep = "")
    }
    cat("- Active assay:", object@active.assay, "\n")
    cat("- Reductions:", names(object@reduction), "\n")
    cat("- Object size:", format(object.size(object), "Gb", digits = 2))
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
      if (!is.null(object@image.data.trans)) {
        object@image.data.trans <- object@image.data.trans[unique(object@image.ids), , drop = F]
      }
      if (!is.null(object@image.data.trans)) {
        object@image.data.norm <- object@image.data[unique(object@image.ids), , drop = F]
      }

      # Filter dimension reductions
      if (length(object@reduction) >= 1) {
        for (k in seq_along(object@reduction)) {
          object@reduction[[k]] <- object@reduction[[k]][i, , drop = F]
        }
      }

      # Set the graph to NULL to avoid issues
      if (!is.null(object@graph)) {
        warning("Removing graph from object to avoid bugs. Subsetting graph is on TODO list.")
        object@graph <- NULL
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
#' Get image data and features per object (cell).
#'
#' @description Select columns from assay, assay.image, image.meta from 'data' or 'scale.data' slots
#' and return them as a data.frame.
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j Character with column names from assay, assay.image, image.meta to select
#' @param assay The assay to select from
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm"
#' @param slot Slot to fetch features from: "data" or "scale.data"
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setMethod(
  "getDataByObject", signature("TglowDataset"),
  function(object, j, assay, assay.image, slot, drop) {
    # Check the imputs
    check_dataset_assay_slot(object, assay = assay, slot = slot, assay.image = assay.image)

    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) > 0) {
      stop("j can not have NA")
    }

    # Figure out from which table to retieve which features
    j <- as.character(j)

    # Image level features
    is.image <- (j %in% colnames(object@image.meta)) | (j %in% colnames(object@image.data@data))
    j.image  <- j[is.image]

    # Object level features
    j.object <- j[!is.image]
    is.meta <- (j %in% colnames(object@meta))
    j.meta <- j.object[is.meta]
    j.feature <- j.object[!is.meta]

    # Find and get feature data
    data <- NULL
    if (!is.null(assay)) {
      if (!is.null(slot)) {
        if (length(j.feature) >= 1) {
          data <- data.frame(slot(object@assays[[assay]], slot))[, j.feature, drop = F]
        }
      } else {
        stop("Must provide a slot")
      }
    } else {
      if (length(j.feature) >= 1) {
        stop("Object level features but no assay was provided")
      }
    }

    # Find and get metadata
    if (length(j.meta) >= 1) {
      cur.meta <- object@meta[, j.meta, drop = F]

      if (!is.null(data)) {
        data <- cbind(data, cur.meta)
      } else {
        data <- cur.meta
      }
    }

    # Find and get image data
    image.data <- NULL
    if (length(j.image) >= 1) {
      image.data <- data.frame(getImageDataByObject(object, j.image, assay.image, slot = slot, drop = F))
    }

    if (!is.null(data) && !is.null(image.data)) {
      data <- cbind(data, image.data)
    } else if (!is.null(image.data)) {
      data <- image.data
    }

    return(data[, j, drop = drop])
  }
)

#-------------------------------------------------------------------------------
#' Fetch image data or meta data from a tglow object per object (cell).
#'
#' @description Select columns from assay.image, image.meta and return them
#' as a data.frame per object (cell).
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j Character with column names from image.meta or assay.image to select
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm"
#' @param slot Slot to fetch features from data or scale.data
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setMethod(
  "getImageDataByObject", signature("TglowDataset"),
  function(object, j, assay.image, slot, drop) {
    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) > 0) {
      stop("j can not have NA")
    }

    data <- getImageData(object, j, assay.image, slot, drop = F)[object@image.ids, , drop = F]
    rownames(data) <- object@object.ids
    return(data[, j, drop = drop])
  }
)

#-------------------------------------------------------------------------------
#' Fetch image data or meta data from a tglow object
#'
#' @description Select columns from assay.image, image.meta and return them
#' as a data.frame per image.
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j character with column names from image.meta or assay.image to select
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm"
#' @param slot slot to fetch features fromdata or scale.data
#' @param drop should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setMethod(
  "getImageData", signature("TglowDataset"),
  function(object, j, assay.image, slot, drop) {
    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) > 0) {
      stop("j can not have NA")
    }

    is.meta <- j %in% colnames(object@image.meta)

    if (!is.null(assay.image)) {
      is.data <- j %in% colnames(slot(object, assay.image)@data)
    }

    meta <- data.frame(object@image.meta[, j[is.meta], drop = F])

    if (!is.null(assay.image)) {
      data <- data.frame(slot(slot(object, assay.image), slot)[, j[is.data], drop = F])
      output <- cbind(meta, data)
    } else {
      output <- meta
    }

    return(output[, j, drop = drop])
  }
)
