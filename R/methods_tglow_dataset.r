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
  function(x, i, j, drop = F, na.check = T) {
    object <- x


    if (!missing(i) && na.check) {
      if (sum(is.na(i)) != 0 ) {
        stop("i has NA values, this is not allowed")
      }
    }

    if (!missing(j) && na.check) {
      if (sum(is.na(j)) != 0 ) {
        stop("j has NA values, this is not allowed")
      }
    }

    # Select rows
    if (!missing(i)) {
      # Filter assays
      for (assay in seq_along(object@assays)) {
        object@assays[[assay]] <- object@assays[[assay]][i, , drop = F, na.check = na.check]
      }

      # Filter meta
      object@meta <- object@meta[i, , drop = F]

      if (is.null(names(object@object.ids))) {
        names(object@object.ids) <- object@object.ids
      }

      object@object.ids <- object@object.ids[i]

      # Select images
      object@image.ids <- object@image.ids[i, drop = F]
      object@image.data <- object@image.data[unique(object@image.ids), , drop = F, na.check = na.check]
      object@image.meta <- object@image.meta[unique(object@image.ids), , drop = F]

      if (!is.null(object@image.data.trans)) {
        object@image.data.trans <- object@image.data.trans[unique(object@image.ids), , drop = F, na.check = na.check]
      }
      if (!is.null(object@image.data.norm)) {
        object@image.data.norm <- object@image.data.norm[unique(object@image.ids), , drop = F, na.check = na.check]
      }

      # Filter dimension reductions
      if (length(object@reduction) >= 1) {
        for (k in seq_along(object@reduction)) {
          object@reduction[[k]] <- object@reduction[[k]][i, , drop = F, na.check = na.check]
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
        object@assays[[assay]] <- object@assays[[assay]][, j, drop = F, na.check = na.check]
      }
    }

    object
  }
)


#-------------------------------------------------------------------------------
setMethod(
  "objectIds", signature("TglowDataset"),
  function(object) {
    return(object@object.ids)
  }
)

setMethod(
  "objectIds<-", signature("TglowDataset"),
  function(object, value) {
    if (length(value) != nrow(object)) {
      stop("New id's must have same length as object")
    }

    if (length(unique(value)) != length(value)) {
      stop("New id's must be unique")
    }

    object@object.ids <- value
    names(object@object.ids) <- object@object.ids
    
    names(object@image.ids) <- value
    rownames(object@meta) <- value

    for (assay in names(object@assays)) {
      objectIds(object@assays[[assay]]) <- value
    }

    for (reduction in names(object@reduction)) {
      objectIds(object@reduction[[reduction]]) <- value
    }

    object
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "isAvailable", signature("TglowDataset"),
  function(object, j, assay, assay.image, slot, return.names) {
    # Check the inputs
    check_dataset_assay_slot(object, assay = assay, slot = slot, assay.image = assay.image)

    if (class(j) != "character") {
      stop("j must be character vector with column names in meta or assay.")
    }
    if (sum(is.na(j)) != 0 ) {
      stop("j has NA values, this is not allowed")
    }
    # Image level features
    if (!is.null(assay.image)) {
      is.image <- j %in% colnames(slot(object, assay.image))
    } else {
      is.image <- rep(TRUE, length(j))
    }

    is.image.meta <- j %in% colnames(object@image.meta)
    is.meta <- j %in% colnames(object@meta)

    if (!is.null(assay)) {
      is.assay <- j %in% colnames(object[[assay]])
    } else {
      is.assay <- rep(TRUE, length(j))
    }

    exists <- is.meta | is.image.meta | is.meta | is.assay

    if (is.null(assay.image)) {
      is.image <- rep(NA, length(is.image))
    }

    if (is.null(assay)) {
      is.assay <- rep(NA, length(is.assay))
    }

    exists.sum <- rowSums(cbind(is.image, is.image.meta, is.meta, is.assay), na.rm = T)

    if (any(exists.sum > 1)) {
      warning("Collumn names are not unique accross image.data, image.meta, meta, assay")
    }

    if (!return.names) {
      return(sum(exists) == length(j))
    } else {
      return(list(
        is.available = sum(exists) == length(j),
        exists = exists,
        exists.sum = exists.sum,
        is.in.slot = list(image.data = is.image, image.meta = is.image.meta, meta = is.meta, assay = is.assay)
      ))
    }
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "getDataByObject", signature("TglowDataset"),
  function(object, j, assay, assay.image, slot, drop) {
    # Check the inputs
    if (!isAvailable(object, j, assay, assay.image, slot, return.names = FALSE)) {
      stop("One or more items in j was not available. You can check which with isAvailable(...,return.names=T)")
    }

    # Image level features
    is.image <- (j %in% colnames(object@image.meta)) | (j %in% colnames(object@image.data@data))
    j.image <- j[is.image]

    # Object level features
    j.object <- j[!is.image]
    is.meta <- (j.object %in% colnames(object@meta))
    j.meta <- j.object[is.meta]
    j.feature <- j.object[!is.meta]

    # Find and get feature data
    data <- NULL
    if (!is.null(assay)) {
      if (!is.null(slot)) {
        if (length(j.feature) >= 1) {
          #data <- data.frame(slot(object@assays[[assay]], slot)@.Data[, j.feature, drop = F])
          data <- data.frame(slot(object@assays[[assay]], slot)[, j.feature, drop = F])
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
setMethod(
  "getImageData", signature("TglowDataset"),
  function(object, j, assay.image, slot, drop) {
    # Check the inputs
    if (!isAvailable(object, j, assay = NULL, assay.image = assay.image, slot, return.names = FALSE)) {
      stop("One or more items in j was not available. You can check which with isAvailable(...,return.names=T)")
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

#-------------------------------------------------------------------------------
setMethod(
  "isValid", signature("TglowDataset"),
  function(object, object.names) {
    #-------------------------------------
    # Check object IDs
    #-------------------------------------
    if (is.null(object@object.ids) || !is.character(object@object.ids)) {
      warning("@object.ids on TglowDataset may not be NULL and must be character")
      return(FALSE)
    }

    if (is.null(object@image.ids) || !is.character(object@image.ids)) {
      warning("@image.ids on TglowDataset may not be NULL and must be character")
      return(FALSE)
    }

    # If object id's not provided, use the slot
    if (is.null(object.names)) {
      object.names <- object@object.ids
    }

    #-------------------------------------
    # Check object and image ids
    #-------------------------------------
    if (length(object.names) != length(object@object.ids)) {
      warning("object.names and @object.ids on TglowDataset are not the same length")
      return(FALSE)
    }

    if (length(object.names) != length(object@image.ids)) {
      warning("object.names and @image.ids on TglowDataset are not the same length")
      return(FALSE)
    }

    if (sum(object@object.ids == object.names) != length(object.names)) {
      warning("@object.ids on TglowDataset are in a different order then provided object.names")
      return(FALSE)
    }
    
    if (is.null(names(object@image.ids))) {
      warning("names of @image.ids on TglowDataset are NULL. Affects slicing with strings. Run objectIds(object) <- object@object.ids to fix")
      return(FALSE)
    }
    
    if (is.null(names(object@object.ids))) {
      warning("names of @object.ids on TglowDataset are NULL. Affects slicing with strings. Run objectIds(object) <- object@object.ids to fix")
      return(FALSE)
    }

    #-------------------------------------
    # Check metadata
    #-------------------------------------
    if (is.null(object@meta)) {
      warning("@meta on TglowDataset cannot be NULL")
      return(FALSE)
    }

    if (!is(object@meta, "data.frame")) {
      warning("@meta on TglowDataset must be data.frame")
      return(FALSE)
    }

    if (is.null(rownames(object@meta))) {
      warning("@meta on TglowDataset must have rownames")
      return(FALSE)
    }

    if (sum(rownames(object@meta) == object.names) != length(object.names)) {
      warning("Rows of @meta on TglowDataset are in a different order then provided object.names")
      return(FALSE)
    }

    #-------------------------------------
    # Check image meta data
    #-------------------------------------
    if (is.null(object@image.meta)) {
      warning("@image.meta on TglowDataset on TglowDataset cannot be NULL")
      return(FALSE)
    }

    if (!is(object@image.meta, "data.frame")) {
      warning("@image.meta on TglowDataset must be data.frame")
      return(FALSE)
    }

    if (is.null(rownames(object@image.meta))) {
      warning("@image.meta on TglowDataset must have rownames")
      return(FALSE)
    }

    if (sum(rownames(object@image.meta) %in% object@image.ids) != nrow(object@image.meta)) {
      warning("Not all rows in @image.meta of TglowDataset found in @image.ids")
      return(FALSE)
    }

    if (sum(unique(object@image.ids) %in% rownames(object@image.meta)) != nrow(object@image.meta)) {
      warning("Not all unqiue values of @image.ids of TglowDataset found in rownames of @image.meta")
      return(FALSE)
    }

    #-------------------------------------
    # Check image data
    #-------------------------------------
    if (is.null(object@image.data)) {
      warning("@image.data on TglowDataset must not be NULL")
      return(FALSE)
    }

    if (!isValid(object@image.data, rownames(object@image.meta))) {
      warning("@image.data on TglowDataset is invalid")
      return(FALSE)
    }

    if (!is.null(object@image.data.norm)) {
      if (!isValid(object@image.data.norm, rownames(object@image.meta))) {
        warning("@image.data.norm on TglowDataset is invalid")
        return(FALSE)
      }
    }

    if (!is.null(object@image.data.trans)) {
      if (!isValid(object@image.data.trans, rownames(object@image.meta))) {
        warning("@image.data.trans on TglowDataset is invalid")
        return(FALSE)
      }
    }

    #-------------------------------------
    # Check assays
    #-------------------------------------
    if (!is(object@assays, "list")) {
      warning("@assays on TglowDataset must be list")
      return(FALSE)
    }

    if (is.null(names(object@assays))) {
      warning("@assays on TglowDataset must have names set. Use names(object@assays) <- c('names', 'go', 'here') to fix")
      return(FALSE)
    }

    if (length(object@assays) < 1) {
      warning("@assays on TglowDataset must have at least one assay")
      return(FALSE)
    }

    for (assay in names(object@assays)) {
      if (!is(object@assays[[assay]], "TglowAssay")) {
        warning(paste0("@assays$", assay, " on TglowDataset is not a TglowAssay"))
        return(FALSE)
      }

      if (!isValid(object@assays[[assay]], object@object.ids)) {
        warning(paste0("@assays$", assay, " on TglowDataset is invalid"))
        return(FALSE)
      }
    }

    #-------------------------------------
    # Check reductions
    #-------------------------------------
    if (!is(object@reduction, "list")) {
      warning("@reduction on TglowDataset must be list")
      return(FALSE)
    }

    if (length(object@reduction) >= 1) {
      if (is.null(names(object@reduction))) {
        warning("@reduction on TglowDataset must have names set. Use names(object@reduction) <- c('names', 'go', 'here') to fix")
        return(FALSE)
      }

      for (reduction in names(object@reduction)) {
        if (!is(object@reduction[[reduction]], "TglowReduction")) {
          warning(paste0("@reduction$", reduction, " on TglowDataset is not a TglowReduction"))
          return(FALSE)
        }

        if (!isValid(object@reduction[[reduction]], object@object.ids)) {
          warning(paste0("@reduction$", reduction, " on TglowDataset is invalid"))
          return(FALSE)
        }
      }
    }


    return(TRUE)
  }
)
