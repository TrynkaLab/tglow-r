#' @include generics.r
NULL

#-------------------------------------------------------------------------------
setMethod(
  "$", signature("TglowAssay"),
  function(x, name) {
    dx <- x@data@.Data[, name]
    sx <- NULL

    if (!is.null(x@scale.data)) {
      sx <- x@scale.data@.Data[, name]
    }
    return(list(data = dx, scale.data = sx))
  }
)

#' @export
.DollarNames.TglowAssay <- function(x, pattern = "") {
  # x is the .CompletionEnv
  return(x@features$id)
}

#-------------------------------------------------------------------------------
setMethod(
  "[[", signature("TglowAssay"),
  function(x, i) {
    dx <- x@data[, i]
    sx <- NULL

    if (!is.null(x@scale.data)) {
      sx <- x@scale.data[, i]
    }
    return(list(data = dx, scale.data = sx))
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "show", signature("TglowAssay"),
  function(object) {
    cat(
      "TglowAssay with:",
      nrow(object@data), "objects and",
      ncol(object@data), "features.",
      "Size:", format(object.size(object), "Gb", digits = 2), "\n"
    )
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "nrow", signature("TglowAssay"),
  function(x) {
    nrow(x@data)
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "ncol", signature("TglowAssay"),
  function(x) {
    ncol(x@data)
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "colnames", signature("TglowAssay"),
  function(x) {
    colnames(x@data)
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "rownames", signature("TglowAssay"),
  function(x) {
    rownames(x@data)
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "[",
  "TglowAssay",
  function(x, i, j, drop = F) {
    object <- x

    if (is.null(object@data)) {
      return(object)
    }

    if (nrow(object) == 0 || ncol(object) == 0) {
      return(object)
    }


    # Select rows
    if (!missing(i)) {
      # Filter main assay
      object@data <- object@data[i, , drop = F]

      # Filter scale data
      if (!is.null(object@scale.data)) {
        if (nrow(object@scale.data) >= 1) {
          object@scale.data <- object@scale.data[i, , drop = F]
        }
      }
    }

    # Select columns
    if (!missing(j)) {
      # Filter main assay
      object@data <- object@data[, j, drop = F]

      # Filter scale data
      if (!is.null(object@scale.data)) {
        if (nrow(object@scale.data) >= 1) {
          object@scale.data <- object@scale.data[, j, drop = F]
        }
      }

      # Filter features
      object@features <- object@features[j, , drop = F]
    }

    object
  }
)
