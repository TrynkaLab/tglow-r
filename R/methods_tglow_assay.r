#' @include generics.r
#-------------------------------------------------------------------------------
setMethod(
  "show", signature("TglowAssay"),
  function(object) {
    cat(
      "TglowAssay with: ",
      nrow(object@data), " cells and ",
      ncol(object@data), " features \n"
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
  "[",
  "TglowAssay",
  function(x, i, j, drop = F) {
    object <- x

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
