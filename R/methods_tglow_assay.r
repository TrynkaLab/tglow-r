#' @include generics.r
NULL

#-------------------------------------------------------------------------------
setMethod(
  "$", signature("TglowAssay"),
  function(x, name) {
    #dx <- x@data@.Data[, name]
    dx <- x@data[, name]
    sx <- NULL

    if (!is.null(x@scale.data)) {
      #sx <- x@scale.data@.Data[, name]
      sx <- x@scale.data[, name]
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

#-------------------------------------------------------------------------------
setMethod(
  "objectIds", signature("TglowAssay"),
  function(object) {
    return(rownames(object@data))
  }
)


setMethod(
  "objectIds<-", signature("TglowAssay"),
  function(object, value) {
    if (length(value) != nrow(object)) {
      stop("New id's must have same length as object")
    }

    if (length(unique(value)) != length(value)) {
      stop("New id's must be unique")
    }

    if (!is.null(object@data)) {
      rownames(object@data) <- value
    }

    if (!is.null(object@scale.data)) {
      rownames(object@scale.data) <- value
    }
    object
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "isValid", signature("TglowAssay"),
  function(object, object.names) {
    if (is.null(object.names)) {
      stop("object.names should be supplied to properly validate the TglowAssay")
    }
    # Check data slot
    if (!isValid(object@data, object.names)) {
      warning("@data slot is not valid")
      return(FALSE)
    }

    # Check scale.data slot
    if (!is.null(object@scale.data)) {
      if (!isValid(object@scale.data, object.names)) {
        warning("@scale.data slot is not valid")
        return(FALSE)
      }
    }

    # Check features slot
    if (is.null(object@features)) {
      warning("Features slot is NULL")
      return(FALSE)
    } else {
      if (is.null(rownames(object@features))) {
        warning("Rownames on @features are not set")
        return(FALSE)
      }

      if (sum(rownames(object@features) %in% colnames(object@data)) != nrow(object@features)) {
        warning("Not all features in @features found in colnames of @data")
        return(FALSE)
      }

      if (sum(colnames(object@data) %in% rownames(object@features)) != ncol(object@data)) {
        warning("Not all features in columns in @data found in rownames of @features")
        return(FALSE)
      }

      if (!is.null(object@scale.data)) {
        if (sum(rownames(object@features) %in% colnames(object@scale.data)) != nrow(object@features)) {
          warning("Not all features in @features found in colnames of @scale.data")
          return(FALSE)
        }
        if (sum(colnames(object@scale.data) %in% rownames(object@features)) != ncol(object@scale.data)) {
          warning("Not all features in columns in @scale.data found in rownames of @features")
          return(FALSE)
        }
      }
    }

    return(TRUE)
  }
)
