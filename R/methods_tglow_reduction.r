#-------------------------------------------------------------------------------
setMethod(
    "[",
    "TglowReduction",
    function(x, i, j, drop = F) {
        x@x <- x@x[i, , drop = F]
        return(x)
    }
)

#-------------------------------------------------------------------------------
setMethod(
  "objectIds", signature("TglowReduction"),
  function(object) {
    return(rownames(object@data))
  }
)


setMethod(
  "objectIds<-", signature("TglowReduction"),
  function(object, value) {
    if (length(value) != nrow(object@x)) {
      stop("New id's must have same length as object")
    }

    if (length(unique(value)) != length(value)) {
      stop("New id's must be unique")
    }

    if (!is.null(object@x)) {
      rownames(object@x) <- value
    }
    object
  }
)
