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


#-------------------------------------------------------------------------------
setMethod(
  "isValid", signature("TglowReduction"),
  function(object, object.names) {
    if (is.null(object@x)) {
      warning("@x on TglowReduction may not be NULL")
      return(FALSE)
    }

    if (nrow(object@x) != length(object.names)) {
      warning("Number of rows of @x on TglowReduction does not match length of object.names")
      return(FALSE)
    }

    if (is.null(rownames(object@x))) {
      warning("Rownames of @x on TglowReduction are NULL")
      return(FALSE)
    }

    if (sum(rownames(object@x) == object.names) != length(object.names)) {
      warning("Rows of @x on TglowReduction are in a different order then object.names")
      return(FALSE)
    }

    if (!is.null(object@var)) {
      if (!is.numeric(object@var)) {
        warning("@var on TglowReduction must be numeric")
        return(FALSE)
      }

      if (length(object@var) != ncol(object@x)) {
        warning("Lenght of @var on TglowReduction must be equal ncol of @x")
        return(FALSE)
      }
    }

    if (!is.null(object@var_total)) {
      if (!is.numeric(object@var_total)) {
        warning("@var_total on TglowReduction must be numeric")
        return(FALSE)
      }
    }

    return(TRUE)
  }
)
