#' @include generics.r
#'
NULL

#-------------------------------------------------------------------------------
setMethod(
  "$", signature("TglowMatrix"),
  function(x, name) {
    return(as.numeric(x[, name]))
  }
)

#' @export
.DollarNames.TglowMatrix <- function(x, pattern = "") {
  # x is the .CompletionEnv
  return(colnames(x))
}

#-------------------------------------------------------------------------------
# setMethod(
#   "[[", signature("TglowMatrix"),
#   function(x, name) {
#     return(x[, name])
#   }
# )

#-------------------------------------------------------------------------------
setMethod(
  "[",
  "TglowMatrix",
  function(x, i, j, drop = F) {
    
    if (!missing(i)) {
      if (sum(is.na(i)) != 0 ) {
        stop("i has NA values, this is not allowed")
      }
    }

    if (!missing(j)) {
      if (sum(is.na(j)) != 0 ) {
        stop("j has NA values, this is not allowed")
      }
    }
    
    # if (drop) {
    #   if ((!missing(i) && !missing(j)) && (length(i) == 1 || length(j) == 1)) {
    #     return(x@.Data[i, j, drop = drop])
    #   }

    #   if (!missing(i) && length(i) == 1 && missing(j)) {
    #     return(x@.Data[i, , drop = drop])
    #   }

    #   if (!missing(j) && length(j) == 1 && missing(i)) {
    #     return(x@.Data[, j, drop = drop])
    #   }
    # }
    #x@.Data <- x@.Data[i, j, drop = FALSE]

    result <- callNextMethod()
    
    if (is(result, "matrix")) {
      new("TglowMatrix", result)
    } else {
      result
    }  
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "isValid", signature("TglowMatrix"),
  function(object, object.names) {
    if (is.null(object.names)) {
      stop("object.names should be supplied to properly validate the TglowMatrix")
    }

    if (is.null(rownames(object))) {
      warning("Rownames of TglowMatrix are NULL")
      return(FALSE)
    }

    if (is.null(colnames(object))) {
      warning("Colnames of TglowMatrix are NULL")
      return(FALSE)
    }

    if (nrow(object) != length(object.names)) {
      warning("Number of rows on TglowMatrix does not match length of object.names")
      return(FALSE)
    }

    if (sum(rownames(object) == object.names) != length(object.names)) {
      warning("Rows of TglowMatrix are in a different order then object.names")
      return(FALSE)
    }

    return(TRUE)
  }
)
