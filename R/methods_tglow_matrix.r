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
setMethod(
  "[[", signature("TglowMatrix"),
  function(x, name) {
    return(x[, name])
  }
)

#-------------------------------------------------------------------------------
setMethod(
  "[",
  "TglowMatrix",
  function(x, i, j, drop = F) {
    
    if (drop) {
      if ((!missing(i) && !missing(j)) && (length(i) == 1 || length(j) == 1)) {
        return(x@.Data[i, j, drop = drop])
      }

      if (!missing(i) && length(i) == 1 && missing(j)) {
        return(x@.Data[i, , drop = drop])
      }

      if (!missing(j) && length(j) == 1 && missing(i)) {
        return(x@.Data[, j, drop = drop])
      }
    }

    x@.Data <- x@.Data[i, j, drop = FALSE]
    return(x)
  }
)
