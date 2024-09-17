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
    return(TglowMatrix(x@.Data[i, j, drop = drop]))
  }
)
