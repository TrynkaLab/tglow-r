#' @include generics.r
#'
#'
NULL

#-------------------------------------------------------------------------------
setMethod(
  "$", signature("TglowMatrix"),
  function(x, name) {
    return(x[, name])
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
