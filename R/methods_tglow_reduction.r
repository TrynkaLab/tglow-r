#-------------------------------------------------------------------------------
setMethod(
    "[",
    "TglowReduction",
    function(x, i, j, drop = F) {
        x@x <- x@x[i, , drop = F]
        return(x)
    }
)
