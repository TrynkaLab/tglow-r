#-------------------------------------------------------------------------------
# Filters
#-------------------------------------------------------------------------------
# All filters are inclusive!

#-------------------------------------------------------------------------------
#' Filters to apply to TglowDatasets
#'
#' @description Filters come in 2 variants, the base function, which is applied over a vector or value.
#' The base then has a <filter>.sum function which sums over multiple vectors or values.
#'
#' @param vec Input vector or value
#' @param thresh Filter threshold value
#' @param grouping Optional grouping for filters applied on a subgroup of vec, not all filters respect this!
#' @returns A logical where TRUE should be kept and FALSE values should be removed.
#'
#' @rdname tglow_filters
#' @export f.min
f.min <- function(vec, thresh, grouping = NULL) {
    return(vec >= thresh)
}
#' @rdname tglow_filters
#' @export
f.min.sum <- function(...) {
    f.sum(..., func = f.min)
}

#-------------------------------------------------------------------------------
#' Max filter
#' @rdname tglow_filters
#' @export
f.max <- function(vec, thresh, grouping = NULL) {
    return(vec <= thresh)
}
#' @rdname tglow_filters
#' @export
f.max.sum <- function(...) {
    f.sum(..., func = f.max)
}

#-------------------------------------------------------------------------------
#' NA filter
#' @rdname tglow_filters
#' @export
f.na <- function(vec, thresh, grouping = NULL) {
    return((sum(is.na(vec)) / length(vec)) <= thresh)
}
#' @rdname tglow_filters
#' @export
f.na.multicol <- function(vec, thresh, grouping) {
    res <- apply(vec, 2, f.na, thresh = thresh, grouping = grouping)
    return(res)
}

#-------------------------------------------------------------------------------
#' Caret near zero variance filter
#' @rdname tglow_filters
#' @export
f.near.zero.var <- function(vec, thresh = NULL) {
    return(length(nearZeroVar(vec)) == 0)
}
#' @rdname tglow_filters
#' @export
f.near.zero.var.sum <- function(...) {
    f.sum(..., func = f.near.zero.var)
}

#-------------------------------------------------------------------------------
#' Zero variance filter
#' @rdname tglow_filters
#' @export
f.zero.var <- function(vec, thresh = 0) {
    return(Rfast::Var(vec[!is.na(vec)]) > thresh)
}

#' @rdname tglow_filters
#' @export
f.zero.var.sum <- function(...) {
    f.sum(..., func = f.zero.var)
}

#-------------------------------------------------------------------------------
#' Coefficient of variation filter
#' @rdname tglow_filters
#' @export
f.coef.var <- function(vec, thresh) {
    return((sqrt(Rfast::Var(vec[!is.na(vec)])) / mean(vec[!is.na(vec)])) > thresh)
}
#' @rdname tglow_filters
#' @export
f.coef.var.sum <- function(...) {
    f.sum(..., func = f.zero.var)
}

#-------------------------------------------------------------------------------
#' Minimum number of unique values
#' @rdname tglow_filters
#' @export
f.unique.val <- function(vec, thresh = NULL) {
    return(length(unique(vec)) > thresh)
}
#' @rdname tglow_filters
#' @export
f.unique.val.sum <- function(...) {
    f.sum(..., func = f.unique.val)
}

#-------------------------------------------------------------------------------
#' Infinite median filter
#' @rdname tglow_filters
#' @export
f.inf.median <- function(vec, thresh = NULL) {
    return(!is.infinite(median(vec, na.rm = T)))
}
#' @rdname tglow_filters
#' @export
f.inf.median.sum <- function(...) {
    f.sum(..., func = f.inf.median)
}

#-------------------------------------------------------------------------------
#' Infinite sum filter
#' @rdname tglow_filters
#' @export
f.inf <- function(vec, thresh = NULL, grouping = NULL) {
    return(sum(is.infinite(vec)) <= thresh)
}
#' @rdname tglow_filters
#' @export
f.inf.mutlicol <- function(vec, thresh, grouping) {
    res <- apply(vec, 2, f.inf, thresh = thresh, grouping = grouping)
    return(res)
}

#-------------------------------------------------------------------------------
#' Modified z-score filter
#' @rdname tglow_filters
#' @export
f.mod.z <- function(vec, thresh, grouping = NULL, absolute = T, method = "mod.z") {
    if (is.null(grouping)) {
        grouping <- rep(1, length(vec))
    }
    mod.z <- grouped.scale(vec, grouping = grouping, method = method)

    if (absolute) {
        mod.z <- abs(mod.z)
    }

    return(mod.z < thresh)
}
#' @rdname tglow_filters
#' @export
Met <- function(...) {
    f.sum(..., func = f.mod.z)
}
#' @rdname tglow_filters
#' @export
f.mod.z.perc <- function(vec, thresh, thresh2, grouping = NULL) {
    i <- 0
    j <- ncol(vec)
    cat("\n[INFO] Normalizing features per group.\n")

    res <- apply(vec, 2, function(x) {
        i <<- i + 1
        cat("\r[INFO]", round((i / j) * 100, digits = 2), "%")
        return(f.mod.z(x, thresh = thresh, grouping = grouping))
    })
    cat("\n[INFO] Done normalizing per group. Calculating percentages per cell\n")
    return((rowSums(res, na.rm = T) / ncol(vec)) >= thresh2)
}

#-------------------------------------------------------------------------------
#' If data is multicolumn, take the sum over all collumns, if one is false, exlcude
#' @rdname tglow_filters
#' @export
f.sum <- function(vec, thresh, grouping, func) {
    res <- apply(vec, 2, func, thresh = thresh, grouping = grouping)

    res[is.na(res)] <- T

    return(rowSums(res, na.rm = T) == ncol(vec))
}
