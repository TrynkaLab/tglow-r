#-------------------------------------------------------------------------------
# Filters
#-------------------------------------------------------------------------------
# All filters are inclusive!

#-------------------------------------------------------------------------------
#' Filters to apply to TglowDatasets
#'
#' @description Filters come in 2 variants, the base function, which is applied over a vector or value
#' The base then has a <filter>.sum function which sums over multiple vectors or values
#'
#' @details
#' Uses Rfast whose variance calculations are not exactly zero, so set the zero threshold to 1e-10 by default
#'
#' @param vec Input vector or value
#' @param thresh Filter threshold value
#' @param grouping Optional grouping for filters applied on a subgroup of vec, not all filters respect this!
#' @returns A logical where TRUE should be kept and FALSE values should be removed
#'
#' @rdname tglow_filters
#' @export filter_min
filter_min <- function(vec, thresh, grouping = NULL) {
    return(vec >= thresh)
}
#' @rdname tglow_filters
#' @export
filter_min_sum <- function(...) {
    filter_sum(..., func = filter_min)
}

#-------------------------------------------------------------------------------
#' Max filter
#' @rdname tglow_filters
#' @export
filter_max <- function(vec, thresh, grouping = NULL) {
    return(vec <= thresh)
}
#' @rdname tglow_filters
#' @export
filter_max_sum <- function(...) {
    filter_sum(..., func = filter_max)
}

#-------------------------------------------------------------------------------
#' NA filter
#' @rdname tglow_filters
#' @export
filter_na <- function(vec, thresh, grouping = NULL) {
    return((sum(is.na(vec)) / length(vec)) <= thresh)
}
#' @rdname tglow_filters
#' @export
filter_na_multicol <- function(vec, thresh, grouping) {
    res <- apply(vec, 2, filter_na, thresh = thresh, grouping = grouping)
    return(res)
}

#-------------------------------------------------------------------------------
#' Caret near zero variance filter
#' @rdname tglow_filters
#' @export
filter_near_zero_var <- function(vec, thresh = NULL) {
    return(length(nearZeroVar(vec)) == 0)
}
#' @rdname tglow_filters
#' @export
filter_near_zero_var_sum <- function(...) {
    filter_sum(..., func = filter_near_zero_var)
}

#-------------------------------------------------------------------------------
#' Zero variance filter
#' @rdname tglow_filters
#' @export
filter_zero_var <- function(vec, thresh = 0) {
    # return(Rfast::Var(vec[!is.na(vec)]) > thresh)
    return(var(vec[!is.na(vec)]) > thresh)
}

#' @rdname tglow_filters
#' @export
filter_zero_var_sum <- function(...) {
    filter_sum(..., func = filter_zero_var)
}

#-------------------------------------------------------------------------------
#' Absolute coefficient of variation filter
#' @rdname tglow_filters
#' @export
filter_coef_var <- function(vec, thresh) {
    #cur.var <- var(vec[!is.na(vec)])
    # cur.var[cur.var < 1e-10] <- 0
    return(abs((sd(vec[!is.na(vec)]) / mean(vec[!is.na(vec)]))) > thresh)
}

#' @rdname tglow_filters
#' @export
filter_coef_var_sum <- function(...) {
    filter_sum(..., func = filter_coef_var)
}

#-------------------------------------------------------------------------------
#' Minimum number of unique values
#' @rdname tglow_filters
#' @export
filter_unique_val <- function(vec, thresh = NULL) {
    return(length(unique(vec)) > thresh)
}
#' @rdname tglow_filters
#' @export
filter_unique_val_sum <- function(...) {
    filter_sum(..., func = filter_unique_val)
}

#-------------------------------------------------------------------------------
#' Infinite median filter
#' @rdname tglow_filters
#' @export
filter_inf_median <- function(vec, thresh = NULL) {
    return(!is.infinite(median(vec, na.rm = T)))
}
#' @rdname tglow_filters
#' @export
filter_inf_median_sum <- function(...) {
    filter_sum(..., func = filter_inf_median)
}

#-------------------------------------------------------------------------------
#' Infinite sum filter
#' @rdname tglow_filters
#' @export
filter_inf <- function(vec, thresh = NULL, grouping = NULL) {
    return(sum(is.infinite(vec)) <= thresh)
}
#' @rdname tglow_filters
#' @export
filter_inf_mutlicol <- function(vec, thresh, grouping) {
    res <- apply(vec, 2, filter_inf, thresh = thresh, grouping = grouping)
    return(res)
}

#-------------------------------------------------------------------------------
#' Modified z-score filter
#' @rdname tglow_filters
#' @export
filter_mod_z <- function(vec, thresh, grouping = NULL, absolute = T, method = "mod.z") {
    if (is.null(grouping)) {
        grouping <- rep(1, length(vec))
    }
    mod.z <- grouped_scale(vec, grouping = grouping, method = method)

    if (absolute) {
        mod.z <- abs(mod.z)
    }

    return(mod.z < thresh)
}
#' @rdname tglow_filters
#' @export
filter_mod_z_sum <- function(...) {
    filter_sum(..., func = filter_mod_z)
}

#-------------------------------------------------------------------------------
#' @rdname tglow_filters
#' @export
filter_mod_z_perc <- function(vec, thresh, thresh2, grouping = NULL) {
    i <- 0
    j <- ncol(vec)
    cat("\n[INFO] Normalizing features per group.\n")

    res <- apply(vec, 2, function(x) {
        i <<- i + 1
        cat("\r[INFO]", round((i / j) * 100, digits = 2), "%")
        return(filter_mod_z(x, thresh = thresh, grouping = grouping))
    })
    cat("\n[INFO] Done normalizing per group. Calculating percentages per cell\n")
    return((rowSums(res, na.rm = T) / ncol(vec)) >= thresh2)
}

#-------------------------------------------------------------------------------
#' If data is multicolumn, take the sum over all collumns, if one is false, exlcude
#' @rdname tglow_filters
#' @export
filter_sum <- function(vec, thresh, grouping, func) {
    
    if (is.null(ncol(vec))) {
        stop("Data must be a matrix or data.frame with mutlipel columns")
    }
    
    res <- apply(vec, 2, func, thresh = thresh, grouping = grouping)

    res[is.na(res)] <- T

    return(rowSums(res, na.rm = T) == ncol(vec))
}


#-------------------------------------------------------------------------------
#' Blacklist filter, always returns false
#' @rdname tglow_filters
#' @export
filter_blacklist <- function(vec, thresh, grouping = NULL) {
    return(FALSE)
}
