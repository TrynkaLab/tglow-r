#-------------------------------------------------------------------------------
# Filters
#-------------------------------------------------------------------------------
# All filters are inclusive!

#-------------------------------------------------------------------------------
#' Filters to apply to TglowDatasets
#'
#' @description Filters come in filter_vec which returns a vector, and filter_agg which takes a vector and retruns a single value.
#' They respectively have modifiers _sum and _multicol which apply them over multicolumn data. See details.
#' 
#' @param vec Input vector or matrix
#' @param thresh Filter threshold value
#' @param grouping Optional grouping for filters applied on a subgroup of vec, not all filters respect this!
#' @returns A logical where TRUE should be kept and FALSE values should be removed
#'
#' @details
#' Filters can be easily configured based on a filter table, making it easy to template sets of operations.
#' Filters are NOT applied seqeuntially, but run independently. If you do want to run filters in seqeuntially,
#'  you will have to run successive iterations, but this is easy enough to do. An easy way to maintain filters
#' and edit them is to store them in a google sheet and load them into R. Then using the function `tglow_filters_from_table`
#' to create the filter objects.
#'
#' The filter table should have the following columns, and one sheet for feature level filters, and one for object level filters.
#' Exact layouts are customizable, see the help of `tglow_filters_from_table`

#' There are two flavors of filters:
#' - filter_vec_x: Accepts a vector and returns a logical vector of the same length (i.e. 'which objects for this feature are > 0')
#' - filter_agg_x: Accepts a vector and aggregates on a statistic and returns a single logical (i.e 'is the variance of this feature > 0')
#'
#' Then there are the filter modifiers
#' - filter_vec_x_sum: Applies the filter to multiple columns, returning a logical of nrow(input), where T only if all columns for that row are T, otherwise F
#' - filter_agg_x_mutlticol: Applies a filter to data with multiple columns and returns a logical vector of ncol(input). If you want to apply these at the object level (i.e. 'filter objects with >x% of NA features'), make sure to set `transpose=T` in the filter definition, if you want to filter features (i.e. 'filter features with >x% of NA objects') leave `transpose=F`.
#'
#' @rdname tglow_filters
#' @export filter_vec_min
filter_vec_min <- function(vec, thresh, grouping = NULL) {
    return(vec >= thresh)
}
#' @rdname tglow_filters
#' @export
filter_vec_min_sum <- function(...) {
    filter_sum(..., func = filter_vec_min)
}

#-------------------------------------------------------------------------------
#' Max filter
#' @rdname tglow_filters
#' @export
filter_vec_max <- function(vec, thresh, grouping = NULL) {
    return(vec <= thresh)
}
#' @rdname tglow_filters
#' @export
filter_vec_max_sum <- function(...) {
    filter_sum(..., func = filter_vec_max)
}

#-------------------------------------------------------------------------------
#' Modified z-score filter
#' @rdname tglow_filters
#' @export
filter_vec_mod_z <- function(vec, thresh, grouping = NULL, absolute = T, method = "mod.z") {
    if (is.null(grouping)) {
        grouping <- rep(1, length(vec))
    }
    mod.z <- tglowr::grouped_scale(vec, grouping = grouping, method = method)

    if (absolute) {
        mod.z <- abs(mod.z)
    }

    return(mod.z < thresh)
}
#' @rdname tglow_filters
#' @export
filter_vec_mod_z_sum <- function(...) {
    filter_sum(..., func = filter_vec_mod_z)
}

#-------------------------------------------------------------------------------
#' @rdname tglow_filters
#' @export
filter_vec_mod_z_perc <- function(vec, thresh, thresh2, grouping = NULL) {
    i <- 0
    j <- ncol(vec)
    cat("\n[INFO] Normalizing features per group.\n")

    res <- apply(vec, 2, function(x) {
        i <<- i + 1
        cat("\r[INFO]", round((i / j) * 100, digits = 2), "%")
        return(filter_vec_mod_z(x, thresh = thresh, grouping = grouping))
    })
    cat("\n[INFO] Done normalizing per group. Calculating percentages per cell\n")
    return((rowSums(res, na.rm = T) / ncol(vec)) >= thresh2)
}

#-------------------------------------------------------------------------------
#' NA filter
#' @rdname tglow_filters
#' @export
filter_agg_na <- function(vec, thresh, grouping = NULL) {
    return((sum(is.na(vec)) / length(vec)) <= thresh)
}
#' @rdname tglow_filters
#' @export
filter_agg_na_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_na)
}

#-------------------------------------------------------------------------------
#' Caret near zero variance filter
#' rdname tglow_filters
#' importFrom caret nearZeroVar
#' export
#filter_agg_near_zero_var <- function(vec, thresh = NULL, grouping = NULL) {
#    return(length(caret::nearZeroVar(vec)) == 0)
#}
#' rdname tglow_filters
#' export
#filter_agg_near_zero_var_multicol <- function(...) {
#    filter_multicol(..., func = filter_agg_near_zero_var)
#}

#-------------------------------------------------------------------------------
#' Zero variance filter
#' @rdname tglow_filters
#' @export
filter_agg_zero_var <- function(vec, thresh = 0) {
    # return(Rfast::Var(vec[!is.na(vec)]) > thresh)
    return(var(vec[!is.na(vec)]) > thresh)
}

#' @rdname tglow_filters
#' @export
filter_agg_zero_var_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_zero_var)
}

#-------------------------------------------------------------------------------
#' Absolute coefficient of variation filter
#' @rdname tglow_filters
#' @export
filter_agg_coef_var <- function(vec, thresh) {
    # cur.var <- var(vec[!is.na(vec)])
    # cur.var[cur.var < 1e-10] <- 0
    return(abs((sd(vec[!is.na(vec)]) / mean(vec[!is.na(vec)]))) > thresh)
}

#' @rdname tglow_filters
#' @export
filter_agg_coef_var_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_coef_var)
}
#-------------------------------------------------------------------------------
#' Absolute skewness filter
#' @rdname tglow_filters
#' @export
filter_agg_skewness <- function(vec, thresh) {
    # cur.var <- var(vec[!is.na(vec)])
    # cur.var[cur.var < 1e-10] <- 0
    return(abs(skewness(vec, na.rm=T)) > thresh)
}

#' @rdname tglow_filters
#' @export
filter_agg_skewness_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_skewness)
}

#-------------------------------------------------------------------------------
#' Absolute kurotsis filter
#' @rdname tglow_filters
#' @export
filter_agg_kurtosis <- function(vec, thresh) {
    # cur.var <- var(vec[!is.na(vec)])
    # cur.var[cur.var < 1e-10] <- 0
    return(abs(kurtosis(vec, na.rm=T)) > thresh)
}

#' @rdname tglow_filters
#' @export
filter_agg_kurtosis_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_skewness)
}

#-------------------------------------------------------------------------------
#' Minimum number of unique values
#' @rdname tglow_filters
#' @export
filter_agg_unique_val <- function(vec, thresh = NULL) {
    return(length(unique(vec)) > thresh)
}
#' @rdname tglow_filters
#' @export
filter_agg_unique_val_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_unique_val)
}

#-------------------------------------------------------------------------------
#' Infinite median filter
#' @rdname tglow_filters
#' @export
filter_agg_inf_median <- function(vec, thresh = NULL) {
    return(!is.infinite(median(vec, na.rm = T)))
}
#' @rdname tglow_filters
#' @export
filter_agg_inf_median_sum <- function(...) {
    filter_sum(..., func = filter_agg_inf_median)
}

#-------------------------------------------------------------------------------
#' Infinite sum filter
#' @rdname tglow_filters
#' @export
filter_agg_inf <- function(vec, thresh = NULL, grouping = NULL) {
    return(sum(is.infinite(vec)) <= thresh)
}

#' @rdname tglow_filters
#' @export
filter_agg_inf_multicol <- function(...) {
    filter_multicol(..., func = filter_agg_inf)
}

#-------------------------------------------------------------------------------
#' Blacklist filter, always returns false
#' @rdname tglow_filters
#' @export
filter_agg_blacklist <- function(vec, thresh, grouping = NULL) {
    return(FALSE)
}

#-------------------------------------------------------------------------------
#' If data is multicolumn, take the sum over all collumns, if one is false, exclude
#' @rdname tglow_filters
#' @export
filter_sum <- function(vec, thresh, grouping, func) {
    if (is.null(ncol(vec))) {
        stop("Data must be a matrix or data.frame with mutliple columns")
    }

    res <- apply(vec, 2, func, thresh = thresh, grouping = grouping)

    res[is.na(res)] <- T

    if (!is(res, "matrix")) {
        res <- matrix(res, ncol = 1)
    }

    return(rowSums(res, na.rm = T) == ncol(vec))
}

#-------------------------------------------------------------------------------
#' If data is multicolumn, apply filter over all columns
#' @rdname tglow_filters
#' @export
filter_multicol <- function(vec, thresh, grouping, func) {
    if (is.null(ncol(vec))) {
        stop("Data must be a matrix or data.frame with mutliple columns")
    }

    res <- apply(vec, 2, func, thresh = thresh, grouping = grouping)
    return(res)
}
