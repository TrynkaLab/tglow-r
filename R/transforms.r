#-------------------------------------------------------------------------------
#' Faster alternative to scale a matrix
#'
#' @description
#' Adapted from https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' Yields a 10-20% speedup and adds modfied z-score scaling
#'
#' @param x A matrix
#' @param center Should data be centered
#' @param scale Should data be scaled
#' @param add_attr Should means and variances be added to output
#' @param na.rm Should NA's be removed during mean and variance calculations
#' @param scale.method Scaling method, either `mean` for z-score or `median` modified z-score
#' @param reference.group Calculate means/medians, sds/mads only on these samples but apply it to all
#'
#' @returns The scaled matrix
#'
#' @examples
#' # Generate a matrix
#' r <- matrix(runif(100000 * 1000), nrow = 100000, ncol = 1000)
#'
#' # Run using scale
#' system.time({
#'     r.scale <- scale(r)
#' })
#'
#' # Run using fast_colscale
#' system.time({
#'     r.fast_colscale <- fast_colscale(r, add_attr = F)
#' })
#'
#' # Compare results
#' diff <- abs(r.scale - r.fast_colscale)
#' all(diff < 1e-15)
#' @importFrom matrixStats colMedians colSds colMads
#' @export
fast_colscale <- function(x,
                          center = TRUE,
                          scale = TRUE,
                          add_attr = TRUE,
                          na.rm = F,
                          scale.method = "mean",
                          reference.group = NULL) {

    if (is(x, "TglowMatrix")) {
        x <- x@.Data
    }                            
                            
    if (is.null(reference.group)) {
        x.ref <- x
    } else {
        x.ref <- x[reference.group, ]
    }

    if (scale.method == "mean") {
        # Get the column means
        cm <- colMeans(x.ref, na.rm = na.rm)
    } else if (scale.method == "median") {
        cm <- matrixStats::colMedians(x.ref, na.rm = na.rm)
    } else {
        stop("Method must be `mean` or `median`")
    }

    # Get the column sd or mad
    if (scale) {
        if (scale.method == "mean") {
            csd <- matrixStats::colSds(x.ref, center = cm, na.rm = na.rm)
        } else if (scale.method == "median") {
            csd <- matrixStats::colMads(x.ref, center = cm, na.rm = na.rm)
        }
    } else {
        # just divide by 1 if not
        csd <- rep(1, length = length(cm))
    }

    if (!center) {
        # just subtract 0
        cm <- rep(0, length = length(cm))
    }

    # apply the centering and scaling
    if (scale.method == "mean") {
        x <- t((t(x) - cm) / csd)
    } else if (scale.method == "median") {
        x <- t((0.6745 * (t(x) - cm)) / csd)
    }

    if (add_attr) {
        if (center) {
            attr(x, "scaled:center") <- cm
        }
        if (scale) {
            attr(x, "scaled:scale") <- csd
        }
    }
    return(x)
}


#-------------------------------------------------------------------------------
#' Modified z-score
#'
#' @param x A numeric vectors
#' @return The vector as a modfied z-score
#'
#' @export
mod_zscore <- function(x) {
    return((0.6745 * (x - median(x, na.rm = TRUE))) / mad(x, na.rm = TRUE))
}

#-------------------------------------------------------------------------------
#' Transform assay data using boxcox transform
#'
#' @description Transform a vector using boxcox transform
#' If value is not positive it is offset to make everything positive
#'
#' @param x A numeric vector with values to transform
#' @param return.lambda return the lambda value instead
#' @param limit The absolute value of lambda's to test
#' @param fudge Fudgefactor to apply either log transform, or just leave it as is
#' @param downsample Calculate labmda's on a random sample of the data. Either NULL (no downsampling), an integer (selects x values), or a selection vector
#' @param add.lambda Return a list with two items instead, one transformed x, two, the lamda used
#'
#' @returns A numeric vector with the boxcox transformed data or the optimal lambda value
#' @export
boxcox_transform <- function(x, return.lambda = FALSE, limit = 5, fudge = 0.1, downsample = NULL, add.lambda = FALSE) {
    # If is not positive - then offset to make everything positive
    if (any(x <= 0, na.rm = T) == T) {
        x <- abs(min(x, na.rm = T)) + x + 1
    }

    # Create data to estimate transformation parameters
    # Optionally downsample
    if (!is.null(downsample)) {
        if (class(downsample) == "integer" || class(downsample) == "numeric") {
            if (length(downsample) > 1) {
                y <- x[downsample]
            } else {
                y <- x[sample(seq_along(x), downsample)]
            }
        } else {
            stop("downsample must be a single integer, or a vector representing values to use for lambda estimation")
        }
    } else {
        y <- x
    }

    # Estimate transformation parameters
    bc <- MASS::boxcox(y ~ 1, plotit = F, lambda = seq(-limit, limit, 1 / 10)) # don't plot lambda outcome
    lambda <- bc$x[which.max(bc$y)]

    # If return lambda or transformed values
    if (return.lambda) {
        return(lambda)
    } else {
        # Define fudge: uncertainty value so that anything within that range becomes the closest value

        # Transform data
        if ((lambda < fudge) && (lambda > -fudge)) { # If the data is in between -0.1 & 0.1, then just do a log
            t <- log(x)
            lambda <- 0
        } else if ((lambda < (1 + fudge)) && (lambda > (1 - fudge))) { # If the data is in still between -1.2 & 1.2, then just leav it
            t <- x
            lambda <- 1
        } else { # Otherwise use the calculated lambda
            t <- (x^lambda - 1) / lambda
        }

        if (add.lambda) {
            return(list(x = t, lambda = lambda))
        } else {
            return(t)
        }
    }
}


#-------------------------------------------------------------------------------
#' Transform assay data using boxcox transform
#'
#' @description Transform a assay using boxcox transform. If assay="image.data", result is stored in slot image.data.trans
#' Lambda values are added to the features data.frame: 0=log, 1=no transformation
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param assay.out The name of the output assay. Defaults to paste0(assay, "-trans")
#' @param trim Logical indicating to trim zero variance columns and NA's after applying boxcox transform
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param verbose Raise a warning if columns are trimmed
#' @param rfast.zerotol threshold to consider variances zero in Rfast. See details
#' @param ... Arguments passed to \code{\link{boxcox_transform}}
#'
#' @details
#'
#' `rfast.zerotol`
#' NOTE: Currently NOT using \code{\link{Rfast::colVars()}} as there are inconsistencies with \code{\link{matrixStats::colVars()}}
#' Rfast's does not return zero variances, so we have to set some minial value to consider as 0
#' https://github.com/RfastOfficial/Rfast/issues/62
#' https://github.com/RfastOfficial/Rfast/issues/115
#'
#'
#' @returns \linkS4class{TglowDataset} with new boxcox assay
#' @importFrom progress progress_bar
#' @export
apply_boxcox <- function(dataset, assay, assay.out = NULL, trim = TRUE, slot = "data", verbose = TRUE, rfast.zerotol = 1e-10, ...) {
    # Checks for input
    check_dataset_assay_slot(dataset, assay, slot)

    if (assay == "image.data") {
        mat <- slot(dataset@image.data, slot)@.Data
        features <- dataset@image.data@features
    } else {
        mat <- slot(dataset[[assay]], slot)@.Data
        features <- dataset[[assay]]@features

        if (is.null(assay.out)) {
            assay.out <- paste0(assay, "-trans")
        }
    }
    ncol.orig <- ncol(mat)

    # Apply transform
    pb <- progress_bar$new(format = paste0("[INFO] Transforming [:bar] :current/:total (:percent) eta :eta"), total = ncol(mat))

    lambdas <- c()
    mat <- apply(mat, 2, function(x) {
        pb$tick()
        res <- boxcox_transform(x, add.lambda = T, ...)

        lambdas <<- c(lambdas, res$lambda)
        return(res$x)
    })

    features$lambda <- lambdas
    if (trim) {
        vars <- matrixStats::colVars(mat, na.rm = T)
        # vars[vars < rfast.zerotol] <- 0
        mat <- mat[, vars != 0, drop = F]
        features <- features[vars != 0, ]
    }

    if ((ncol(mat) != ncol.orig) && verbose) {
        warning(paste0("Trimmed ", ncol.orig - ncol(mat), " cols with zero variance or NA's post transform."))
    }

    res <- new("TglowAssay",
        data = TglowMatrix(mat),
        scale.data = TglowMatrix(fast_colscale(mat, na.rm = T)),
        features = features
    )

    if (assay == "image.data") {
        dataset@image.data.trans <- res
    } else {
        dataset@assays[[assay.out]] <- res
    }

    return(dataset)
}


#-------------------------------------------------------------------------------
#' Scale a TglowAssay
#'
#' @description Scale the data slot in a tglow assay and override the scale.data slot
#' Defaults to z-score, but can also return modified z-score when providing scale.method
#'
#' @param assay TglowAssay
#' @param grouping Vector with a grouping variable of length nrow(assay)
#' @param reference.group Scale not to the vector mean/median, sd/mad but to the objects indiciated here
#' @param ... Arguments passed to \code{\link{fast_colscale}}
#' @returns The asssay with the scale.data slot populated
#'
#' @details
#'
#' `grouping`
#'
#' When this is provided, scaling is done within each unique value of this vector
#'
#' `reference.group`
#'
#' When this is provided, the data is scaled to the mean/median, sd/mad of only
#' the samples set to TRUE in this vector. It much be a logical vector indicating
#' if the object is a reference object. This is usefull to scaling imaging features
#' relative to a control sample.
#'
#' `grouping` and `reference.group`
#'
#' All selection vectors are assumed to have the same row order as assay
#' groupung and reference.group interact, so for example if you want to scale to
#' the controls on a plate you would set grouping to getDataByObject(object, "plate")
#' and reference.group to getDataByObject(object, "is_control")
#'
#' @export
scale_assay <- function(assay, grouping = NULL, reference.group = NULL, ...) {
    if (!is(assay, "TglowAssay")) {
        stop("Assay must be TglowAssay")
    }

    if (is.null(grouping)) {
        assay@scale.data <- TglowMatrix(fast_colscale(assay@data@.Data, reference.group = reference.group, ...))
    } else {
        if (length(grouping) != nrow(assay)) {
            stop("Length of grouping must equal nrow(assay)")
        }

        if (!is.null(reference.group)) {
            if (length(grouping) != length(reference.group)) {
                stop("Length of reference.group must equal length of grouping")
            }
            if (!is.logical(reference.group)) {
                stop("Currently only supporting logical reference.group selector to ensure interaction with grouping works")
            }
        }
        mat <- assay@data@.Data
        mat.scale <- matrix(NA, ncol = ncol(mat), nrow = nrow(mat), dimnames = dimnames(mat))

        for (group in unique(grouping)) {
            selector <- grouping == group
            mat.scale[selector, ] <- fast_colscale(mat[selector, ], reference.group = reference.group[selector], ...)
        }

        assay@scale.data <- TglowMatrix(mat.scale)
    }
    return(assay)
}
