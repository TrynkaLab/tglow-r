#-------------------------------------------------------------------------------
#' Faster alternative to scale a matrix
#'
#' @description
#' Sourced from https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' Yields a 10-20% speedup
#'
#' @param x A matrix
#' @param center Should data be centered
#' @param scale Should data be scaled
#' @param add_attr Should means and variances be added to output
#' @param na.rm Should NA's be removed during mean and variance calculations
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
#' @export
fast_colscale <- function(x,
                          center = TRUE,
                          scale = TRUE,
                          add_attr = TRUE,
                          na.rm = F) {
    # Get the column means
    cm <- colMeans(x, na.rm = na.rm)

    # Get the column sd
    if (scale) {
        csd <- matrixStats::colSds(x, center = cm, na.rm = na.rm)
    } else {
        # just divide by 1 if not
        csd <- rep(1, length = length(cm))
    }

    if (!center) {
        # just subtract 0
        cm <- rep(0, length = length(cm))
    }

    # apply the centering and scaling
    x <- t((t(x) - cm) / csd)
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

#' @param x A numeric vector with values to transform
#' @param return.lambda return the lambda value instead
#' @param limit The absolute value of lambda's to test
#' @param fudge Fudgefactor to apply either log transform, or just leave it as is
#' @param downsample Calculate labmda's on a random sample of the data. Either NULL (no downsampling), an integer (selects x values), or a selection vector.
#' @param add.lambda Return a list with two items instead, one transformed x, two, the lamda used

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
#' @description Transform a assay using boxcox transform. If assay="image.data", result is stored in slot image.data.trans.
#' Lambda values are added to the features data.frame: 0=log, 1=no transformation
#'
#' @param dataset A tglow dataset
#' @param assay The assay to use.
#' @param assay.out The name of the output assay. Defaults to paste0(assay, "-trans")
#' @param trim Logical indicating to trim zero variance columns and NA's after applying boxcox transform
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param verbose Raise a warning if columns are trimmed
#' @param rfast.zerotol threshold to consider variances zero in Rfast see details.
#' @param ... Arguments passed to \code{\link{boxcox_transform}}
#'
#' @details
#' NOTE: Currently NOT using Rfast::colVars() as there are inconsistencies with matrixStats::colVars()
#' Rfast's does not return zero variances, so we have to set some minial value to consider as 0
#' https://github.com/RfastOfficial/Rfast/issues/62
#' https://github.com/RfastOfficial/Rfast/issues/115
#'
#'
#'
#' @returns Tglow object with new boxcox assay
#' @export
apply_boxcox <- function(dataset, assay, assay.out = NULL, trim = TRUE, slot = "data", verbose = TRUE, rfast.zerotol = 1e-10, ...) {
    # Checks for input
    if (!is(dataset, "TglowDataset")) {
        stop("dataset must be of class TglowDataset")
    }

    if (assay == "image.data") {
        mat <- slot(dataset@image.data, slot)
        features <- dataset@image.data@features
    } else {
        mat <- slot(dataset[[assay]], slot)
        features <- dataset[[assay]]@features

        if (is.null(assay.out)) {
            assay.out <- paste0(assay, "-trans")
        }
    }
    ncol.orig <- ncol(mat)

    # Apply transform
    pb <- txtProgressBar(min = 0, max = ncol(mat), style = 3)
    lambdas <- c()
    mat <- apply(mat, 2, function(x) {
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
        res <- boxcox_transform(x, add.lambda = T, ...)

        lambdas <<- c(lambdas, res$lambda)
        return(res$x)
    })
    close(pb)

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
        scale.data = TglowMatrix(fast_colscale(mat)),
        features = features
    )

    if (assay == "image.data") {
        dataset@image.data.trans <- res
    } else {
        dataset@assays[[assay.out]] <- res
    }

    return(dataset)
}
