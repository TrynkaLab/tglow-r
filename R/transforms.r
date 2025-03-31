#-------------------------------------------------------------------------------
#' Faster alternative to scale a matrix
#'
#' @description
#' Adapted from https://www.r-bloggers.com/2016/02/a-faster-scale-function/
#' Yields a 10-20 percent speedup and adds modfied z-score scaling
#'
#' @param x A matrix
#' @param center Should data be centered
#' @param scale Should data be scaled
#' @param add_attr Should means and variances be added to output
#' @param na.rm Should NA's be removed during mean and variance calculations
#' @param scale.method Scaling method, either `mean` for z-score or `median` modified z-score
#' @param reference.group Calculate means/medians, sds/mads only on these samples but apply it to all
#' @param na.zero Set values that are exactly 0 to NA when calculating means and variances. See details
#'
#' @returns The scaled matrix
#'
#' @details
#'
#' `scale.method - mean`
#'
#'  Standard centering/scaling: (x-mean(x)) / sd(x)
#'
#' `scale.method - median`
#'
#' Median adjusted deviation centering/scaling: (x-median(x)) / mad(x)
#'
#' `na.zero`
#' This removes values which are exactly zero from mean and variance calculations. Some features
#' might have a lot of zeroes, which can skew the mean and variance estimates and give issues when
#' scaling to reference samples. This sets them to NA when calculating the scaling factors.
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
                          reference.group = NULL,
                          na.zero = F) {
                            
    # This should no longer be needed and was bad
    # if (is(x, "TglowMatrix")) {
    #     x <- x@.Data
    # }

    if (is.null(reference.group)) {
        x.ref <- x
    } else {
        x.ref <- x[reference.group, ]
    }

    if (na.zero) {
        x.ref[x.ref == 0] <- NA
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
#'
#' As trafo's implementation
#' https://github.com/cran/trafo/blob/master/R/trafos.R
#' 
#' @rdname modul
#' @export
modul <- function(y, lambda = lambda) {
  u <- abs(y) + 1L
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <-  sign(y)*log(u)
  } else {
    yt <- sign(y)*(u^lambda - 1L)/lambda
  }
  return(y = yt)
}

#' Standardized transformation: Modulus
#' @rdname modul
#' @export
modul_std <- function(y, lambda) {
  u <- abs(y) + 1L
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    zt <-  sign(y) * log(u) * geometric.mean(u) 
  } else {
    zt <- sign(y)*(u^lambda - 1L)/lambda * (1/geometric.mean(u)^(lambda - 1))
  }
  y <- zt
  
  return(y)
}

#-------------------------------------------------------------------------------
#' Transform vector using modulus transform, finding optimal lambda
#'
#' Modelled after MASS:boxcox() and shares the same basic signature
#' So see [MASS:boxcox()] for help
#' 
#' 
#' @details 
#' Implements modulus transformation described here:
#' https://scales.r-lib.org/reference/boxcox_trans.html
#' John, J. A., & Draper, N. R. (1980). An alternative family of transformations. Applied Statistics, 190-197. http://www.jstor.org/stable/2986305
#' 
#' 
#' @export 
modulus <- function(object, lambda = seq(-2, 2, 1/10), plotit =  TRUE,
         interp = (plotit && (m < 100)), eps = 1/50,
         xlab = expression(lambda), ylab = "log-Likelihood", ...) {
            
    m <- length(lambda)
    object <- lm(object, y = TRUE, qr = TRUE, ...)
    
    if(is.null(object$y) || is.null(object$qr))
        object <- update(object, y = TRUE, qr = TRUE, ...)

    if(is.null(y <- object$y) || is.null(xqr <- object$qr))
        stop(gettextf("%s does not have both 'qr' and 'y' components",
                      sQuote(deparse(substitute(object)))), domain = NA)
   

    n <- length(y)
    ## scale y[]  {for accuracy in  y^la - 1 }:
    #y    <- y / exp(mean(log(y)))
    #logy <- log(y) # now  ydot = exp(mean(log(y))) == 1
    
    xl   <- loglik <- as.vector(lambda)
    m    <- length(xl)
    
    for(i in 1L:m) {
        if(abs(la <- xl[i]) > eps) {
            # This is the boxcox, replace it with modulus transform
            #yt <- (y^la - 1)/la
            # https://scales.r-lib.org/reference/boxcox_trans.html
            
            # Raw modulus transform
            #yt <- sign(y) * ((abs(y) + 1)^la - 1) / la
            # I don't 100% get why, but weighting by the genometric mean is important
            # otherwise the likelihood just increases with higher lambda's and you cannot
            # find the actual best value. I assume this ensures the actual value is close to the mean
            yt <- modul_std(y, la)
        } else {
            #yt <- logy * (1 + (la * logy)/2 *
            #                (1 + (la * logy)/3 * (1 + (la * logy)/4)))
            #yt <- sign(y) * log(abs(y) + 1)
            yt <- modul_std(y, 0)
        }
        loglik[i] <- - n/2 * log(sum(qr.resid(xqr, yt)^2))
    }
    if(interp) {
        sp <- spline(xl, loglik, n = 100)
        xl <- sp$x
        loglik <- sp$y
        m <- length(xl)
    }
    if(plotit) {
        mx <- (1L:m)[loglik == max(loglik)][1L]
        Lmax <- loglik[mx]
        lim <- Lmax - qchisq(19/20, 1)/2
        dev.hold(); on.exit(dev.flush())
        plot(xl, loglik, xlab = xlab, ylab = ylab, type
             = "l", ylim = range(loglik, lim))
        plims <- par("usr")
        abline(h = lim, lty = 3)
        y0 <- plims[3L]
        scal <- (1/10 * (plims[4L] - y0))/par("pin")[2L]
        scx <- (1/10 * (plims[2L] - plims[1L]))/par("pin")[1L]
        text(xl[1L] + scx, lim + scal, " 95%", xpd = TRUE)
        la <- xl[mx]
        if(mx > 1 && mx < m)
            segments(la, y0, la, Lmax, lty = 3)
        ind <- range((1L:m)[loglik > lim])
        if(loglik[1L] < lim) {
            i <- ind[1L]
            x <- xl[i - 1] + ((lim - loglik[i - 1]) *
                              (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
            segments(x, y0, x, lim, lty = 3)
        }
        if(loglik[m] < lim) {
            i <- ind[2L] + 1
            x <- xl[i - 1] + ((lim - loglik[i - 1]) *
                              (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
            segments(x, y0, x, lim, lty = 3)
        }
    }
    list(x = xl, y = loglik)
}


#-------------------------------------------------------------------------------
#' Transform assay data using boxcox transform
#'
#' @description Transform a vector using boxcox transform
#' If value is not positive it is offset to make everything positive
#'
#' @param x A numeric vector with values to transform
#' @param return.lambda return the lambda value instead
#' @param mode Either 'boxcox' or 'modulus', see details
#' @param limit The absolute value of lambda's to test
#' @param fudge Fudgefactor to apply either log transform, or just leave it as is
#' @param downsample Calculate labmda's on a random sample of the data. Either NULL (no downsampling), an integer (selects x values), or a selection vector
#' @param add.lambda Return a list with two items instead, one transformed x, two, the lamda used
#' @param filter.iqr Should values outside the 2xIQR be removed prior to estimating lambda
#'
#' @returns A numeric vector with the boxcox transformed data or the optimal lambda value
#' 
#' @details 
#' In mode boxcox [MASS::boxcox()] is fit on the offset values if there are any values < 0. Values are made positive by adding the min absolute value
#' In mode modulus the modulus transform is applied instead, which is a generalization of boxcox. See [modulus_transform()]
#' 
#' @export
boxcox_transform <- function(x, return.lambda = FALSE, mode="boxcox", limit = 5, fudge = 0.2, downsample = NULL, add.lambda = FALSE, filter.iqr=FALSE) {
    
    
    if (mode == "boxcox") {
        # If is not positive - then offset to make everything positive
        if (any(x <= 0, na.rm = T) == T) {
            x <- abs(min(x, na.rm = T)) + x + 1
        }
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

    # Filter 2x IQR outliers
    if (filter.iqr) {
        q1 <- quantile(y, 0.25, na.rm=T) 
        q3 <- quantile(y, 0.75, na.rm=T)
        iqr   <- q3 - q1
        lower <- q1 - 2 * iqr
        upper <- q3 + 2 * iqr
        y <- y[y >= lower & y <= upper]   
        #x[x >= lower & x <= upper] <- NA
    }

    # Estimate transformation parameters
    if (mode == "boxcox") {
        bc <- MASS::boxcox(y ~ 1, plotit = F, lambda = seq(-limit, limit, 1 / 10), interp=F) # don't plot lambda outcome
    } else if (mode == "modulus") {
        # Use modulus transform instrad
        bc <- modulus(y ~ 1, plotit = F, lambda = seq(-limit, limit, 1 / 10), interp=F) # don't plot lambda outcome
    } else {
        stop(paste0(mode, " is not a valid mode, bust be boxcox or"))
    }

    lambda <- bc$x[which.max(bc$y)]

    # If return lambda or transformed values
    if (return.lambda) {
        return(lambda)
    } else {
        # Define fudge: uncertainty value so that anything within that range becomes the closest value

        # Transform data
        if ((lambda <= fudge) && (lambda >= -fudge)) { # If the data is in between -0.1 & 0.1, then just do a log
            
            if (mode == "boxcox") {t <- log(x)}
            if (mode == "modulus") {t <- modul(x, 0)}
            lambda <- 0
        } else if ((lambda <= (1 + fudge)) && (lambda >= (1 - fudge))) { # If the data is in still between -1.2 & 1.2, then just leav it
            t <- x
            lambda <- 1
        } else { # Otherwise use the calculated lambda
            if (mode == "boxcox") {t <- (x^lambda - 1) / lambda}
            if (mode == "modulus") {t <- modul(x, lambda)}
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
#' @param filter.iqr Remove values outside +- 2x IQR to reduce impact of outliers
#' @param ... Arguments passed to \code{\link{boxcox_transform}}
#'
#' @details
#' Calls boxcox_transform on an assay, please read  \code{\link{boxcox_transform}} for more relvant details and options
#'
#'
#' @returns \linkS4class{TglowDataset} with new boxcox assay
#' @importFrom progress progress_bar
#' @export
apply_boxcox <- function(dataset, assay, assay.out = NULL, trim = TRUE, slot = "data", verbose = TRUE, filter.iqr=F, ...) {
    # Checks for input
    check_dataset_assay_slot(dataset, assay, slot)

    if (assay == "image.data") {
        #mat <- slot(dataset@image.data, slot)@.Data
        mat <- slot(dataset@image.data, slot)
        features <- dataset@image.data@features
    } else {
        #mat <- slot(dataset[[assay]], slot)@.Data
        mat <- slot(dataset[[assay]], slot)
        features <- dataset[[assay]]@features

        if (is.null(assay.out)) {
            assay.out <- paste0(assay, "-trans")
        }
    }
    ncol.orig <- ncol(mat)

    # Apply transform
    pb <- progress::progress_bar$new(format = paste0("[INFO] Transforming [:bar] :current/:total (:percent) eta :eta"), total = ncol(mat))

    lambdas <- c()
    mat <- apply(mat, 2, function(x, filter.iqr) {
        pb$tick()
        
        res <- tglowr::boxcox_transform(x, add.lambda = T, filter.iqr=filter.iqr, ...)
        lambdas <<- c(lambdas, res$lambda)
        return(res$x)
    }, filter.iqr=filter.iqr)

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
        data = tglowr::TglowMatrix(mat),
        scale.data = tglowr::TglowMatrix(tglowr::fast_colscale(mat, na.rm = T)),
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
#' @param assay \linkS4class{TglowAssay}
#' @param grouping Vector with a grouping variable of length nrow(assay)
#' @param reference.group Scale not to the vector mean/median, sd/mad but to the objects indiciated here
#' @param ... Arguments passed to \code{\link{fast_colscale}}. Strongly reccomend looking at these!
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
        #assay@scale.data <- tglowr::TglowMatrix(tglowr::fast_colscale(assay@data@.Data, reference.group = reference.group, ...))
        assay@scale.data <- tglowr::TglowMatrix(tglowr::fast_colscale(assay@data, reference.group = reference.group, ...))
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
        #mat <- assay@data@.Data
        mat <- assay@data
        mat.scale <- matrix(NA, ncol = ncol(mat), nrow = nrow(mat), dimnames = dimnames(mat))

        for (group in unique(grouping)) {
            selector <- grouping == group
            mat.scale[selector, ] <- tglowr::fast_colscale(mat[selector, ], reference.group = reference.group[selector], ...)
        }

        assay@scale.data <- tglowr::TglowMatrix(mat.scale)
    }
    return(assay)
}

#-------------------------------------------------------------------------------
#' Scale a tglow dataset
#'
#' @description  For TglowDataset, scales the provided assay, or all assays when assay is null
#' @param dataset \linkS4class{TglowDataset}
#' @param assay Character string indicating assay, if NULL all assays are scaled.
#' @param grouping Colname on dataset or Vector with a grouping variable of length nrow(assay)
#' @param ... Remaining options passed to [scale_assay()]
#'
#' @returns The \linkS4class{TglowDataset} with assay scale.data slots populated
#' @export
scale_dataset <- function(dataset, assay = NULL, grouping = NULL, ...) {
    if (is.character(grouping) && length(grouping) == 1) {
        if (tglowr::isAvailable(dataset, grouping, assay = assay, slot = "data")) {
            grouping <- tglowr::getDataByObject(dataset, assay = assay)[, 1]
        } else {
            stop(paste0(grouping, " is not available meta, image.meta, or assay item on dataset"))
        }
    }

    if (!is.null(assay)) {
        check_dataset_assay_slot(dataset, assay = assay, slot = "data")
        dataset@assays[[assay]] <- tglowr::scale_assay(dataset@assays[[assay]], grouping = grouping, ...)
    } else {
        for (assay in names(dataset@assays)) {
            cat("[INFO] Scaling assay: ", assay, "\n")
            dataset@assays[[assay]] <- tglowr::scale_assay(dataset@assays[[assay]], grouping = grouping, ...)
        }
    }

    return(dataset)
}



#-------------------------------------------------------------------------------
#' Scale a TglowAssay between a min and max value for each feature
#'
#' @description Scale the data slot between min and max and return a new assay
#'
#' @param assay \linkS4class{TglowAssay}
#' @param slot The slot on assay to draw from. Should be 'data' in most cases.
#' @param min The minimal value for each feature
#' @param max The maximum value for each feature
#' @param scale.by.var Multiply the result by the orignal variance
#' @param scaling.factors Vector of scaling factors per feature
#'
#' @returns A \linkS4class{TglowAssay} with the data slot populated with the rescaled data.
#' @export
scale_assay_min_max <- function(assay, slot = "data", min = 0, max = 1, scale.by.var = FALSE, scaling.factors = NULL) {
    if (!is(assay, "TglowAssay")) {
        stop("Assay must be TglowAssay")
    }

    #m <- slot(assay, slot)@.Data
    m <- slot(assay, slot)

    # return (maxAllowed - minAllowed) * (unscaledNum - min) / (max - min) + minAllowed;

    m <- apply(m, 2, function(x) {
        xx <- (((max - min) * (x - min(x))) / (max(x) - min(x))) + min
        return(xx)
    })

    features <- assay@features

    # Scale each column by scaling factors or the orignal variance
    if (!is.null(scaling.factors)) {
        if (length(scaling.factors) != ncol(m)) {
            stop("scaling.factors must equal ncol(assay) in length")
        }
        m <- m %*% diag(scaling.factors)
        features[,"scaling_factors"] <- scaling.factors
    } else if (scale.by.var) {
        #scaling.factors <- matrixStats::colVars(slot(assay, slot)@.Data)
        scaling.factors <- matrixStats::colVars(slot(assay, slot))
        m <- m %*% diag(scaling.factors)
        features[,"scaling_factors"] <- scaling.factors
    }

    return(new("TglowAssay",
        data = tglowr::TglowMatrix(m),
        features = features
    ))
}
