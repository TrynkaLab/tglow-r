#-------------------------------------------------------------------------------
#' Find marker features using a t-test
#'
#' @description
#' Find marker features using a two sample t-test. By default top 10 positive results
#' for a class are returned
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param ident A column in meta, image.meta, assay, assay.image to use as class labels
#' @param assay The assay to use
#' @param assay.image The image assay to use for grabbing ident, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param slot The slot to use for finding markers. Can be "data" or "scale.data"
#' @param return.top Return the top x values by t-statistic. Defaults to 10
#' @param ref.classes A vector of values to compare class mean against. Default NULL uses all other classes as reference
#' @param na.rm Should NA's be removed during t-test
#'
#' @returns A data.frame with t-test results
#' @importFrom progress progress_bar
#' @export
find_markers <- function(dataset, ident, assay, slot, assay.image = NULL, return.top = 10, ref.classes = NULL, na.rm = T) {
    # Check input
    check_dataset_assay_slot(dataset, assay, slot)

    cur.assay <- slot(dataset[[assay]], slot)@.Data
    cur.ident <- as.character(getDataByObject(dataset, ident, assay = assay, assay.image = assay.image, slot = slot))

    if (is.null(cur.ident)) {
        stop("Ident not valid")
    }

    classes <- unique(cur.ident)

    if (length(classes) <= 1) {
        stop("Ident only has one class")
    }

    res <- data.frame(matrix(NA, nrow = ncol(cur.assay) * length(classes), ncol = 9))
    colnames(res) <- c("class", "feature", "t-stat", "df", "pval", "mean.diff", "mean.se", "mean.ref", "mean.class")
    i <- 1

    pb <- progress_bar$new(format = "[INFO] Finding markers [:bar] :current/:total (:percent) eta :eta", total = ncol(cur.assay) * length(classes))
    for (class in classes) {
        ident.is.class <- cur.ident == class

        if (is.null(ref.classes)) {
            ident.is.ref <- !ident.is.class
        } else {
            ident.is.ref <- (cur.ident != class) & (cur.ident %in% ref.classes)
        }

        for (col in colnames(cur.assay)) {
            pb$tick()

            tmp <- t.test(cur.assay[ident.is.class, col], cur.assay[ident.is.ref, col], na.rm = na.rm)

            res[i, "class"] <- class
            res[i, "feature"] <- col
            res[i, "t-stat"] <- tmp$statistic
            res[i, "df"] <- tmp$parameter
            res[i, "pval"] <- tmp$p.value
            res[i, "mean.diff"] <- tmp$estimate[1] - tmp$estimate[2]
            res[i, "mean.se"] <- tmp$stderr
            res[i, "mean.ref"] <- tmp$estimate[2]
            res[i, "mean.class"] <- tmp$estimate[1]

            i <- i + 1
        }
    }


    if (!is.null(return.top)) {
        subsets <- list()
        for (class in classes) {
            tmp <- res[res$class == class, ]
            tmp <- tmp[order(tmp$`t-stat`, decreasing = T), ]
            subsets[[class]] <- tmp[1:return.top, ]
        }

        res <- do.call(rbind, subsets)
    }

    return(res)
}


#-------------------------------------------------------------------------------
#' Linearly correct for a set of covariates
#'
#' @description Fit a linear model using OLS and correct an assay for specified covariates
#' @param object A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for regressing against. Can be "data" or "scale.data"
#' @param covariates Character vector of covariates to correct for
#' @param slot.covar The slot to grab covariates from. Can be "data" or "scale.data". Defaults to slot
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param formula The formula to use for regression. Defaults to additive model. See details
#' @param assay.out Name of the output assay. Defaults to <assay>.lm.corrected
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param covariates.dont.use Vector of covariate names to NOT use when calculating residuals. See detaills
#'
#' @details
#'
#' `grouping`
#'
#'  If this is provided, scaling for populating scale.data slot is done over ALL residuals, not per group to ensure the mean and sd of the whole vector is as expected
#'
#' `formula`
#'
#'  If NULL an additive model of all covariates is performed. Otherwise should be a string interpretable by \code{\link{base::as.formula}}
#'
#' `covariates.dont.correct`
#'
#' The beta's for these variables are removed when calculating the residuals. When specifying more complex models in formula,
#' use the term names, with for instance, interaction terms for example, if it has a form of '~ a + b + c + b:c' and you don't
#' want to consider the interaction term, add 'b:c'. To remove the intercept, add '(Intercept)'
#'
#' @returns The \linkS4class{TglowDataset} with a corrected assay
#' @importFrom progress progress_bar
#' @export
correct_lm <- function(object, assay, slot, covariates, slot.covar = NULL, assay.image = NULL, formula = NULL, assay.out = NULL, grouping = NULL, covariates.dont.use = NULL) {
    check_dataset_assay_slot(object, assay, slot)

    if (is.null(slot.covar)) {
        slot.covar <- slot
    }

    data <- getDataByObject(object, covariates, assay, assay.image, slot.covar, drop = F)

    if (is.null(formula)) {
        design <- model.matrix(~., data = data)
    } else {
        design <- model.matrix(formula, data = data)
    }

    response <- slot(object@assays[[assay]], slot)@.Data
    residuals <- matrix(NA, nrow = nrow(response), ncol = ncol(response), dimnames = dimnames(response))

    if (is.null(grouping)) {
        grouping <- rep(1, nrow(response))
    }

    for (group in unique(grouping)) {
        selector <- grouping == group
        residuals[selector, ] <- lm_matrix(response[selector, ], design[selector, ], covariates.dont.use = covariates.dont.use, residuals.only = TRUE)
    }

    if (is.null(assay.out)) {
        assay.out <- paste0(assay, ".lm.corrected")
    }

    cat("[INFO] Regressions done. Scaling residuals\n")

    object@assays[[assay.out]] <- new("TglowAssay",
        data = TglowMatrix(residuals),
        scale.data = TglowMatrix(fast_colscale(residuals)),
        features = object@assays[[assay]]@features
    )

    return(object)
}

#-------------------------------------------------------------------------------
#' Linearly correct for a set of covariates
#'
#' @description Fit a linear model using OLS and correct an assay for specified covariates
#' @param object A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for regressing against. Can be "data" or "scale.data"
#' @param covariates.group List with specific covariates for groups of features. See detaills
#' @param slot.covar The slot to grab covariates from. Can be "data" or "scale.data". Defaults to slot
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param assay.out Name of the output assay. Defaults to <assay>.lm.corrected
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param covariates.dont.use Vector of covariate names to NOT use when calculating residuals. See detaills
#'
#' @details
#'
#' `covariates.group`
#'
#' This must be a list where the name of the list item is a grep pattern of the features to apply the covariates to, and the item is a vector of feature names to correct for
#'
#' `grouping`
#'
#'  If this is provided, scaling for populating scale.data slot is done over ALL residuals, not per group to ensure the mean and sd of the whole vector is as expected
#'
#' `covariates.dont.correct`
#'
#' The beta's for these variables are removed when calculating the residuals. When specifying more complex models in formula,
#' use the term names, with for instance, interaction terms for example, if it has a form of '~ a + b + c + b:c' and you don't
#' want to consider the interaction term, add 'b:c'. To remove the intercept, add '(Intercept)'
#'
#' @returns The \linkS4class{TglowDataset} with a corrected assay
#' @importFrom progress progress_bar
#' @export
correct_lm_per_featuregroup <- function(object, assay, slot, covariates.group, slot.covar = NULL, assay.image = NULL, assay.out = NULL, grouping = NULL, covariates.dont.use = NULL) {
    check_dataset_assay_slot(object, assay, slot)

    if (!is.list(covariates.group)) {
        stop("grouping.features must be a list")
    }

    if (is.null(names(covariates.group))) {
        stop("Names attribute of grouping features must be a grep compatible pattern")
    }

    for (fgroup in names(covariates.group)) {
        if (!is.character(covariates.group[[fgroup]])) {
            stop(paste0(fgroup, " must be a character vector indicating the covariates to correct group for"))
        }
    }

    feature.colnames <- list()
    for (fgroup in names(covariates.group)) {
        feature.colnames[[fgroup]] <- grep(fgroup, colnames(object@assays[[assay]]), value = TRUE)
    }

    if (list_has_overlap(feature.colnames)) {
        stop("Overlap found between features selected by covariates.group, make sure your grep patterns resolve to unique sets!")
    }

    if (is.null(slot.covar)) {
        slot.covar <- slot
    }

    response <- slot(object@assays[[assay]], slot)@.Data
    residuals <- matrix(NA, nrow = nrow(response), ncol = ncol(response))
    rownames(residuals) <- rownames(response)
    colnames(residuals) <- colnames(response)

    if (is.null(grouping)) {
        grouping <- rep(1, nrow(response))
    }

    for (group in unique(grouping)) {
        selector <- grouping == group

        for (fgroup in names(covariates.group)) {
            data <- getDataByObject(object, covariates.group[[fgroup]], assay, assay.image, slot.covar, drop = F)
            design <- model.matrix(~., data = data[selector, , drop = FALSE])
            group.features <- feature.colnames[[fgroup]]
            cat("[INFO] Running object group: ", group, " for ", fgroup, " and selected ", length(group.features), " features\n")

            res <- lm_matrix(response[selector, group.features],
                design,
                covariates.dont.use = covariates.dont.use,
                residuals.only = TRUE
            )

            residuals[selector, group.features] <- res
        }
    }

    if (is.null(assay.out)) {
        assay.out <- paste0(assay, ".lm.corrected.featuregroup")
    }

    cat("[INFO] Regressions done. Scaling residuals\n")

    object@assays[[assay.out]] <- new("TglowAssay",
        data = TglowMatrix(residuals),
        scale.data = TglowMatrix(fast_colscale(residuals)),
        features = object@assays[[assay]]@features
    )

    return(object)
}


#-------------------------------------------------------------------------------
#' Calculate linear coefficients
#'
#' @description Fit a linear model using OLS and find coefficients
#' @param object A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param slot The slot to use for regressing against. Can be "data" or "scale.data"
#' @param slot.covar The slot to grab covariates from. Can be "data" or "scale.data"
#' @param covariates Character vector of in the model
#' @param formula The formula to use for regression. Defaults to additive model. See details
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param covariates.dont.use Only use if you understand the implications. See detaills
#'
#' @details
#' `grouping`
#'
#'  Runs seperate regressions per group in the indicated grouping variable. length(grouping) must be nrow(assay)
#'
#' `formula`
#'
#'  If NULL an additive model of all covariates is performed. Otherwise should be a string interpretable by \code{\link{base::as.formula}}
#'
#' `covariates.dont.use`
#'
#' Fit the model with these covariates included, but don't use them for calculating residuals (affects r2, se etc)
#' The covariates for these variables are removed in the output. When specifying more complex models in formula
#' use the term names, with for instance, interaction terms for example, if it has a form of '~ a + b + c + b:c' and you don't
#' want to consider the interaction term, add 'b:c'. To remove the intercept, add '(Intercept)'
#'
#' @returns A list of regression results. If grouping != NULL, there is one list per group
#' @export
calculate_lm <- function(object, assay, assay.image = NULL, slot, slot.covar = "scale.data", covariates, formula = NULL, grouping = NULL, covariates.dont.use = NULL) {
    check_dataset_assay_slot(object, assay, slot)

    data <- getDataByObject(object, covariates, assay, assay.image, slot.covar, drop = F)

    if (ncol(data) <= 0) {
        stop("data (covariates) cannot be empty. Are your collumn names correct?")
    }

    if (is.null(formula)) {
        design <- model.matrix(~., data = data)
    } else {
        design <- model.matrix(formula, data = data)
    }

    response <- slot(object@assays[[assay]], slot)@.Data

    if (nrow(design) != nrow(response)) {
        stop("nrow(design) must equal nrow(assay)")
    }

    if (is.null(grouping)) {
        res <- lm_matrix(response, design, covariates.dont.use = covariates.dont.use)
        return(res)
    } else {
        if (length(grouping) != nrow(response)) {
            stop("grouping must have the same length as nrow(assay)")
        }
        results <- list()
        for (group in unique(grouping)) {
            cat("[INFO] Starting regressions for group: ", group, "\n")
            selector <- grouping == group
            results[[group]] <- lm_matrix(response[selector, ], design[selector, ], covariates.dont.use = covariates.dont.use)
        }
        return(results)
    }
}

#-------------------------------------------------------------------------------
#' Calculate linear coefficients
#'
#' @description Fit a linear model using OLS and find coefficients
#' @param response A matrix with response variables
#' @param design The common design matrix to regress response againsts
#' @param covariates.dont.use Only use if you understand the implications. See detaills
#' @param residuals.only Only return the residual matrix
#' @param return.residuals Return residual matrix in the output list. Defaults to T if residual.only = TRUE
#'
#' @details
#' This is more efficient as the b component of the design matrix can be re-used between regressions
#' b <- chol2inv(chol(crossprod(design)))
#'
#' `covariates.dont.use`
#' Fit the model with these covariates included, but don't use them for calculating residuals (affects r2, se etc)
#' The covariates for these variables are removed in the output. When specifying more complex models in formula
#' use the term names, with for instance, interaction terms for example, if it has a form of '~ a + b + c + b:c' and you don't
#' want to consider the interaction term, add 'b:c'. To remove the intercept, add '(Intercept)'
#'
#' @returns A list with regression results
#' @importFrom progress progress_bar
#' @export
lm_matrix <- function(response, design, covariates.dont.use = NULL, residuals.only = FALSE, return.residuals = FALSE) {
    if (residuals.only && !return.residuals) {
        return.residuals <- TRUE
    }

    # Matrices to save model coefficients
    if (!residuals.only) {
        coef <- matrix(NA, nrow = ncol(response), ncol = sum(!colnames(design) %in% covariates.dont.use))
        se <- matrix(NA, nrow = ncol(response), ncol = sum(!colnames(design) %in% covariates.dont.use))
        model.stats <- matrix(NA, nrow = ncol(response), 5)
        rownames(coef) <- colnames(response)
        colnames(coef) <- colnames(design)[!colnames(design) %in% covariates.dont.use]
        rownames(se) <- colnames(response)
        colnames(se) <- colnames(design)[!colnames(design) %in% covariates.dont.use]
        rownames(model.stats) <- colnames(response)
        colnames(model.stats) <- c("r2", "adj.r2", "f-stat", "p-value", "df")
    }

    # Matrix to save residuals
    if (return.residuals) {
        residuals <- matrix(NA, nrow = nrow(response), ncol = ncol(response))
        rownames(residuals) <- rownames(response)
        colnames(residuals) <- colnames(response)
    }

    # Calculate the beta's. Use chol2inv on the cholesky decomposition
    # to invert rather then solve as this is faster. This only works on
    # positive-definite matrices, but I think this should always be true in this case
    # (dont quote me on this). If not it throws an error so you can use solve instead

    # Pre-calculate b component on design matrix
    b <- chol2inv(chol(crossprod(design)))
    b.tmp <- b[!colnames(design) %in% covariates.dont.use, !colnames(design) %in% covariates.dont.use]
    cat("[INFO] Starting regressions\n")

    pb <- progress_bar$new(format = paste0("[INFO] Regressing [:bar] :current/:total (:percent) eta :eta"), total = ncol(response))

    for (col in seq_len(ncol(response))) {
        pb$tick()
        # Calculate effectiszes
        beta <- crossprod(b, crossprod(design, response[, col]))

        # Residuals
        if (is.null(covariates.dont.use)) {
            beta.tmp <- beta
            design.tmp <- design
        } else {
            # Include covariates for beta fitting, but don't correct for them
            beta.tmp <- beta[!colnames(design) %in% covariates.dont.use, ]
            design.tmp <- design[, !colnames(design) %in% covariates.dont.use]
        }
        rs <- response[, col] - (design.tmp %*% beta.tmp)

        # Optionally save residuals
        if (return.residuals) {
            residuals[, col] <- rs
        }

        # Optionally save model statistics
        if (!residuals.only) {
            # Residual and total sum of squares
            rss <- crossprod(rs)
            tss <- crossprod(response[, col] - mean(response[, col]))

            # RSS / df
            df <- (length(rs) - ncol(design))
            mse <- as.numeric(rss / df)
            se[col, ] <- as.numeric(sqrt(diag(mse * b)))[!colnames(design) %in% covariates.dont.use]
            coef[col, ] <- as.numeric(beta.tmp)

            r2 <- 1 - (rss / tss)
            r2.adj <- 1 - (1 - r2) * (length(rs) - 1) / df

            # Calculate F-statistic
            msr <- (tss - rss) / (ncol(design.tmp) - 1)
            f.stat <- msr / mse
            p <- 1 - pf(f.stat, ncol(design.tmp) - 1, df)

            model.stats[col, ] <- c(r2, r2.adj, f.stat, p, df)
        }
    }

    # Return only residual matrix
    if (residuals.only) {
        return(residuals)
    }

    # Return output list
    out.list <- list(coef = coef, se = se, model.stats = model.stats, df = df, df.m = ncol(design.tmp) - 1, residuals = NULL)
    if (return.residuals) {
        out.list$residuals <- residuals
    }
    return(out.list)
}
