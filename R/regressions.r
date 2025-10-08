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

    #cur.assay <- slot(dataset[[assay]], slot)@.Data
    cur.assay <- slot(dataset[[assay]], slot)
    cur.ident <- as.character(tglowr::getDataByObject(dataset, ident, assay = assay, assay.image = assay.image, slot = slot))

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

    pb <- progress::progress_bar$new(format = "[INFO] Finding markers [:bar] :current/:total (:percent) eta :eta", total = ncol(cur.assay) * length(classes))
    for (class in classes) {
        ident.is.class <- cur.ident == class

        if (is.null(ref.classes)) {
            ident.is.ref <- !ident.is.class
        } else {
            ident.is.ref <- (cur.ident != class) & (cur.ident %in% ref.classes)
        }

        skip <- FALSE
        if (class %in% ref.classes) {
            msg <- paste0("Skipping class ", class, " as it is in ref.classes, will return NA for this class")
            cat("[WARN] ", msg, "\n")
            warning(msg)
            skip <- TRUE
        } else if (sum(ident.is.ref) == 0) {
            warning(paste0("No reference items found. Consider setting ref.classes. skipping class ", class))
            skip <- TRUE
        } else if (sum(ident.is.class) == 0) {
            warning(paste0("No class items found, skipping class ", class))
            skip <- TRUE
        }

        if (skip) {
            res[i:(i + (ncol(cur.assay) - 1)), "class"] <- class
            res[i:(i + (ncol(cur.assay) - 1)), "feature"] <- colnames(cur.assay)
            i <- i + ncol(cur.assay)
            pb$tick(ncol(cur.assay))
            next()
        }

        for (col in colnames(cur.assay)) {
            pb$tick()

            x <- cur.assay[ident.is.class, col]
            y <- cur.assay[ident.is.ref, col]

            if (is.null(x) || is.null(y)) {
                msg <- paste0("Null after subsetting, skipping class:feature (", class, ":", col, ")")
                warning(msg)
                #cat("\n[WARN] ", msg, "\n")
                res[i, "class"] <- class
                res[i, "feature"] <- col
                i <- i + 1
                next()
            }

            # Double check for NA values
            if ((sum(is.na(x)) >= (length(x) - 3)) || (sum(is.na(y)) >= (length(y) - 3)) || var(x, na.rm = TRUE) < 1e-16 || var(y, na.rm = TRUE) < 1e-16) {
                msg <- paste0("Need at least 3 non NA's in class or reference groups and need variance > 1e-16, skipping class:feature (", class, ":", col, ")")
                warning(msg)
                #cat("\n[WARN] ", msg, "\n")
                res[i, "class"] <- class
                res[i, "feature"] <- col
                i <- i + 1
                next()
            }

            tmp <- t.test(x, y, na.rm = na.rm)

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
#' @param assay.covar The assay to grab covariates from. Defaults to assay argument
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param formula The formula to use for regression. Defaults to additive model. See details
#' @param assay.out Name of the output assay. Defaults to <assay>.lm.corrected
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param covariates.dont.use Vector of covariate names to NOT use when calculating residuals. See detaills
#' @param rescale.group When grouping is active, should the subset be re-centered and scaled prior to regressing
#' @param ... Remaining arguments passed to [lm_matrix()]. Use with caution.
#' @details
#'
#' `grouping`
#'
#' This fits a model and calculates residuals per group seperately. 
#' This can introduce a bias which can make groups uncomparable depending on the covariate structure. 
#' Generally we do not reccomend using this unless you understand the implications.
#' If this is provided, scaling for populating scale.data slot is done over ALL residuals, 
#' not per group to ensure the mean and sd of the whole vector is as expected
#'
#' `formula`
#'
#' If NULL an additive model of all covariates is performed. Otherwise should be a string interpretable by [stats::formula()]
#'
#' `covariates.dont.use`
#'
#' The beta's for these variables are removed when calculating the residuals. When specifying more complex models in formula,
#' use the term names, with for instance, interaction terms for example, if it has a form of '~ a + b + c + b:c' and you don't
#' want to consider the interaction term, add 'b:c'. To remove the intercept, add '(Intercept)'
#'
#' @returns The \linkS4class{TglowDataset} with a corrected assay
#' @importFrom progress progress_bar
#' @export
correct_lm <- function(object, assay, slot, covariates, slot.covar = NULL, assay.covar = NULL, assay.image = NULL, formula = NULL, assay.out = NULL, grouping = NULL, covariates.dont.use = NULL, rescale.group = FALSE, ...) {
    check_dataset_assay_slot(object, assay, slot)

    if (is.null(slot.covar)) {
        slot.covar <- slot
    }
    
    if (is.null(assay.covar)) {
        assay.covar <- assay
    }
    
    if (!is.null(grouping)) {
        if (length(grouping) != nrow(object)) {
        stop("grouping must be a vector of nrow(object)")
        }
    }


    data <- getDataByObject(object, covariates, assay.covar, assay.image, slot.covar, drop = F)

    covariates.dont.use <- check_unused_covar(data, covariates.dont.use)

    if (is.null(formula)) {
        design <- model.matrix(~., data = data)
    } else {
        
        if (is.character(formula)) {
            formula <- as.formula(formula)
            warning(paste0("Formula is character, converting to formula object: ", paste0(as.character(formula), collapse=" ")))
        }
        design <- model.matrix(formula, data = data)
    }

    #response <- slot(object@assays[[assay]], slot)@.Data
    response  <- slot(object@assays[[assay]], slot)
    residuals <- matrix(NA, nrow = nrow(response), ncol = ncol(response), dimnames = dimnames(response))

    if (is.null(grouping)) {
        grouping <- rep(1, nrow(response))
    }

    for (group in unique(grouping)) {
        selector <- grouping == group

        if (rescale.group) {
            cat("[INFO] Rescaling within ", group, ". Set rescale.group=F if this is not what you want\n")
            response.cur <- fast_colscale(response[selector, ])
        } else {
            response.cur <- response[selector, ]
        }

        residuals[selector, ] <- lm_matrix(response.cur, design[selector, ], covariates.dont.use = covariates.dont.use, residuals.only = TRUE, ...)
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
#' @param assay.covar The assay to grab covariates from. Defaults to assay argument
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param assay.out Name of the output assay. Defaults to <assay>.lm.corrected
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param covariates.dont.use Vector of covariate names to NOT use when calculating residuals. See detaills
#' @param rescale.group When grouping is active, should the group be re-centered and scaled prior to regressing
#' @param na.rm Remove NA's from the design matrix prior to regressing (does not remove NAs in the response)
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
correct_lm_per_featuregroup <- function(object, assay, slot, covariates.group, slot.covar = NULL, assay.covar = NULL, assay.image = NULL, assay.out = NULL, grouping = NULL, covariates.dont.use = NULL, rescale.group = FALSE, na.rm=F) {
    check_dataset_assay_slot(object, assay, slot)

    if (!is.list(covariates.group)) {
        stop("grouping.features must be a list")
    }

    if (is.null(names(covariates.group))) {
        stop("Names attribute of grouping features must be a grep compatible pattern")
    }

    
    if (!is.null(grouping)) {
        if (length(grouping) != nrow(object)) {
        stop("grouping must be a vector of nrow(object)")
        }
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

    if (is.null(assay.covar)) {
        assay.covar <- assay
    }


    #response <- slot(object@assays[[assay]], slot)@.Data
    response <- slot(object@assays[[assay]], slot)
    residuals <- matrix(as.numeric(NA), nrow = nrow(response), ncol = ncol(response))
    rownames(residuals) <- rownames(response)
    colnames(residuals) <- colnames(response)

    if (is.null(grouping)) {
        grouping <- rep(1, nrow(response))
    }

    for (group in unique(grouping)) {
        selector <- grouping == group

        if (rescale.group) {
            cat("[INFO] Rescaling within ", group, ". Set rescale.group=F if this is not what you want\n")
            response.cur <- fast_colscale(response[selector, ])
        } else {
            response.cur <- response[selector, ]
        }

        for (fgroup in names(covariates.group)) {
            data <- getDataByObject(object, covariates.group[[fgroup]], assay.covar, assay.image, slot.covar, drop = F)

            has.na <- rowSums(is.na(data)) != 0

            if (na.rm) {
                # Optionally remove NAs from the design matrix    
                if (sum(has.na) > 0 ) {
                    warning(paste0("Dropping ", sum(has.na), " cells from design matrix because of NAs in: ", fgroup))
                }
                selector.final <- selector & !has.na
            } else {
                
                if (sum(has.na) > 0 ) {
                    stop("Detected NA's in design matrix. Either set na.rm=T or remove them manually.")
                }                
                selector.final <- selector
            }

            covariates.dont.use.cur <- check_unused_covar(data, covariates.dont.use)
            design                  <- model.matrix(~., data = data[selector.final, , drop = FALSE])
            group.features          <- feature.colnames[[fgroup]]
            
            # Remove the features being regressed so its not regressing features against themselves
            group.features <- group.features[!group.features %in% covariates.group[[fgroup]]]
            
            if (length(group.features) == 0) {
                warning(paste0("Skipped group, found no features for group", group, " and pattern ", fgroup))
                next()
            }
            
            cat("[INFO] Running object group: ", group, " for ", fgroup, " and selected ", length(group.features), " features\n")

            res <- lm_matrix(response.cur[selector.final, group.features],
                design,
                covariates.dont.use = covariates.dont.use.cur,
                residuals.only = TRUE
            )

            residuals[selector.final, group.features] <- res
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
#' @param slot The slot to use for regressing against. Can be "data" or "scale.data"
#' @param covariates Character vector of independent variables to use in the model
#' @param formula The formula to use for regression, string or formula. Defaults to additive model. See details
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param assay.covar The assay to grab covariates from. Defaults to assay argument
#' @param slot.covar The slot to grab covariates from. Can be "data" or "scale.data"
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param covariates.dont.use Only use if you understand the implications. See detaills
#' @param rescale.group When grouping is active, should the group be re-centered and scaled prior to regressing
#' @param ... Remaining arguments passed to [lm_matrix()]. Use with caution.
#'
#' @details
#' `grouping`
#'
#'  Runs seperate regressions per group in the indicated grouping variable. length(grouping) must be nrow(assay)
#'
#' `formula`
#'
#'  If NULL an additive model of all covariates is performed. Otherwise should be a string interpretable by [stats::formula()]
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
calculate_lm <- function(object, assay, slot, covariates, formula = NULL, grouping = NULL, assay.covar = NULL, slot.covar = NULL, assay.image = NULL, covariates.dont.use = NULL, rescale.group = FALSE, ...) {
    check_dataset_assay_slot(object, assay, slot)

    if (is.null(slot.covar)) {
        slot.covar <- slot
    }

    if (is.null(assay.covar)) {
        assay.covar <- assay
    }

    if (!is.null(grouping)) {
        if (length(grouping) != nrow(object)) {
        stop("grouping must be a vector of nrow(object)")
        }
    }

    data <- getDataByObject(object, covariates, assay.covar, assay.image, slot.covar, drop = F)

    if (ncol(data) <= 0) {
        stop("data (covariates) cannot be empty. Are your collumn names correct?")
    }

    covariates.dont.use <- check_unused_covar(data, covariates.dont.use)

    if (is.null(formula)) {
        design <- model.matrix(~., data = data)
    } else {
        
        if (is.character(formula)) {
            formula <- as.formula(formula)
            warning(paste0("Formula is character, converting to formula object: ", paste0(as.character(formula), collapse=" ")))
        }
        
        design <- model.matrix(formula, data = data)
    }
    response <- slot(object@assays[[assay]], slot)

    # Remove NA's from the design matrix
    rows.to.keep <- rowSums(is.na(data))==0
    
    if (sum(rows.to.keep) != nrow(design)) {
        warning(paste0(sum(!rows.to.keep), " NA's in the design matrix, removing these"))
    }
    
    response <- response[rows.to.keep,]
    grouping <- grouping[rows.to.keep]
    
    if (nrow(design) != nrow(response)) {
        stop("nrow(design) must equal nrow(assay)")
    }
    
    if (is.null(grouping)) {
        res <- lm_matrix(response, design, covariates.dont.use = covariates.dont.use, ...)
        return(res)
    } else {
        
        
        if (length(grouping) != nrow(response)) {
            stop("grouping must have the same length as nrow(assay)")
        }
        results <- list()
        for (group in unique(grouping)) {
            cat("[INFO] Starting regressions for group: ", group, "\n")
            selector <- grouping == group

            if (rescale.group) {
                cat("[INFO] Rescaling within ", group, ". Set rescale.group=F if this is not what you want\n")
                response.cur <- fast_colscale(response[selector, ])
            } else {
                response.cur <- response[selector, ]
            }

            results[[group]] <- lm_matrix(response.cur, design[selector, ], covariates.dont.use = covariates.dont.use, ...)
        }
        return(results)
    }
}

#-------------------------------------------------------------------------------
#' Check covariates.dont.use and expand factor levels
#'
#' @keywords internal
check_unused_covar <- function(data, covariates.dont.use) {
    if (!is.null(covariates.dont.use)) {
        tmp <- c()
        for (covar in covariates.dont.use) {
            if (!covar %in% colnames(data)) {
                tmp <- c(tmp, covar)
                next()
            }

            if (is.character(data[, covar])) {
                tmp <- c(tmp, paste0(covar, unique(data[, covar])))
            } else if (is.factor(data[, covar])) {
                tmp <- c(tmp, paste0(covar, levels(data[, covar])))
            } else {
                tmp <- c(tmp, covar)
            }
        }
        covariates.dont.use <- tmp

        cat("[INFO] Not using beta's of following covariates when calculating residuals: ", covariates.dont.use, "\n")
    }

    return(covariates.dont.use)
}


#-------------------------------------------------------------------------------
#' Calculate linear coefficients using a linear mixed model
#'
#' @description Fit a linear mixed model using REML / MLE (lme4) and find coefficients
#' @param object A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for regressing against. Can be "data" or "scale.data"
#' @param covariates Character vector of independent variables to use in the model
#' @param formula The formula to use for regression. Defaults to additive model. See details
#' @param grouping Vector with grouping variable if residuals be calculated per group of objects. See details
#' @param assay.covar The assay to grab covariates from. Defaults to assay argument
#' @param slot.covar The slot to grab covariates from. Can be "data" or "scale.data"
#' @param assay.image The image assay to use for grabbing covariates, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param rescale.group When grouping is active, should the subset be re-centered and scaled prior to regressing
#' @param ... Remaining arguments passed to [tglowr::lmm_matrix()] and then to [lmerTest::lmer()]. See detaills
#' @details
#'
#' Strongly reccomend reading extra parameters available in [tglowr::lmm_matrix()]
#'
#' `grouping`
#'
#'  Runs seperate regressions per group in the indicated grouping variable. length(grouping) must be nrow(assay)
#'
#' `formula`
#'
#'  Should be a string interpretable by [stats::formula()] and [lmerTest::lmer()] but should NOT
#'  have response in the forumla. E.g. when regressing `y ~ x + (1|donor)` forumula should be `~ x + (1|donor)`
#'  Defaults to a regular additive linear model.
#'  To specify derrived terms you might need to add I(), for polynomials for example use either poly() or I(x^2)
#'
#' @returns A list of regression results. If grouping != NULL, there is one list per group
#' @export
calculate_lmm <- function(object, assay, slot, covariates, formula = NULL, grouping = NULL, assay.covar = NULL, slot.covar = NULL, assay.image = NULL, rescale.group = TRUE, ...) {
    check_dataset_assay_slot(object, assay, slot)

    if (is.null(slot.covar)) {
        slot.covar <- slot
    }

    if (is.null(assay.covar)) {
        assay.covar <- assay
    }

    data <- getDataByObject(object, covariates, assay.covar, assay.image, slot.covar, drop = F)

    if (ncol(data) <= 0) {
        stop("data (covariates) cannot be empty. Are your collumn names correct?")
    }

    if (is.null(formula)) {
        formula <- paste("~", paste(colnames(data), collapse = " + "))
    }
    
    response <- slot(object@assays[[assay]], slot)


    #response <- response[rows.to.keep,]
    #grouping <- grouping[rows.to.keep]
    #data     <- data[rows.to.keep]


    if (is.null(grouping)) {
        res <- lmm_matrix(response, data, formula = formula, ...)
        return(res)
    } else {
        if (length(grouping) != nrow(response)) {
            stop("grouping must have the same length as nrow(assay)")
        }
        results <- list()
        for (group in unique(grouping)) {
            cat("[INFO] Starting regressions for group: ", group, "\n")
            selector <- grouping == group

            if (rescale.group) {
                cat("[INFO] Rescaling within ", group, ". Set rescale.group=F if this is not what you want\n")
                response.cur <- fast_colscale(response[selector, ])
            } else {
                response.cur <- response[selector, ]
            }


            results[[group]] <- lmm_matrix(response.cur, data[selector, ], formula = formula, ...)
        }
        return(results)
    }
}


#' Calculate linear coefficients
#'
#' @description Fit a linear mixed model using OLS and find coefficients
#' @param response A matrix with response variables
#' @param design The common design matrix to regress response againsts
#' @param covariates.dont.use Only use if you understand the implications. See detaills
#' @param residuals.only Only return the residual matrix
#' @param return.residuals Return residual matrix in the output list. Defaults to T if residual.only = TRUE
#' @param keep.zerocol Keep columns in design matrix that sum to zero
#' @param tol Minimal eigenvalue to keep columns in the design matrix
#' @param calculate_ll Calculate the log likelihood of the model
#' @param use.qr Use the qr decomposition instead of the cholesky decomposition on the design matrix. See details
#'
#' @details
#' This function re-uses the decompostion of the design matrix between different responses, which makes it much faster compared to looped [lm()] calls.
#' It also has extra options such as ignoring covariates when calculating residuals.
#' 
#' `use.qr`
#' If TRUE, the qr decomposition is used to calculate the effectsizes. This avoids a matrix inversion that can be unstable with highly correlated designs.
#' Singularites are removed based on the provided tol and the absolute diagonal of the QR decomp.
#' It uses the same basic implementation as [lm()]. This is the reccomended way
#' 
#' If FALSE, the following is computed. In most cases, it procudes identical results, but there are cases where this is unstable.
#' Singularities are removed based on the eigenvalue of crossprod(design) and the provided tol.
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
lm_matrix <- function(response, design, covariates.dont.use = NULL, residuals.only = FALSE, return.residuals = FALSE, keep.zerocol = FALSE, tol = 1e-7, calculate_ll=TRUE, use.qr=TRUE) {
  if (residuals.only && !return.residuals) {
    return.residuals <- TRUE
  }
  
  if (!is.null(covariates.dont.use)) {
    if (sum(covariates.dont.use %in% colnames(design)) != length(design)) {
      warning("Not all covariates specified in covariates.dont.use found. Check the output carefully if all is expected. This can happen if covariates.dont.use contains factors")
    }
  }
  
  if (!keep.zerocol) {
    design.colsums <- colSums(design)
    if (any(design.colsums == 0)) {
      msg <- "Found columns which sum to zero, dropped these as usually they indicate an empty factor level\n"
      msg <- paste0(msg, "If you do want to keep them, set keep.zerocol = TRUE\n")
      msg <- paste0(msg, "Offending collumns: ", colnames(design)[design.colsums == 0])
      warning(msg)
      
      design <- design[, design.colsums != 0]
    }
  }
  
  # Pre-calculate b component on design matrix
  if (use.qr) {
    cat("[INFO] Using QR decomp\n")
    
    # Find near singlular values
    qr.decomp <- qr(design)
    ev        <- abs(diag(qr.decomp$qr))
    offenders <- qr.decomp$pivot[ev < tol]
    
    if (any(ev < tol)) {
      msg <- paste0("qr diag < ", tol, " detected, dropping these variables from the design matrix. Set tol=NULL to skip this check\n")
      msg <- paste0(msg, "Offending collumns: ", paste0(colnames(design)[offenders], collapse=", "))
      warning(msg)
      cat("[WARN] ", msg, "\n")
      design <- design[, qr.decomp$pivot[ev > tol]]
    }
    
    # Final QR on the removed values
    qr.decomp <- qr(design)
    
    # Used later for calculating SEs
    p1           <- 1:qr.decomp$rank
    R            <- chol2inv(qr.decomp$qr[p1, p1, drop = FALSE])
    b            <- R
  } else {
    cat("[INFO] Using cholesky inverse of crossprod(design)\n")
    # Check for near zero eigenvalues in the design matrix
    a <- crossprod(design)
    if (!is.null(tol)) {
      ev <- eigen(a, only.values = T)$values
      if (any(ev < tol)) {
        msg <- paste0("Eigenvalues < ", tol, " detected, dropping these variables from the design matrix. Set tol=NULL to skip this check\n")
        msg <- paste0(msg, "Offending collumns: ", paste0(colnames(design)[ev < tol], collapse=", "))
        warning(msg)
        cat("[WARN] ", msg, "\n")
        design <- design[, ev > tol]
        a <- crossprod(design)
      }
    }
    
    # Use chol2inv on the cholesky decomposition
    # to invert rather then solve as this is faster. This only works on
    # positive-definite matrices, but I think this should always be true in this case
    # (dont quote me on this). If not it throws an error so you can use solve instead
    # Update: Turns out there are rare cases where this procudes non-sensical results
    # due to inversion not being perfect. The QR approach matches R's implementation
    # and is to be preferred.
    b <- chol2inv(chol(a))
  }
  
  # Matrices to save model coefficients
  if (!residuals.only) {
    coef <- matrix(NA, nrow = ncol(response), ncol = sum(!colnames(design) %in% covariates.dont.use))
    se <- matrix(NA, nrow = ncol(response), ncol = sum(!colnames(design) %in% covariates.dont.use))
    model.stats <- matrix(NA, nrow = ncol(response), ncol=8)
    rownames(coef) <- colnames(response)
    colnames(coef) <- colnames(design)[!colnames(design) %in% covariates.dont.use]
    rownames(se) <- colnames(response)
    colnames(se) <- colnames(design)[!colnames(design) %in% covariates.dont.use]
    rownames(model.stats) <- colnames(response)
    colnames(model.stats) <- c("r2", "adj.r2", "f-stat", "p-value", "df", "rss", "tss", "ll")
    mse.vec <- rep(NA, ncol(response))
    names(mse.vec) <- colnames(response)
    
    if (!is.null(covariates.dont.use)) {
      warning("Specified covariates.dont.use while returning model stats, this is not reccomended unless you understand the implications.")
    }
  }
  
  # Matrix to save residuals
  if (return.residuals) {
    residuals <- matrix(NA, nrow = nrow(response), ncol = ncol(response))
    rownames(residuals) <- rownames(response)
    colnames(residuals) <- colnames(response)
  }
  
  cat("[INFO] Starting regressions\n")
  
  pb <- progress::progress_bar$new(format = paste0("[INFO] Regressing [:bar] :current/:total (:percent) eta :eta"), total = ncol(response))
  
  for (col in seq_len(ncol(response))) {
    pb$tick()
    
    # Calculate effectiszes
    if (use.qr) {
      beta <- matrix(qr.coef(qr.decomp, response[, col]), ncol=1)
    } else {
      beta <- crossprod(b, crossprod(design, response[, col]))
    }
    
    # Residuals
    if (is.null(covariates.dont.use)) {
      beta.tmp   <- beta
      design.tmp <- design
    } else {
      # Include covariates for beta fitting, but don't correct for them
      beta.tmp   <- beta[!colnames(design) %in% covariates.dont.use, ]
      design.tmp <- design[, !colnames(design) %in% covariates.dont.use]
    }
    
    ypred <- (design.tmp %*% beta.tmp)
    rs    <- response[, col] - ypred
    
    # Optionally save residuals
    if (return.residuals) {
      residuals[, col] <- rs
    }
    
    # Optionally save model statistics
    if (!residuals.only) {
      if (!is.null(covariates.dont.use)) {
        # Make sure residuals are centered when calculating stats
        # This is not the case if dropping covariates
        rs <- rs - mean(rs)
      }
      
      # Residual and total sum of squares
      rss <- crossprod(rs)
      tss <- crossprod(response[, col] - mean(response[, col]))
      
      # RSS / df
      df           <- (length(rs) - ncol(design))
      mse          <- as.numeric(rss / df)
      mse.vec[col] <- mse
      
      # SE
      if (use.qr) {
        # This is how summary.lm calculates SEs
        se[col,]     <- sqrt(diag(R) * mse)[!colnames(design) %in% covariates.dont.use]
      } else {
        se[col, ]    <- as.numeric(sqrt(diag(mse * b)))[!colnames(design) %in% covariates.dont.use]
      }
      
      # Beta
      coef[col, ]  <- as.numeric(beta.tmp)
      
      # R2
      r2           <- 1 - (rss / tss)
      r2.adj       <- 1 - (1 - r2) * (length(rs) - 1) / df
      
      # Calculate F-statistic
      msr <- (tss - rss) / (ncol(design.tmp) - 1)
      f.stat <- msr / mse
      p <- 1 - pf(f.stat, ncol(design.tmp) - 1, df)
      
      # Calculate ll
      if (calculate_ll) {
        ll <- sum(dnorm(response[, col], mean = ypred, sd = sd(rs), log = TRUE))
      } else {
        ll <- NA
      }
      
      model.stats[col, ] <- c(r2, r2.adj, f.stat, p, df, rss, tss, ll)
    }
  }
  
  # Return only residual matrix
  if (residuals.only) {
    return(residuals)
  }
  
  # Return output list
  out.list <- list(coef = coef,
                   se = se,
                   pval = 2 * pt(-abs(coef/se), df = df),
                   model.stats = model.stats,
                   df = df,
                   df.m = ncol(design.tmp) - 1,
                   residuals = NULL,
                   mse = mse.vec,
                   cov.unscaled = b)
  
  class(out.list) <- "tglowlm"
  if (return.residuals) {
    out.list$residuals <- residuals
  }
  return(out.list)
}



#-------------------------------------------------------------------------------
#' Calculate linear coefficients using a Linear mixed model
#'
#' @description Fit a linear mixed model using REML / MLE (lme4) and find coefficients
#' @param response A matrix with response variables
#' @param design The matrix with predictor variables
#' @param formula The latter component of the formula as a string. I.e. `~ x + (1|donor)`
#' @param formula.null The null formula for a LRT as a string. The latter component of the formula. I.e. `~ x + (1|donor)`
#' @param residuals.only Only return the residual matrix
#' @param return.residuals Return residual matrix in the output list. Defaults to T if residual.only = TRUE
#' @param refit Refit the models using MLE during the LRT (anova)
#' @param ... Remaining parameters passed to [lmerTest::lmer()]
#'
#' @returns A list with regression results
#' @importFrom progress progress_bar
#' @importFrom lmerTest lmer
#' @importFrom performance model_performance
#' @importFrom lme4 lFormula
#' @export
lmm_matrix <- function(response, design, formula, formula.null = NULL, residuals.only = FALSE, return.residuals = FALSE, refit = FALSE, ...) {
    
    if (is(formula, "formula")) {
        formula <- paste0(as.character(formula), collapse=" ")
        warning(paste0("Provided formula is formula object, converted to character for parsing: ", formula))
    }
    
    if (is(formula.null, "formula")) {
        formula.null <- paste0(as.character(formula.null), collapse=" ")
        warning(paste0("Provided formula.null is formula object, converted to character for parsing: ", formula.null))
    }
    
    if (!is(formula, "character")) {
        stop("formula must be of class character")
    }

    if (!startsWith(formula, "~")) {
        stop("formula must be only the 2nd half, starting with ~")
    }

    if (!is.null(formula.null)) {
        if (!is(formula.null, "character")) {
            stop("formula.null must be of class character")
        }

        if (!startsWith(formula.null, "~")) {
            stop("formula.null must be only the 2nd half, starting with ~")
        }
    }

    if (residuals.only && !return.residuals) {
        return.residuals <- TRUE
    }

    # Matrices to save model coefficients
    if (!residuals.only) {
        tmp <- lme4::lFormula(paste0(colnames(response)[1], formula), data = cbind(response[, 1, drop = FALSE], design))
        n.coefs <- ncol(tmp$X)

        coef <- matrix(NA, nrow = ncol(response), ncol = n.coefs)
        rownames(coef) <- colnames(response)
        colnames(coef) <- colnames(tmp$X)
        se <- coef
        pval <- coef
        df <- coef

        if (!is.null(formula.null)) {
            model.stats <- matrix(NA, nrow = ncol(response), 7)
            rownames(model.stats) <- colnames(response)
            colnames(model.stats) <- c("r2_cond", "r2_marg","r2_cond_null", "r2_marg_null", "lrt_chisqr", "lrt_p-value", "lrt_df")
        } else {
            model.stats <- matrix(NA, nrow = ncol(response), 2)
            rownames(model.stats) <- colnames(response)
            colnames(model.stats) <- c("r2_cond", "r2_marg")
        }

    }

    # Matrix to save residuals
    if (return.residuals) {
        residuals <- matrix(NA, nrow = nrow(response), ncol = ncol(response))
        rownames(residuals) <- rownames(response)
        colnames(residuals) <- colnames(response)
    }

    cat("[INFO] Starting regressions\n")
    pb <- progress::progress_bar$new(format = paste0("[INFO] Regressing [:bar] :current/:total (:percent) eta :eta"), total = ncol(response))

    for (col in seq_len(ncol(response))) {
        pb$tick()

        data.cur <- cbind(response[, col, drop = F], design)
        form.cur <- as.formula(paste(colnames(data.cur)[1], formula))
        
        if (any(colSums(is.na(data.cur))==nrow(data.cur))) {
            warning(paste0("No non-NA cases, skipping col: ", col))
            next()
        }
        
        m <- lmerTest::lmer(form.cur, data = data.cur, ...)

        if (!is.null(formula.null)) {
            form.null.cur <- as.formula(paste(colnames(data.cur)[1], formula.null))
            m.null <- lmerTest::lmer(form.null.cur, data = data.cur, ...)
        }

        if (return.residuals) {
            residuals[, col] <- residuals(m)
        }

        if (!residuals.only) {
            tmp <- coefficients(summary(m))
            coef[col, ] <- as.numeric(tmp[, 1])
            se[col, ] <- as.numeric(tmp[, 2])
            df[col, ] <- as.numeric(tmp[, 3])
            pval[col, ] <- as.numeric(tmp[, 5])

            # Model R2
            perf <- performance::model_performance(m)
            model.stats[col, "r2_cond"] <- perf$R2_conditional
            model.stats[col, "r2_marg"] <- perf$R2_marginal

            if (!is.null(formula.null)) {
                lrt <- anova(m.null, m, refit = refit)
                perf.null <- performance::model_performance(m.null)
                model.stats[col, "r2_cond_null"] <- perf.null$R2_conditional
                model.stats[col, "r2_marg_null"] <- perf.null$R2_marginal
                
                model.stats[col, "lrt_chisqr"]   <- lrt$`Chisq`[2]
                model.stats[col, "lrt_p-value"]  <- lrt$`Pr(>Chisq)`[2]
                model.stats[col, "lrt_df"]       <- lrt$`Df`[2]
            }
        }
    }

    # Return only residual matrix
    if (residuals.only) {
        return(residuals)
    }

    # Return output list
    out.list <- list(coef = coef, se = se, pval = pval, df = df, model.stats = model.stats, residuals = NULL)
    class(out.list) <- "tglowlm"
    if (return.residuals) {
        out.list$residuals <- residuals
    }
    return(out.list)
}
