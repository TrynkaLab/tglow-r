#-------------------------------------------------------------------------------
#' Create TglowFilters from a filter table
#'
#' @description See  \linkS4class{TglowFilter} for detaills
#'
#' @param filter.table The table with filters
#' @param name.col Column id with filter name
#' @param col.col Column id with collumn pattern to apply filter too
#' @param func.col Column id with name of filter functions
#' @param thresh.col Column id with col with threshold
#' @param trans.col Column id with col indicating if transpose is applied before running func
#' @param active.col Column id with boolean setting filter to active or inactive
#' @returns A list of \linkS4class{TglowFilter}
#' @export
tglow_filters_from_table <- function(filter.table, name.col = 1, col.col = 2, func.col = 3, thresh.col = 4, trans.col = NULL, active.col = NULL) {
    filters <- list()

    for (row in seq(1, nrow(filter.table))) {
        f <- new("TglowFilter",
            name = filter.table[row, name.col],
            column_pattern = filter.table[row, col.col],
            func = filter.table[row, func.col],
            threshold = filter.table[row, thresh.col],
            transpose = ifelse(is.null(trans.col), FALSE, filter.table[row, trans.col]),
            active = ifelse(is.null(active.col), TRUE, filter.table[row, active.col])
        )
        filters[[row]] <- f
    }

    return(filters)
}


#-------------------------------------------------------------------------------
#' Find which feature filters to a tglow dataset
#'
#' @description
#' Calculate a feature filter table on a \linkS4class{TglowDataset}
#' All filters are inclusive
#'
#' @param dataset A tglow dataset
#' @param filters A list with TglowFilter objects to apply
#' @param assay The assay to use for calculating filters
#' @param slot The slot to use for calculating filters: "data" or "scale.data"
#' @param features An optional subset of features to use for calculation
#' @param na.fail Should NA be treated as fail, defaults to yes
#'
#' @returns A logical matrix of ncol(dataset[[assay]]) x length(filters) where T indicates filter pass and F indicates filter fail
#' @importFrom progress progress_bar
#' @export
calculate_feature_filters <- function(dataset, filters, assay, slot, features = NULL, na.fail = TRUE) {
    # Check inputs
    check_dataset_assay_slot(dataset, assay, slot)
    check_filter_list(filters)

    if (is.null(features)) {
        features <- colnames(slot(dataset@assays[[assay]], slot))
    }

    data <- slot(dataset[[assay]], slot)[, features]

    res <- matrix(TRUE, nrow = ncol(data), ncol = length(filters))
    rownames(res) <- colnames(data)
    colnames(res) <- sapply(filters, function(x) {
        x@name
    })

    for (i in seq_along(filters)) {
        filter <- filters[[i]]
        if (!filter@active) {
            cat("[INFO] Skipping, ", filter@name, ": ", filter@func, " as its not active\n")
            next()
        }

        cat("[INFO] ", filter@name, ": ", filter@func, "\n")

        # Select the features to apply to
        if (filter@column_pattern == "all") {
            cur.features <- colnames(data)
        } else {
            cur.features <- grep(filter@column_pattern, colnames(data), value = T)
        }

        cat("[INFO] Applying pattern: ", filter@column_pattern, " and selected: ", length(cur.features), " features \n")

        pb <- progress_bar$new(format = paste0("[INFO] ", filter@name, " [:bar] :current/:total (:percent) eta :eta"), total = length(cur.features))
        res[cur.features, filter@name] <- apply(data[, cur.features, drop=FALSE], 2, function(x) {
            pb$tick()
            if (is.na(filter@threshold) || is.null(filter@threshold)) {
                return(do.call(filter@func, list(vec = x)))
            } else {
                return(do.call(filter@func, list(vec = x, thresh = filter@threshold)))
            }
        })
    }

    if (na.fail) {
        res[is.na(res)] <- FALSE
    }

    return(res)
}


#-------------------------------------------------------------------------------
#' Apply a cell level filter to a tglow dataset
#'
#' @description
#' Calculate an object filter table on a \linkS4class{TglowDataset}
#' All filters are inclusive
#'
#' @param dataset A tglow dataset
#' @param filters A list with TglowFilter objects to apply
#' @param assay The assay to use for calculating filters
#' @param slot The slot to use for calculating filters: "data" or "scale.data"
#' @param grouping A vector specifying the grouping for some filters
#' @param features An optional subset of features to use for calculation
#' @param na.fail Should NA be treated as fail, defaults to yes
#'
#' @returns A logical matrix of nrow(dataset) x length(filters) where T indicates filter pass and F indicates filter fail
#' @export
calculate_object_filters <- function(dataset, filters, assay, slot = "data", grouping = NULL, features = NULL, na.fail=TRUE) {
    # Check inputs
    check_dataset_assay_slot(dataset, assay, slot)
    check_filter_list(filters)

    if (is.null(features)) {
        features <- colnames(slot(dataset@assays[[assay]], slot))
    }

    data <- slot(dataset[[assay]], slot)[, features]

    res <- matrix(TRUE, nrow = nrow(data), ncol = length(filters))
    rownames(res) <- rownames(data)
    colnames(res) <- sapply(filters, function(x) {
        x@name
    })

    for (i in seq_along(filters)) {
        cur.filter <- filters[[i]]

        cat("[INFO] ", cur.filter@name, ": ", cur.filter@func, "\n")

        # Select the features to apply to
        if (cur.filter@column_pattern == "all") {
            cur.features <- colnames(data)
        } else {
            cur.features <- grep(cur.filter@column_pattern, colnames(data), value = T)
        }

        cat("[INFO] Applying pattern: ", cur.filter@column_pattern, " and selected: ", length(cur.features), " features \n")

        cur.data <- data[, cur.features, drop = F]
        if (cur.filter@transpose) {
            cur.data <- t(cur.data)
        }

        if (is.na(cur.filter@threshold) || is.null(cur.filter@threshold)) {
            res[, cur.filter@name] <- do.call(cur.filter@func, list(
                vec = cur.data,
                grouping = grouping
            ))
        } else {
            res[, cur.filter@name] <- do.call(cur.filter@func, list(
                vec = cur.data,
                thresh = cur.filter@threshold,
                grouping = grouping
            ))
        }
    }

    if (na.fail) {
        res[is.na(res)] <- FALSE
    }

    return(res)
}


#-------------------------------------------------------------------------------
#' Apply feature filters to a tglow dataset
#'
#' @description
#' Apply a logical matrix or a logical array with filters on the columns and features on the rows
#' to a tglow dataset. Each assay is seperately subsetted to remove features
#' which are not selected
#'
#' If filter.res is a matrix, all rows must be true
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param filter.res The matrix of filter output from \code{\link{calculate_feature_filters}} or a named logical vector
#' @param assays Only apply filteres to these assays
#'
#' @returns The filtered \linkS4class{TglowDataset}
#' @export
apply_feature_filters <- function(dataset, filter.res, assays = NULL) {
    # Check inputs
    check_dataset_assay_slot(dataset, NULL, NULL)

    # Can be a single logical as well
    if (is(filter.res, "logical")) {
        if (is.null(names(filter.res))) {
            stop("If providing logical it must have names")
        }

        rn <- names(filter.res)
        filter.res <- matrix(filter.res, ncol = 1)
        rownames(filter.res) <- rn
    }

    if (!is(filter.res, "matrix")) {
        stop("filter.res must be of class matrix")
    }

    if (is.null(rownames(filter.res))) {
        stop("filter.res must have feature names as rownames")
    }

    # NA's are considered false
    filter.res[is.na(filter.res)] <- F
    keep <- rowSums(filter.res) == ncol(filter.res)

    if (is.null(assays)) {
        assays <- seq_along(dataset@assays)
    }

    for (assay in assays) {
        cur.keep <- colnames(dataset[[assay]])[colnames(dataset[[assay]]) %in% names(keep)[keep]]
        dataset@assays[[assay]] <- dataset[[assay]][, cur.keep]
    }

    return(dataset)
}


#-------------------------------------------------------------------------------
#' Filter a TglowDataset at the image level
#'
#' @description
#' Apply image level filter to a TglowDataset. All objects in that image are
#' removed. If filter.res is a matrix, all rows must be true to keep the image
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param filter.res A logical matrix or vector indicating which images to keep
#'
#' @returns The filtered \linkS4class{TglowDataset}
#' @export
apply_image_filters <- function(dataset, filter.res) {
    # Checks for input
    check_dataset_assay_slot(dataset, NULL, NULL)

    # Can be a single logical as well
    if (is(filter.res, "logical")) {
        filter.res <- matrix(filter.res, ncol = 1)
    }

    if (!is(filter.res, "matrix")) {
        stop("filter.res must be of class matrix")
    }

    if (nrow(filter.res) != nrow(dataset@image.meta)) {
        stop("filter.res must have the same number of rows as dataset@image.meta")
    }

    selector <- rowSums(filter.res) == ncol(filter.res)

    keep.imgs <- rownames(dataset@image.meta[selector, ])
    keep <- dataset@image.ids %in% keep.imgs
    dataset <- dataset[keep, ]

    return(dataset)
}
