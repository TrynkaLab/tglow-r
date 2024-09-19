#-------------------------------------------------------------------------------
#' Aggregate tglow objects
#'
#' @description Aggregate tglow objects using a grouping variable
#'
#' @details Aggregates objects using mean, median or sum based on a grouping variable.
#' Strings are returned if there is a single unique value, otherwise they are dropped by default.
#' Set drop.multival.string = FALSE to retain multival strings, which are pasted together using sep.
#' drop.na.col removes any columns if drop.multival.string = TRUE and the result has only NA's left.
#'
#' Both assays and image.data, image.data.trans and image.data.norm are aggregated.
#'
#' @param input Input: either TglowMatrix, TglowAssay, data.frame or TglowDataset
#' @param grouping Grouping variable of nrow input. For aggregate_by_imagecol, grouping must be a column name on image.meta
#' @param method Aggregation method: "mean", "median", "sum"
#' @param group.order The row order of output
#' @param na.rm Should NA's be removed
#' @param drop.multival.string When concatenating character metadata, should strings with multiple values be set to NA. This is reccomended for large datasets.
#' @param drop.na.col If after aggregating, a column is fully NA, should it be dropped?
#' @param sep Seperator string when concatenating multiple strings
#' @returns A logical where TRUE should be kept and FALSE values should be removed.
#'
#' @rdname tglow_aggregate
#' @export
aggregate_tglow_matrix <- function(matrix, grouping, method = "mean", group.order = NULL, na.rm = TRUE) {
    if (!is(matrix, "TglowMatrix")) {
        stop("input must be TglowMatrix")
    }

    dt <- as.data.table(matrix@.Data)
    dt[, group := grouping]

    if (method == "mean") {
        result <- dt[, lapply(.SD, mean, na.rm = na.rm), by = group]
    } else if (method == "median") {
        result <- dt[, lapply(.SD, median, na.rm = na.rm), by = group]
    } else if (method == "sum") {
        result <- dt[, lapply(.SD, sum, na.rm = na.rm), by = group]
    } else {
        stop(paste0(method, " is not a valid method"))
    }


    result <- as.data.frame(result)
    rownames(result) <- result[, "group"]

    if (!is.null(group.order)) {
        result <- result[group.order, ]
    }

    return(as.matrix(result[, -which(colnames(result) == "group")]))
}


#-------------------------------------------------------------------------------
#' Aggregate a TglowAssay
#' @rdname tglow_aggregate
#' @export
aggregate_assay <- function(assay, grouping, method, group.order = NULL, na.rm = TRUE) {
    if (!is(assay, "TglowAssay")) {
        stop("input must be TglowAssay")
    }
    if (is.null(group.order)) {
        group.order <- unique(grouping)
    }

    new.data <- aggregate_tglow_matrix(assay@data, grouping, method, group.order, na.rm)
    new.data <- TglowMatrix(new.data)

    if (!is.null(assay@scale.data)) {
        new.scale.data <- aggregate_tglow_matrix(assay@scale.data, grouping, method, group.order, na.rm)
        new.scale.data <- TglowMatrix(new.scale.data)
    } else {
        new.scale.data <- NULL
    }

    new.assay <- new("TglowAssay",
        data = new.data,
        scale.data = new.scale.data,
        features = assay@features
    )

    return(new.assay)
}

#-------------------------------------------------------------------------------
#' Aggregate a a dataframe with metadata
#' @rdname tglow_aggregate
#' @export
aggregate_metadata <- function(data, grouping, method, group.order = NULL, na.rm = TRUE, drop.multival.string = TRUE, drop.na.col = TRUE, sep = "; ") {
    if (!is(data, "data.frame")) {
        stop("input must be data.frame")
    }


    result <- aggregate(data, by = list(grouping), FUN = function(x, method) {
        if (is.numeric(x)) {
            if (method == "mean") {
                return(mean(x, na.rm = na.rm))
            } else if (method == "median") {
                return(median(x, na.rm = na.rm))
            } else if (method == "sum") {
                return(sum(x, na.rm = na.rm))
            } else {
                stop(paste0(method, " is not a valid method"))
            }
        } else if (is.character(x)) {
            x <- unique(x)
            if (length(x) > 1) {
                if (!drop.multival.string) {
                    return(paste(x, collapse = sep))
                } else {
                    return(NA)
                }
            } else {
                return(x)
            }
        } else if (is.factor(x)) {
            x <- as.character(x)
            x <- paste(unique(x), collapse = sep)
            return(factor(x))
        }
    }, method = method)


    if (drop.na.col) {
        result <- result[, colSums(is.na(result)) != nrow(result)]
    }

    return(result)
}

#-------------------------------------------------------------------------------
#' Aggregate a TglowDataset by an image metadata column
#' @rdname tglow_aggregate
#' @export
aggregate_by_imagecol <- function(object, grouping, method, group.order = NULL, na.rm = TRUE, drop.multival.string = TRUE, drop.na.col = TRUE, sep = "; ") {
    if (!is(object, "TglowDataset")) {
        stop("Input must be TglowDataset")
    }


    if (!grouping %in% colnames(object@image.meta)) {
        stop("Grouping must be a single character indicating a column in image.meta")
    }

    if (length(grouping) == 1) {
        if (is.null(group.order)) {
            group.order <- unique(getImageData(object, grouping))
        }
        grouping.image <- getImageData(object, grouping)
        grouping <- getImageDataByObject(object, grouping)
    } else {
        stop("Grouping must be a single character indicating a column in image.meta")
    }


    # Aggreage assays
    new.assays <- list()
    for (assay in names(object@assays)) {
        new.assays[[assay]] <- aggregate_assay(object@assays[[assay]], grouping, method, group.order, na.rm)
    }

    # Assays and image data
    new.image.data <- aggregate_assay(object@image.data, grouping.image, method, group.order, na.rm)

    if (!is.null(object@image.data.trans)) {
        new.image.data.trans <- aggregate_assay(object@image.data.trans, grouping.image, method, group.order, na.rm)
    } else {
        new.image.data.trans <- NULL
    }

    if (!is.null(object@image.data.norm)) {
        new.image.data.norm <- aggregate_assay(object@image.data.norm, grouping.image, method, group.order, na.rm)
    } else {
        new.image.data.norm <- NULL
    }

    new.meta <- aggregate_metadata(object@meta, grouping, method, group.order, na.rm, drop.multival.string, drop.na.col, sep)
    new.image.meta <- aggregate_metadata(object@image.meta, grouping.image, method, group.order, na.rm, drop.multival.string, drop.na.col, sep)

    # Put everything on a new TglowObject
    new.object <- new("TglowDataset",
        assays = new.assays,
        image.data = new.image.data,
        image.meta = new.image.meta,
        meta = new.meta,
        object.ids = group.order,
        image.ids = group.order
    )

    return(new.object)
}
