#-------------------------------------------------------------------------------
#' Aggregate tglow objects
#'
#' @description Aggregate tglow objects using a grouping variable
#'
#' @details Aggregates objects using mean, median or sum based on a grouping variable
#' Strings are returned if there is a single unique value, otherwise they are dropped by default
#' Set drop.multival.string = FALSE to retain multival strings, which are pasted together using sep
#' drop.na.col removes any columns if drop.multival.string = TRUE and the result has only NA's left
#'
#' Both assays and image.data, image.data.trans and image.data.norm are aggregated
#'
#' @param input Input: either TglowMatrix, TglowAssay, data.frame or TglowDataset
#' @param grouping Grouping variable of nrow input. For aggregate_by_imagecol, grouping must be a column name on image.meta
#' @param method Aggregation method: "mean", "median", "sum"
#' @param group.order The row order of output
#' @param na.rm Should NA's be removed
#' @param drop.multival.string When concatenating character metadata, should strings with multiple values be set to NA. This is reccomended for large datasets
#' @param drop.na.col If after aggregating, a column is fully NA, should it be dropped?
#' @param sep Seperator string when concatenating multiple strings
#' @returns A logical where TRUE should be kept and FALSE values should be removed
#'
#' @rdname tglow_aggregate
#' @importFrom data.table as.data.table := data.table
#' @export
aggregate_tglow_matrix <- function(matrix, grouping, method = "mean", group.order = NULL, na.rm = TRUE) {
    if (!is(matrix, "TglowMatrix")) {
        stop("input must be TglowMatrix")
    }

    #dt <- data.table::as.data.table(matrix@.Data)
    dt <- data.table::as.data.table(matrix)
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
        # Re-scale to ensure mean 0 sd 1 of result
        new.scale.data <- TglowMatrix(fast_colscale(new.scale.data))
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

    result <- stats::aggregate(data, by = list(grouping), FUN = function(x, method) {
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
        } else if (is.logical(x)) {
            if (length(unique(x)) == 1) {
                return(as.logical(unique(x)))
            } else {
                return(NA)
            }
        }
    }, method = method)

    rownames(result) <- result[, 1]

    if (drop.na.col) {
        result <- result[, colSums(is.na(result)) != nrow(result)]
    }

    if (!is.null(group.order)) {
        return(result[group.order, ])
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
            group.order <- unique(tglowr::getImageData(object, grouping))
        }
        grouping.image <- tglowr::getImageData(object, grouping)
        grouping <- tglowr::getImageDataByObject(object, grouping)
    } else {
        stop("Grouping must be a single character indicating a column in image.meta")
    }

    if (class(grouping) == "logical") {
        warning(paste0("Specified grouping is logical. Converting to character"))
        grouping <- as.character(grouping)
    }

    if (class(grouping.image) == "logical") {
        warning(paste0("Specified image grouping is logical. Converting to character"))
        grouping.image <- as.character(grouping.image)
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

    new.meta       <- aggregate_metadata(object@meta, grouping, method, group.order, na.rm, drop.multival.string, drop.na.col, sep)
    new.image.meta <- aggregate_metadata(object@image.meta, grouping.image, method, group.order, na.rm, drop.multival.string, drop.na.col, sep)

    image.ids <- group.order
    names(image.ids) <- group.order

    # Put everything on a new TglowObject
    new.object <- new("TglowDataset",
        assays = new.assays,
        image.data = new.image.data,
        image.meta = new.image.meta,
        meta = new.meta,
        object.ids = group.order,
        image.ids = image.ids
    )
    
    names(new.object@object.ids)   <- new.object@object.ids

    return(new.object)
}


#-------------------------------------------------------------------------------
#' Convert a feature to a 96/384 plate layout by aggreation function
#' 
#' @param object A \linkS4class{TglowDataset}
#' @param assay The assay to grab the feature from, passed to \code{\link{getDataByObject}}
#' @param slot The slot to add features to. The features are set to NA in the other slot unless preserve.other=TRUE
#' @param feature A single feature accessible by getDataByObject
#' @param feature.well  \linkS4class{TglowFeatureLocation} describing the well feature or NULL (takes it from datasets@feature.map),
#' @param feature.plate  \linkS4class{TglowFeatureLocation} describing the plate feature or NULL (takes it from datasets@feature.map),
#' @param na.rm Should NA's be removed. Passed to method
#' @param format Output plate format, character either '384' or '96'
#' @param method A callable to use as aggeration function, mean, median, unique etc
#' 
#' @details 
#' If feature is a string, method=[base::unique()], which will work if the well only has unique values
#' 
#' @returns A list of plates with features aggeregates
#' @export
#'
aggregate_to_plate <- function(object, assay, slot, feature, feature.well=NULL, feature.plate=NULL, na.rm=T, format="384", method=base::mean) {

    # Grab features from the feature map
    if (is.null(feature.well) || is.null(feature.plate)) {
        if (is.null(object@feature.map)) {
            stop("Must either set object@feature.map or provide feature.well and feature.plate")
        }
        
        if (is.null(feature.well)) {
            feature.well <- object@feature.map@well
        }
        if (is.null(feature.plate)) {
            feature.plate <- object@feature.map@plate
        }   
    } 

    if ((class(feature.well) != "TglowFeatureLocation") || (class(feature.plate) != "TglowFeatureLocation")) {
        stop("feature.well/feature.plate must be of type TglowFeatureLocation")
    }

    if (assay %in% c("image.data", "image.meta", "image.data.norm", "image.data.norm")) {
        # Build the plot df
        data <- cbind(getImageData(object, feature, assay=assay, slot=slot, drop=F),
                    getImageData(object, feature.well@feature, assay=feature.well@assay, slot=feature.well@slot, drop=F),
                    getImageData(object, feature.plate@feature, assay=feature.plate@assay, slot=feature.plate@slot, drop=F))
    } else {
        # Build the plot df
        check_dataset_assay_slot(object, assay=assay, slot=slot)
        
        data <- cbind(getDataByObject(object, feature, assay=assay, slot=slot, drop=F),
                    getDataByObject(object, feature.well@feature, assay=feature.well@assay, slot=feature.well@slot, drop=F),
                    getDataByObject(object, feature.plate@feature, assay=feature.plate@assay, slot=feature.plate@slot, drop=F))
    }

    #data <- getDataByObject(object,j=c(feature.well, feature.plate, feature), assay=assay, slot=slot)

    if (class(data[, feature]) == "character") {
    warning("[WARN] Supplied feature is character, setting unique as default method but this might not work in case wells have multiple values")
    method = unique
    }

    plates <- list()
    for (plate in unique(data[,feature.plate@feature])) {
    cur.data <- data[data[,feature.plate@feature] == plate,]
    tmp <- aggregate(cur.data[,feature, drop=F], by=list(well=cur.data[,feature.well@feature]), FUN=method, na.rm=na.rm)

    if (format == "384") {
        cur_plate <- new_384_plate()
    } else if (format == "96") {
        cur_plate <- new_96_plate()
    } else {
        stop(paste0(format, " is not a valid plate format"))
    }

    for (row in 1:nrow(tmp)) {
        idx = well_to_index(tmp[row, 1])
        cur_plate[idx$row, idx$col] = tmp[row, 2]
    }

    plates[[plate]] <- cur_plate
    }

    return(plates)
}
