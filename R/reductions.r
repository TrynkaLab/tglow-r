#-------------------------------------------------------------------------------
#' Calculate the PCA on an assay
#'
#' @description
#' Calculate a principal component analysis on a \linkS4class{TglowDataset} and add it to
#' the reductions. Does not re-scale by default and assumes scale.data
#' is scaled has mean 0 variance 1, which is not true if using modified zscore
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param pc.n How many PCs to calculate (if not NULL uses \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}} instead of \code{\link[=prcomp]{prcomp()}}
#' @param reduction.name The name to save the PCA results under. Defualts to PCA-<assay>
#' @param ret.prcomp Instead of returning just the PCs and variances, return the whole prcomp object
#' @param use_irlba Logical if \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}} or \code{\link[=prcomp]{prcomp()}} should be used for PCA
#' @param rescale Logical if matrix will be rescaled. This is advisable if using modified z-score

#' @returns \linkS4class{TglowDataset} with populated PCA slot
#' @export
calculate_pca <- function(dataset, assay, slot = "scale.data", pc.n = NULL, reduction.name = NULL, ret.prcomp = FALSE, use_irlba = FALSE, rescale = FALSE) {
    # Checks for input
    # check_dataset_assay_slot(dataset, assay, slot)
    if (is.null(slot(dataset[[assay]], slot))) {
        stop("Provided slot is NULL, have you set scale.data? (In future will implement that here)")
    }

    data <- slot(dataset@assays[[assay]], slot)

    if (rescale) {
        cat("[INFO] Rescaling data to mean 0 variance 1\n")
        #data <- fast_colscale(data@.Data)
        data <- fast_colscale(data)
    }

    # Remove features with ANY NA
    cur.data <- data[, matrixStats::colSums2(is.na(data)) == 0]
    n.rn.na <- ncol(data) - ncol(cur.data)
    if (ncol(cur.data) != ncol(data)) {
        warning(paste0("Removed ", n.rn.na, " features with NA's"))
    }

    if (!is.null(pc.n)) {
        # stop("pc.n not yet implemented #TODO")
        use_irlba <- TRUE
    }

    # PC's either use scale.data slot, which is centered and scaled,
    # or intentionally don't scale
    if (use_irlba) {
        # Calculate PCAs
        pcs <- irlba::prcomp_irlba(cur.data, n = pc.n, center = F, scale = F)
        rownames(pcs$x) <- dataset@object.ids

        res <- new("TglowReduction", x = pcs$x, var = pcs$sdev^2, var_total = pcs$totalvar)
    } else {
        pcs <- prcomp(cur.data, center = F, scale = F)
        rownames(pcs$x) <- dataset@object.ids

        res <- new("TglowReduction", x = pcs$x, var = pcs$sdev^2, var_total = sum(pcs$sdev^2))
    }

    if (ret.prcomp) {
        res@object <- pcs
        rownames(res@object$rotation) <- colnames(cur.data)
    }

    if (is.null(reduction.name)) {
        reduction.name <- paste0("PCA.", assay)
    }

    dataset@reduction[[reduction.name]] <- res

    return(dataset)
}

#-------------------------------------------------------------------------------
#' Calculate the UMAP on an assay
#'
#' @description
#' Calculate UMAP on a \linkS4class{TglowDataset} and add it to the reductions
#'
#' @details
#' Uses an existing reduction if possible. When downsampling objects it returns
#' a matrix in the reductions slot which has nrow(dataset) with objects which were
#' omitted set to NA to ensure downstream ordering is maintained
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param reduction The reduction to use for calculating UMAPs. If NULL and 'PCA.<assay>' is not available it is re-calculated
#' @param assay The assay to use to calculate PCA with, or to grab from reduction 'PCA.<assay>'
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param pc.n How many PC's to calculate
#' @param use_irlba Logical if \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}} or \code{\link[=prcomp]{prcomp()}} should be used for PCA
#' @param reduction.name The name to save the UMAP results under. Defualts to UMAP-<assay>
#' @param downsample Downsample to a random subset of objects prior to running UMAP. Can be an integer or a selection vector of row ids
#' @param ... Arguments passed to \code{\link{uwot::umap()}}
#'
#' @returns \linkS4class{TglowDataset} with populated reduction slot
#' @export
calculate_umap <- function(dataset, reduction = NULL, assay = NULL, slot = "scale.data", pc.n = 30, reduction.name = NULL, downsample = NULL, use_irlba = TRUE, ...) {
    # Check input
    check_dataset_assay_slot(dataset, assay, slot)

    if (!is.null(assay)) {
        if (paste0("PCA.", assay) %in% names(dataset@reduction)) {
            reduction <- paste0("PCA.", assay)
        }
    }

    if (is.null(reduction)) {
        if (is.null(assay)) {
            stop("Assay cannot be null when reduction is null")
        }
        cat("[INFO] Calculating principal components\n")
        dataset <- calculate_pca(dataset, assay = assay, slot = slot, pc.n = pc.n, use_irlba = use_irlba)
        reduction <- paste0("PCA.", assay)
    }

    # Fetch PCA
    pcs <- dataset@reduction[[reduction]]

    if (!is(pcs, "TglowReduction")) {
        stop("Reduction must be of class TglowReduction")
    }

    # Optionally downsample
    sample <- seq(0, nrow(dataset))
    if (!is.null(downsample)) {
        if (is(downsample, "numeric") && length(downsample) == 1) {
            sample <- sample(sample, downsample)
        } else if (is(downsample, "numeric")) {
            sample <- downsample
        } else {
            stop("downsample must be numeric vector with indices or an integer")
        }

        cat("[INFO] Downsampling to ", length(sample), " objects to run UMAP\n")
    }

    # Don't return the output directly from UMAP, but rather make sure when a downsample is
    # done the whole matrix is returned
    res <- matrix(NA, nrow = nrow(dataset), ncol = 2)
    rownames(res) <- dataset@object.ids
    umap <- uwot::umap(pcs@x[sample, 1:pc.n], ...)

    res[sample, ] <- umap

    if (is.null(reduction.name)) {
        reduction.name <- paste0("UMAP.", reduction)
    }
    colnames(res) <- c("UMAP-1", "UMAP-2")

    dataset@reduction[[reduction.name]] <- new("TglowReduction", x = res, object = list(sample = sample))

    return(dataset)
}
