#-------------------------------------------------------------------------------
#' Find outliers in PCA space
#'
#' @description Find outliers in the PCA space of a TglowAssay. Works for both
#' assays and image.data
#'
#' @details
#' To run on image data, just specify assay="image.data"|"image.data.trans"|"image.data.norm" which is implemented
#' as a special case
#'
#'
#' This does not use existing dimension reductions to avoid issues when running with different QC groups
#'
#'
#' NA values are removed on a column basis and a warning message is raised. If you want
#' to handle them differenlty, could make a new assay where the NA's are handled manually
#'
#' By default, estimates the 0.25 * ncol(assay) PC's assuming this is enough PC's to get to 75%. If this is
#' not enough a warning is raised, and you can provide pc.max to override this behaviour or set use_irlba
#' to FALSE to compute all components using \code{\link[=prcomp]{prcomp()}}
#'
#' @param dataset A tglow dataset
#' @param assay The assay to use
#' @param qc.group A vector indicating or column in dataset if PC's should be calculated in subgroups of the data. Default find outliers using all objects at once
#' @param thresh Threshold in absolute PC to consider an outlier
#' @param pc.thresh The percentage of variance of PC's to select for outlier detection
#' @param pc.max The maximum number of components to calculate using \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}}
#' @param pc.n The number of PC's to use. Defaults to the number of PC's that reach pc.thresh or pc.max
#' @param method Method to scale PC's prior to selecting thresh. Value can be 'z' for z-score or 'mod.z' for modified zscore
#' @param return.pcs Should the grouped PC's be returned?
#' @param use_irlba Logical if \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}} or \code{\link[=prcomp]{prcomp()}} should be used for PCA
#' @returns Logical indicating if object is an outlier in pca space, or a list if return.pcs=TRUE
#'
#' @importFrom irlba prcomp_irlba
#' @export
find_outliers_pca <- function(dataset,
                              assay,
                              qc.group = NULL,
                              thresh = 3.5,
                              pc.thresh = 0.75,
                              pc.max = NULL,
                              pc.n = NULL,
                              slot = "data",
                              method = "z",
                              return.pcs = FALSE,
                              use_irlba = TRUE) {
    # Check inputs
    tglowr:::check_dataset_assay_slot(dataset, assay, slot)


    if (is.null(qc.group)) {
        cat("[INFO] No QC group provided, setting all objects to same group\n")
        qc.group <- rep(1, nrow(data))
    }

    if (assay %in% c("image.data", "image.data.trans", "image.data.norm")) {
        cur.assay <- slot(dataset, assay)
        if (is.character(qc.group) && length(qc.group) == 1) {
            qc.group <- getImageData(dataset, qc.group, assay = assay, slot = slot)
        }
    } else {
        cur.assay <- dataset[[assay]]
        if (is.character(qc.group) && length(qc.group) == 1) {
            qc.group <- getDataByObject(dataset, qc.group, assay = assay, slot = slot)
        }
    }

    if (is.null(slot(cur.assay, slot))) {
        stop("Could not find data in slot, ", slot, "\n")
    }

    data <- slot(cur.assay, slot)@.Data

    # Define results matrix
    final.outliers <- rep(NA, nrow(data))
    names(final.outliers) <- rownames(data)
    final.pca <- list()

    # Now for each group
    unique.groups <- unique(qc.group)
    for (group in unique.groups) {
        cat("[INFO] Calculating outliers for for ", group, " ", which(group == unique.groups), "/", length(unique.groups), "\n")

        # Subset data
        cur.data <- data[qc.group == group, ]

        # Check if there enough cells
        if (sum(qc.group == group) < 3) {
            warning("Need at least two cells in group, skipping")
            final.outliers[qc.group == group] <- FALSE
            next
        }

        if (!is.null(pc.n)) {
            if (sum(qc.group == group) <= pc.n) {
                warning("Need more samples then n.pc, skipping")
                final.outliers[qc.group == group] <- FALSE
                next
            }
        }

        cur.data <- fast_colscale(cur.data, add_attr = F)

        # Remove features with ANY NA
        cur.data <- cur.data[, matrixStats::colSums2(is.na(cur.data)) == 0]
        n.rn.na <- ncol(data) - ncol(cur.data)
        if (ncol(cur.data) != ncol(data)) {
            warning(paste0("Removed ", n.rn.na, " features with NA's"))
        }

        # Remove features with zero variance
        cur.data <- cur.data[, matrixStats::colVars(cur.data) > 0]
        if (ncol(cur.data) != (ncol(data) - n.rn.na)) {
            warning(paste0("Removed ", ncol(data) - ncol(cur.data), " features with zero variance"))
        }

        # Again define PCs to less then the smallest dimension
        if (is.null(pc.max)) {
            pc.max <- round(ncol(data) * 0.25)
        }

        if (nrow(cur.data) < pc.max) {
            pc.final <- nrow(cur.data) - 1
        } else {
            pc.final <- pc.max
        }

        if (!is.null(pc.n)) {
            pc.final <- pc.n
        }

        # Define the number of PCs
        if (pc.final > ncol(cur.data)) {
            pc.final <- ncol(cur.data) - 1
        }

        if (!use_irlba) {
            pc.final <- ncol(cur.data)
        }

        cat("[INFO] Calculating  ", pc.final, " pc's on ", nrow(cur.data), " samples and, ", ncol(cur.data), " features\n")
        if (use_irlba) {
            # Calculate PCAs
            pca <- irlba::prcomp_irlba(cur.data, n = pc.final, center = F, scale = F)
            pc.var <- pca$sdev^2 / pca$totalvar
        } else {
            pca <- prcomp(cur.data, center = F, scale = F)
            pc.var <- pca$sdev^2 / sum(pca$sdev^2)
        }

        pc.n.final <- pc.n
        if (is.null(pc.n.final)) {
            if (sum(cumsum(pc.var) > pc.thresh) >= 1) {
                pc.n.final <- min(which(cumsum(pc.var) > pc.thresh))
            } else {
                pc.n.final <- pc.final
            }
        }

        if (cumsum(pc.var)[pc.n.final] < pc.thresh) {
            warning("Last pc doesnt pass pc.thresh, try increasing pc.max")
        }
        cat("[INFO] Selected ", pc.n.final, " pcs explaining ", round(cumsum(pc.var)[pc.n.final], digits = 2) * 100, "% of the variance\n")

        # Now we perform z-scoring on the PCs
        if (method == "mod.z") {
            pcs.norm <- apply(pca$x[, 1:pc.n.final, drop = F], 2, mod_zscore)
        } else if (method == "z") {
            pcs.norm <- fast_colscale(pca$x[, 1:pc.n.final, drop = F])
            # pcs.norm <- apply(pca$x, 2, scale, center = T, scale = T)[, 1:pc.n, drop = F]
        } else {
            stop(paste0("Method ", method, " is not available. Must be 'z' or 'mod.z'"))
        }

        outliers <- rowSums(abs(pcs.norm) < thresh) != pc.n.final

        if (return.pcs) {
            final.pca[[group]] <- list(pcs = pcs.norm, outliers = outliers)
        }

        # Add results to the outlier list
        final.outliers[rownames(cur.data)] <- outliers
    }

    if (return.pcs) {
        return(list(outliers = final.outliers, pcs = final.pca))
    } else {
        return(list(outliers = final.outliers))
    }
}
