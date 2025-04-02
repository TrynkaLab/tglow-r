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
#' This does not use existing dimension reductions to avoid issues when running with different QC groups
#' 
#' NA values are removed on a column basis and a warning message is raised. If you want
#' to handle them differenlty, could make a new assay where the NA's are handled manually
#'
#' 'pc.thresh', 'pc.max', 'pc.n'
#' By default, estimates the 0.25 * ncol(assay) PC's assuming this is enough PC's to get to 75%. If this is
#' not enough a warning is raised, and you can provide pc.max to override this behaviour or set use_irlba
#' to FALSE to compute all components using \code{\link[=prcomp]{prcomp()}}
#'
#' 
#' 'method'
#' When method 'z' or 'mod.z' PC's for a qc.group are scaled by either z-score or modified z-score respectively, and if a cell is an outlier in 
#' any PC it is considered an overall outlier.
#' 
#' When method is 'mahalanobis' an "approximate mahalanobis distance" is calculated.  The pc's passing pc.thresh are taken and the Euclidian distance
#' from each PC's center is calculated. Given cov(data) is positive definite, and all PC's are included, this is equivalent to mahalanobis distance on the data.
#' Otherwise they should be highly correlated, but your milage may vary and this may be dataset specific!
#' 
#' Once the distances are calculated, a p-value is derrived using the chi-sqr distiribution. This pvalue is by default FDR adjusted and any
#' records FDR < 0.05 considered outliers. If the `thresh` parameter is supplied, it applies to the RAW pvalues, not the FDR. FDR is only applied if thresh='auto'
#' Note, the pvalue adjustment is done over all QC groups, not per QC group to ensure the overall FPR is maintained.
#' 
#' @param dataset A tglow dataset
#' @param assay The assay to use
#' @param qc.group A vector indicating or column in dataset if PC's should be calculated in subgroups of the data. Default find outliers using all objects at once
#' @param thresh Threshold in absolute PC to consider an outlier. "auto" defaults to 3.5 for zscore/mad modes, and to bonferoni for mahalanobis 
#' @param pc.thresh The percentage of variance of PC's to select for outlier detection
#' @param pc.max The maximum number of components to calculate using \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}}
#' @param pc.n The number of PC's to use. Defaults to the number of PC's that reach pc.thresh or pc.max
#' @param method Method to scale PC's prior to selecting thresh. Value can be 'z' for z-score, 'mod.z' for modified zscore, 'mahalanobis' for Mahalanobis distance.
#' @param return.pcs Should the grouped PC's be returned?
#' @param use_irlba Logical if \code{\link[=prcomp_irlba]{irlba::prcomp_irlba()}} or \code{\link[=prcomp]{prcomp()}} should be used for PCA
#' @param features Features to include in PCA. By default uses all (NULL)
#' @param padj.method Adjustment to apply to outlier pvalues. Any accepted by [stats::p.adjust()]
#' @returns A list with:
#' - outliers: Boolean indicating outlier status
#' - dist: Distance metric, in case of 'z' or 'mod.z' the number of components that the object is an outlier in
#' - df: Degrees of freedom for chi-sqr (number of PC's used)
#' - pval: Pvalue of chiqr test of dist in case method='mahalanobis'
#' - pcs: In case return.pcs is true, a list of pca object with the pc's for each qc.group
#' 
#' @importFrom irlba prcomp_irlba
#' @export
find_outliers_pca <- function(dataset,
                              assay,
                              qc.group = NULL,
                              thresh = "auto",
                              pc.thresh = 0.75,
                              pc.max = NULL,
                              pc.n = NULL,
                              slot = "data",
                              method = "mahalanobis",
                              return.pcs = FALSE,
                              use_irlba = TRUE,
                              features = NULL,
                              padj.method="fdr") {
    # Check inputs
    check_dataset_assay_slot(dataset, assay, slot)


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

    #data <- slot(cur.assay, slot)@.Data
    data <- slot(cur.assay, slot)
    
    
    if (!is.null(features)) {
        data <- data[,features]
    }

    # Define results matrix
    final.outliers <- rep(NA, nrow(data))
    names(final.outliers) <- rownames(data)
    final.pca <- list()
    
    # Track the distances
    final.distances <- rep(NA, nrow(data))
    names(final.distances) <- rownames(data)
    
    # Track the dfs
    final.df <- rep(NA, nrow(data))
    names(final.df) <- rownames(data)
    
     # Track the pvals
    final.pval <- rep(NA, nrow(data))
    names(final.pval) <- rownames(data)
    
    
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

        # Rescale the data after subsetting
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
            if (thresh == "auto") {thresh <- 3.5}
            pcs.norm  <- apply(pca$x[, 1:pc.n.final, drop = F], 2, mod_zscore)        
            outliers  <- rowSums(abs(pcs.norm) < thresh) != pc.n.final
            distances <- rowSums(abs(pcs.norm) < thresh)
            pval      <- NA
        } else if (method == "z") {
            if (thresh == "auto") {thresh <- 3.5}
            pcs.norm  <- fast_colscale(pca$x[, 1:pc.n.final, drop = F])
            outliers  <- rowSums(abs(pcs.norm) < thresh) != pc.n.final
            distances <- rowSums(abs(pcs.norm) < thresh)
            pval      <- NA
            # pcs.norm <- apply(pca$x, 2, scale, center = T, scale = T)[, 1:pc.n, drop = F]
        } else if (method == "mahalanobis") {
            # Implement it here, you just want to set the vector "outliers" as a TRUE/FALSE logical
            # PCSs are in 'pca$x', data is in 'cur.data'
            #if (thresh == "auto") {thresh <- 0.05 / nrow(cur.data)}
            #cov.matrix <- diag(pca$sdev^2)
            #dist       <- mahalanobis(pca$x[, 1:pc.n.final, drop = F], center=F, diag(pca$sdev[1:pc.n.final]^2))
            
            # mahalanobis distance = euclidian distance in PCA space if data is centered and scaled
            eigenval   <- pca$sdev[1:pc.n.final]^2
            # Weigh the rotation by the eigenvalues
            tmp        <- pca$x[, 1:pc.n.final, drop = F] %*% diag(eigenval^-0.5)

            # Center each col and calc euclidian distance
            tmp         <- tmp - colMeans(tmp)
            dist        <- apply(tmp, 1, function(x){ sum(x^2)})
                        
            # Calculate pvalues on this 
            pval       <- pchisq(dist, df = pc.n.final, lower.tail = FALSE)
            distances  <- dist
            outliers   <- NA
            #if (thresh == "auto") {
            #    thresh <- 0.05
            #    padj   <- p.adjust(pval, method=padj.method)
            #} else {
            #    padj <- pval
            #}
            
            #outliers <- padj < thresh
            
            #stop("Method mahalanobis is not yet implemented")
        } else {
            stop(paste0("Method ", method, " is not available. Must be 'z', 'mod.z', 'mahalanobis'"))
        }

        if (return.pcs) {
            final.pca[[group]] <- list(pcs = pcs.norm, outliers = outliers)
        }

        # Add results to the outlier list
        final.outliers[rownames(cur.data)]  <- outliers
        final.distances[rownames(cur.data)] <- distances
        final.df[rownames(cur.data)]        <- pc.n.final
        final.pval[rownames(cur.data)]      <- pval
    }
    
    
    if (method == "mahalanobis" ){
        if (thresh == "auto") {
           thresh <- 0.05
           padj   <- p.adjust(final.pval, method=padj.method)
        } else {
           padj   <- final.pval
        }
        final.outliers <- padj < thresh
    }

    if (return.pcs) {
        return(list(outliers = final.outliers, dist=final.distances, df=final.df, pval=final.pval, pcs = final.pca))
    } else {
        return(list(outliers = final.outliers, dist=final.distances, df=final.df, pval=final.pval))
    }
}
