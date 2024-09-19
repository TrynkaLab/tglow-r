#-------------------------------------------------------------------------------
#' Find marker features using a t-test
#'
#' @description
#' Find marker features using a two sample t-test. By default top 10 positive results
#' for a class are returned.
#'
#' @param dataset A \linkS4class(TglowDataset)
#' @param ident A column in meta, image.meta, assay, assay.image to use as class labels
#' @param assay The assay to use
#' @param assay.image The image assay to use for grabbing ident, NULL, "image.data", "image.data.trans" or "image.data.norm"
#' @param slot The slot to use for calculating filters. Can be "data" or "scale.data"
#' @param return.top Return the top x values by t-statistic. Defaults to 10
#' @param ref.classes A vector of values to compare class mean against. Default NULL uses all other classes as reference
#'
#' @returns A data.frame with t-test results
#' @export
find_markers <- function(dataset, ident, assay, slot, assay.image = NULL, return.top = 10, ref.classes = NULL) {
    # Check input
    check_dataset_assay_slot(dataset, assay, slot)

    cur.assay <- slot(dataset[[assay]], slot)@.Data
    cur.ident <- as.character(getDataByObject(dataset, ident, assay = assay, assay.image = assay.image, slot = slot))

    if (is.null(cur.ident)) {
        stop("Ident not valid")
    }

    classes <- unique(cur.ident)
    res <- data.frame(matrix(NA, nrow = ncol(cur.assay) * length(classes), ncol = 9))
    colnames(res) <- c("class", "feature", "tstat", "df", "pval", "mean.diff", "mean.se", "mean.ref", "mean.class")
    i <- 1

    pb <- txtProgressBar(min = 0, max = ncol(cur.assay) * length(classes), style = 3)
    for (class in classes) {
        ident.is.class <- cur.ident == class

        if (is.null(ref.classes)) {
            ident.is.ref <- !ident.is.class
        } else {
            ident.is.ref <- (cur.ident != class) & (cur.ident %in% ref.classes)
        }

        for (col in colnames(cur.assay)) {
            setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)

            tmp <- t.test(cur.assay[ident.is.class, col], cur.assay[ident.is.ref, col])

            res[i, "class"] <- class
            res[i, "feature"] <- col
            res[i, "tstat"] <- tmp$statistic
            res[i, "df"] <- tmp$parameter
            res[i, "pval"] <- tmp$p.value
            res[i, "mean.diff"] <- tmp$estimate[1] - tmp$estimate[2]
            res[i, "mean.se"] <- tmp$stderr
            res[i, "mean.ref"] <- tmp$estimate[2]
            res[i, "mean.class"] <- tmp$estimate[1]

            i <- i + 1
        }
    }

    close(pb)

    if (!is.null(return.top)) {
        subsets <- list()
        for (class in classes) {
            tmp <- res[res$class == class, ]
            tmp <- tmp[order(tmp$tstat, decreasing = T), ]
            subsets[[class]] <- tmp[1:return.top, ]
        }

        res <- do.call(rbind, subsets)
    }

    return(res)
}
