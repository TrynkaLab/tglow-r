#-------------------------------------------------------------------------------
#' Construct a feature metadata from cellprofiler feature names
#'
#' @description
#' Constructs a table extracting attributes from the names in the from
#' <object>_<category>_<feature-name>
#'
#' @param featue.names Character vector of feature names
#'
#' @returns A data.frame with feature attributes
#' @export
get_feature_meta_from_names <- function(feature.names) {
  feature.meta <- data.frame(
    id = feature.names,
    object = sapply(strsplit(feature.names, split = "_"), function(x) {
      x[[1]]
    }),
    measurement = sapply(strsplit(feature.names, split = "_"), function(x) {
      paste0(x[-1], collapse = "_")
    })
  )

  feature.meta$category <- sapply(strsplit(feature.meta$measurement, split = "_"), function(x) {
    x[[1]]
  })
  feature.meta$name <- sapply(strsplit(feature.meta$measurement, split = "_"), function(x) {
    paste0(x[-1], collapse = "_")
  })

  rownames(feature.meta) <- feature.meta$id
  return(feature.meta)
}

#-------------------------------------------------------------------------------
#' Merge list of filesets into one
#'
#' @description
#' Data must be a list of lists with outputs from tglow.read.fileset.a/b
#' or have the items, cells, meta, orl,  [children], [features], [cells_norm]
#'
#' @param data List output from \code{\link{tglow.read.fileset}}
#' @param skip.orl Should object relationships be skipped
#' @returns The mered output from \code{\link{tglow.read.fileset}}
#' @importFrom dplyr bind_rows
#' @export
merge_filesets <- function(data, skip.orl = FALSE) {
  if (!is(data, "list")) {
    stop("Data argument must be a list.")
  }

  out <- list()

  ncol.cells <- as.numeric(lapply(data, function(x) {
    ncol(x[["cells"]])
  }))
  ncol.meta <- as.numeric(lapply(data, function(x) {
    ncol(x[["meta"]])
  }))

  if (!skip.orl) {
    ncol.orl <- as.numeric(lapply(data, function(x) {
      ncol(x[["orl"]])
    }))
    selector <- (ncol.cells == as.numeric(names(which.max(table(ncol.cells))))) & (ncol.meta == as.numeric(names(which.max(table(ncol.meta))))) & (ncol.orl == as.numeric(names(which.max(table(ncol.orl)))))
  } else {
    selector <- (ncol.cells == as.numeric(names(which.max(table(ncol.cells))))) & (ncol.meta == as.numeric(names(which.max(table(ncol.meta)))))
  }

  if (sum(selector) != length(selector)) {
    msg <- paste0(
      "Not all filesets have the same collumn number, dropping the ones with least frequent number. Retained ",
      sum(selector), "/", length(selector), " filesets. \n",
      "The following filesets are at issue:\n"
    )

    for (i in seq_len(length(data))[!selector]) {
      msg <- paste0(msg, i, ", ")
    }

    cat(paste0("[WARN] ", msg, "\n"))
    warning(msg)
  }

  out$cells <- dplyr::bind_rows(lapply(data[selector], function(x) {
    x[["cells"]]
  }), )
  out$meta <- dplyr::bind_rows(lapply(data[selector], function(x) {
    x[["meta"]]
  }), )

  if (!skip.orl) {
    out$orl <- dplyr::bind_rows(lapply(data[selector], function(x) {
      x[["orl"]]
    }), )
  } else {
    out$orl <- NULL
  }


  if ("features" %in% names(data[[1]])) {
    out$features <- dplyr::bind_rows(lapply(data[selector], function(x) {
      x[["features"]]
    }))
  }

  if ("cells_norm" %in% names(data[[1]])) {
    out$cells_norm <- dplyr::bind_rows(lapply(data[selector], function(x) {
      x[["cells_norm"]]
    }))
  }

  # Merge the child object matrices
  if ("children" %in% names(data[[1]])) {
    out$children <- list()
    for (obj in names(data[[1]]$children)) {
      out$children[[obj]] <- dplyr::bind_rows(lapply(data[selector], function(x) {
        x[["children"]][[obj]]
      }))
    }
  }

  return(out)
}

#-------------------------------------------------------------------------------
#' Check a dataset, assay, slot
#'
#' @description Check a dataset, assay slot for validity or throw an error
#'
check_dataset_assay_slot <- function(dataset, assay, slot, assay.image = NULL) {
  # Checks for input
  if (!is.null(dataset)) {
    if (!is(dataset, "TglowDataset")) {
      stop("Dataset must be of class TglowDataset")
    }
  }

  if (!is.null(assay)) {
    if (assay %in% c(names(dataset@assays), c("image.data", "image.data.trans", "image.data.norm"))) {
      if (assay %in% c(names(dataset@assays))) {
        if (!is(dataset[[assay]], "TglowAssay")) {
          stop("Assay must be of class TglowAssay")
        }
      }
    } else {
      stop(paste0("Assay ", assay, " not available in dataset"))
    }
  }

  if (!is.null(slot)) {
    if (!slot %in% c("data", "scale.data")) {
      stop("Slot must be 'data' or 'scale.data'")
    }
  }

  if (!is.null(assay.image)) {
    if (!assay.image %in% c("image.data", "image.data.trans", "image.data.norm")) {
      stop("Assay.image must be in c('image.data', 'image.data.trans', 'image.data.norm')")
    }
  }
}



#-------------------------------------------------------------------------------
#' Check a list of TglowFilters for validitiy
#'
#' @description Check a list of TglowFilters for validitiy
#'
check_filter_list <- function(filters) {
  for (filter in filters) {
    if (!is(filter, "TglowFilter")) {
      stop(paste0("Filter must be of class TglowFilter: ", filter))
    }
  }
}

#-------------------------------------------------------------------------------
#' Find nearest index
#'
#' @description
#' Find the index of the value in x that is closest to value
#'
#' @param x A vector X
#' @param value The value to which X should be closest
#'
#' @returns The position in X where value is closest to
#' @export
nearest_index <- function(x, value) {
  which.min(abs(x - value))
}

#-------------------------------------------------------------------------------
#' Retrieve a cell (and its neighbours) based on a feature sumstat
#'
#' @description
#' Gets a objects and its closes neighbours based on a single feature and a
#' sumstat (mean, median, upper.q, lower.q). To customise the quantile used
#' specify q.
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param feature The feature to find a representative objects for
#' @param metric Can be 'mean', 'median', 'upper.q', 'lower.q'
#' @param na.rm Should NA's be removed
#' @param n How many objects either side of the representitative objects should be returend
#' @param subset Look only in a subset of objects. Must be a selection vector
#' @param q Override the quantile when using upper.q (0.75) or lower.q (0.25)
#'
#' @returns vector of indices in objects matrix
#' @export
fetch_representative_object <- function(dataset, assay, slot, feature, metric = "mean", na.rm = F, n = 0, subset = NULL, q = NULL) {
  check_dataset_assay_slot(dataset, assay, slot)

  x <- getDataByObject(dataset, feature, assay = assay, slot = slot, drop = T)
  # x <- slot(dataset@assays[[assay]], slot)[,feature]
  i <- seq_len(length(x))

  if (!is.null(subset)) {
    x <- x[subset]
    i <- i[subset]
  }

  i <- i[order(x)]
  x <- sort(x)

  if (metric == "mean") {
    m <- mean(x, na.rm = na.rm)
    out <- nearest_index(x, m)
  } else if (metric == "median") {
    m <- median(x, na.rm = na.rm)
    out <- nearest_index(x, m)
  } else if (metric == "upper.q") {
    if (is.null(q)) {
      q <- 0.75
    }
    m <- quantile(x, probs = q)
    out <- nearest_index(x, m)
  } else if (metric == "lower.q") {
    if (is.null(q)) {
      q <- 0.25
    }
    m <- quantile(x, probs = q)
    out <- nearest_index(x, m)
  }

  if (n > 0) {
    out <- (out - n):(out + n)
  }

  return(i[out])
}


#-------------------------------------------------------------------------------
#' Retrieve a set of objects (and its neighbours) based on a feature sumstat
#'
#' @description
#' Fetches (n*2)+1 objects arround the 0th, 10th, 25th, 50th, 75th, 90th and 100th quantiles
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param feature The feature to find a representative cell for
#' @param name Prefix to add to the names of the names vector
#' @param n How many objects either side of the representitative cell should be returend
#'
#' @returns A list with row ids and labels for the objects
#' @export
fetch_representative_object_quantiles <- function(dataset, assay, slot, feature, name = NULL, n = 1) {
  check_dataset_assay_slot(dataset, assay, slot)

  if (is.null(name)) {
    name <- feature
  }
  cur.assay <- slot(dataset@assays[[assay]], slot)@.Data

  f <- cur.assay[, feature]
  idx <- dataset@object.ids[!is.na(f)]
  f <- f[!is.na(f)]
  fs <- order(f)

  lower <- which(dataset@object.ids %in% idx[head(fs, n = (n * 2) + 1)])
  upper <- which(dataset@object.ids %in% idx[tail(fs, n = (n * 2) + 1)])

  l10 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.10,
    n = n
  )

  l25 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.25,
    n = n
  )

  l50 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.5,
    n = n
  )

  l75 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.75,
    n = n
  )

  l90 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.90,
    n = n
  )

  cn <- c(
    rep(paste0("q0 - ", name), (n * 2) + 1),
    rep(paste0("q10 - ", name), (n * 2) + 1),
    rep(paste0("q25 - ", name), (n * 2) + 1),
    rep(paste0("q50 - ", name), (n * 2) + 1),
    rep(paste0("q75 - ", name), (n * 2) + 1),
    rep(paste0("q90 - ", name), (n * 2) + 1),
    rep(paste0("q100 - ", name), (n * 2) + 1)
  )

  return(list(ids = c(lower, l10, l25, l50, l75, l90, upper), names = cn))
}



#-------------------------------------------------------------------------------
#' Add features to an existing assay
#'
#' @description Adds features to an existing assay
#'
#' @param assay A \linkS4class{TglowAssay}
#' @param slot The slot to add features to. The features are set to NA in the other slot unless preserve.other=TRUE
#' @param features A vector or matrix of features to append
#' @param names The feature names, default to colnames(features)
#' @param meta Feature metadata matching colnames(assay@features), each row is a feature when features is a matrix
#' @param preserve.other Instead of setting the other slot to NA, skip this step
#'
#' @returns The assay with the extra columns on slot
#' @export
add_features_to_assay <- function(assay, slot, features, names = NULL, meta = NULL, preserve.other = FALSE) {
  if (is.numeric(features) && !is.matrix(features)) {
    if (length(features) != nrow(assay)) {
      stop("Features must have the same length as assay")
    }

    if (is.null(names) && is.null(names(features))) {
      stop("Must provide names argument")
    }

    features <- matrix(features, ncol = 1)
  } else if (is.numeric(features) && is.matrix(features)) {
    if (nrow(features) != nrow(assay)) {
      stop("Features must have the same number of rows as assay")
    }
    if (is.null(names) && is.null(colnames(features))) {
      stop("Must either set colnames(features) or provide names argument")
    }

    if (is.null(names)) {
      names <- colnames(features)
    }
  } else {
    stop("Features must be numeric matrix or numeric vector")
  }

  if (slot == "data") {
    ncol <- ncol(assay@data)
    assay@data@.Data <- cbind(assay@data@.Data, features)
    colnames(assay@data)[(ncol + 1):(ncol + ncol(features))] <- names

    if (!preserve.other && !is.null(assay@scale.data)) {
      assay@scale.data@.Data <- cbind(assay@scale.data@.Data, matrix(NA, nrow(features), ncol(features)))
      colnames(assay@scale.data)[(ncol + 1):(ncol + ncol(features))] <- names
    }
  }

  if (slot == "scale.data") {
    ncol <- ncol(assay@scale.data)
    assay@scale.data@.Data <- cbind(assay@scale.data@.Data, features)
    colnames(assay@scale.data)[(ncol + 1):(ncol + ncol(features))] <- names

    if (!preserve.other && !is.null(assay@data)) {
      assay@data@.Data <- cbind(assay@data@.Data, matrix(NA, nrow(features), ncol(features)))
      colnames(assay@data)[(ncol + 1):(ncol + ncol(features))] <- names
    }
  }

  if (is.null(meta)) {
    meta <- matrix(NA, nrow = ncol(features), ncol = ncol(assay@features))
    meta[, 1] <- names
    colnames(meta) <- colnames(assay@features)
    rownames(meta) <- names
  } else {
    if (!all(colnames(meta) == colnames(assay@features))) {
      stop("Column names in meta must match those in assay@features")
    }
  }

  assay@features <- rbind(assay@features, meta)
  return(assay)
}


#-------------------------------------------------------------------------------
#' Add features to an existing assay
#'
#' @param list A list of vectors to overlap
#' 
#' @returns logical indicating if any items overlap in the list
#' @export
list_has_overlap <- function(list) {
  n <- length(list)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (length(intersect(list[[i]], list[[j]])) > 0) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}
