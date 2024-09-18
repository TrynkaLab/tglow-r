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
#'
#' @returns The mered output from \code{\link{tglow.read.fileset}}
#' @export
merge_filesets <- function(data) {
  if (class(data) != "list") {
    stop("Data argument must be a list.")
  }

  out <- list()

  ncol.cells <- as.numeric(lapply(data, function(x) {
    ncol(x[["cells"]])
  }))
  ncol.meta <- as.numeric(lapply(data, function(x) {
    ncol(x[["meta"]])
  }))
  ncol.orl <- as.numeric(lapply(data, function(x) {
    ncol(x[["orl"]])
  }))
  # ncol.nucl <- as.numeric(lapply(data, function(x){ncol(x[["nucl"]])}))

  selector <- (ncol.cells == as.numeric(names(which.max(table(ncol.cells))))) & (ncol.meta == as.numeric(names(which.max(table(ncol.meta))))) & (ncol.orl == as.numeric(names(which.max(table(ncol.orl)))))

  if (sum(selector) != length(selector)) {
    msg <- paste0(
      "Not all filesets have the same collumn number, dropping the ones with least frequent number. Retained ",
      sum(selector), "/", length(selector), " filesets. \n",
      "The following filesets are at issue:\n"
    )

    for (i in seq_along(data)[!selector]) {
      msg <- paste0(msg, i, ", ")
    }

    warning(msg)
  }

  out$cells <- dplyr::bind_rows(lapply(data[selector], function(x) {
    x[["cells"]]
  }), )
  out$meta <- dplyr::bind_rows(lapply(data[selector], function(x) {
    x[["meta"]]
  }), )
  out$orl <- dplyr::bind_rows(lapply(data[selector], function(x) {
    x[["orl"]]
  }), )
  # out$nucl <- dplyr::bind_rows(lapply(data[selector], function(x){x[["nucl"]]}),)

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
