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

  # Chanded this from dplyr::bind_rows to do.call(rbind())
  out$cells <- do.call(rbind, lapply(data[selector], function(x) {
    x[["cells"]]
  }), )
  out$meta <- do.call(rbind, lapply(data[selector], function(x) {
    x[["meta"]]
  }), )

  if (!skip.orl) {
    out$orl <- do.call(rbind, lapply(data[selector], function(x) {
      x[["orl"]]
    }), )
  } else {
    out$orl <- NULL
  }


  if ("features" %in% names(data[[1]])) {
    out$features <- do.call(rbind, lapply(data[selector], function(x) {
      x[["features"]]
    }))
  }

  if ("cells_norm" %in% names(data[[1]])) {
    out$cells_norm <- do.call(rbind, lapply(data[selector], function(x) {
      x[["cells_norm"]]
    }))
  }

  # Merge the child object matrices
  if ("children" %in% names(data[[1]])) {
    out$children <- list()
    for (obj in names(data[[1]]$children)) {
      out$children[[obj]] <- do.call(rbind, lapply(data[selector], function(x) {
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
#' @keywords internal
check_dataset_assay_slot <- function(dataset, assay, slot, assay.image = NULL) {
  # Checks for input
  if (!is.null(dataset)) {
    if (!is(dataset, "TglowDataset")) {
      stop("Dataset must be of class TglowDataset")
    }
  }

  if (!is.null(assay)) {
    if (assay %in% c(names(dataset@assays), c("image.data", "image.data.trans", "image.data.norm", "image.meta"))) {
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
    if (!assay.image %in% c("image.data", "image.data.trans", "image.data.norm", "image.meta")) {
      stop("Assay.image must be in c('image.data', 'image.data.trans', 'image.data.norm', 'image.meta')")
    }
  }
}



#-------------------------------------------------------------------------------
#' Check a list of TglowFilters for validitiy
#'
#' @description Check a list of TglowFilters for validitiy
#'
#' @keywords internal
check_filter_list <- function(filters) {
  for (filter in filters) {
    if (!is(filter, "TglowFilter")) {
      stop(paste0("Filter must be of class TglowFilter: ", filter))
    }
  }
}

#-------------------------------------------------------------------------------
#' Check if a package is available
#'
#' @keywords internal
check_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(paste0("Package '", package, "' is needed for this function to work. Please install it, also check bioconductor if not in CRAN"),
      call. = FALSE
    )
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
    #assay@data@.Data <- cbind(assay@data@.Data, features)
    assay@data <- TglowMatrix(cbind(assay@data, features))
    colnames(assay@data)[(ncol + 1):(ncol + ncol(features))] <- names

    if (!preserve.other && !is.null(assay@scale.data)) {
      #assay@scale.data@.Data <- cbind(assay@scale.data@.Data, matrix(NA, nrow(features), ncol(features)))
      assay@scale.data <- TglowMatrix(cbind(assay@scale.data, matrix(NA, nrow(features), ncol(features))))
      colnames(assay@scale.data)[(ncol + 1):(ncol + ncol(features))] <- names
    }
  }

  if (slot == "scale.data") {
    ncol <- ncol(assay@scale.data)
    #assay@scale.data@.Data <- cbind(assay@scale.data@.Data, features)
    assay@scale.data <- TglowMatrix(cbind(assay@scale.data, features))
    colnames(assay@scale.data)[(ncol + 1):(ncol + ncol(features))] <- names

    if (!preserve.other && !is.null(assay@data)) {
      #assay@data@.Data <- cbind(assay@data@.Data, matrix(NA, nrow(features), ncol(features)))
      assay@data <- TglowMatrix(cbind(assay@data, matrix(NA, nrow(features), ncol(features))))
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
#' Check if there is any overlap within a list of vectors
#'
#' @param list A list of vectors to overlap
#'
#' @returns logical indicating if any items overlap in the list
#' @export
list_has_overlap <- function(list) {
  n <- length(list)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (length(intersect(list[[i]], list[[j]])) > 0) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}


#-------------------------------------------------------------------------------
#' Utilities for making plate overviews
#' @rdname plate_utils
#' @export
new_384_plate <- function() {
  plate    <- matrix(NA, nrow=16, ncol=24)
  rownames(plate) <- c("A", "B", "C", "D", "E", "F", "G" ,"H", "I", "J", "K", "L", "M", "N", "O", "P")
  colnames(plate) <- 1:24
  return(plate)
}

#' @rdname plate_utils
#' @export
new_96_plate <- function() {
  plate    <- matrix(NA, nrow=8, ncol=12)
  rownames(plate) <- c("A", "B", "C", "D", "E", "F", "G" ,"H")
  colnames(plate) <- 1:12
  return(plate)
}

#' @rdname plate_utils
#' @export
well_to_index <- function(well) {
  row_index <- match(substr(well, 1, 1), LETTERS)
  col_index <- as.numeric(substr(well, 2, 3))
  return(list(row = row_index, col = col_index))
}




##-------------------------------------------------------------------------------
#' Match two objects together based on shared images and nearest neighbour in XY
#'
#' @description Adds features to an existing assay
#'
#' @param a A \linkS4class{TglowDataset} used as reference
#' @param b A \linkS4class{TglowDataset} used as query (added to reference)
#' @param tol Euclidian distance in xy space to consider the NN as valid
#' @param mode 'add' | 'merge', see details
#' @param assay.prefix The prefix to add to assay names when assays in a and b have the same name and mode = merge
#' @param meta.prefix The prefix to add to b@meta columns prior to merging, columns with the same name as a in b are dropped
#' 
#' @details 
#' In mode add all assays from b are added to a as new assays. In mode merge, assays with the same name
#' are merged together, only new features are added. If the @scale.data slot is missing in one of the
#' datasets this is set to NULL
#' 
#' NOTE: Image level data is not merged at this time
#' @returns The object a with the extra data from b
#' @export
match_objects_xy_nn <- function(a, b, tol=2, mode="add", assay.prefix="b_", meta.prefix=NULL) {
  
  if (is.null(a@feature.map)) {
    stop("Object a must have @feature.map set")
  }
  
  if (is.null(b@feature.map)) {
    stop("Object b must have @feature.map set")
  }
  
  
  df.a <- getDataByObject(a, c(a@feature.map@plate@feature,
                               a@feature.map@well@feature,
                               a@feature.map@field@feature,
                               a@feature.map@x@feature,
                               a@feature.map@y@feature),
                          assay=a@feature.map@x@assay,
                          slot=a@feature.map@x@slot)
  colnames(df.a) <- c("plate", "well", "field", "x", "y")
  
  df.b <- getDataByObject(b, c(b@feature.map@plate@feature,
                               b@feature.map@well@feature,
                               b@feature.map@field@feature,
                               b@feature.map@x@feature,
                               b@feature.map@y@feature),
                          assay=b@feature.map@x@assay,
                          slot=b@feature.map@x@slot)
  colnames(df.b) <- c("plate", "well", "field", "x", "y")
  
  if (any(is.na(df.a))) {
    stop("Object a has NA values in plate/well/field/x/y. This is not allowed. Please remove NA's before running")
    #df.a <- na.omit(df.a)
  }
  
  if (any(is.na(df.b))) {
    warning("Object b has NA values in plate/well/field/x/y This shouldn't normally happen and NA's are dropped")
    df.b <- na.omit(df.b)
  }
  
  
  # This is no longer needed as we used distance now
  #df.a$x <- trunc(df.a$x)
  #df.a$y <- trunc(df.a$y)
  #df.b$x <- trunc(df.b$x)
  #df.b$y <- trunc(df.b$y)

  df.a$pwf <- paste0(df.a$plate, ":", df.a$well, ":", df.a$field)
  df.b$pwf <- paste0(df.b$plate, ":", df.b$well, ":", df.b$field)
  
  df.a$idx <- 1:nrow(df.a)
  df.b$idx <- 1:nrow(df.b)
  
  # Use data table as it is much faster to subset then dataframe
  df.a <- data.table(df.a)
  df.b <- data.table(df.b)
  
  pb <- progress::progress_bar$new(format = "[INFO] Finding nearest neighbour [:bar] :current/:total (:percent) eta :eta", total = length(unique(df.a$pwf)))
  pb$tick(0)
  
  nearest.n <- matrix(NA, ncol=2, nrow=nrow(df.a))
  colnames(nearest.n) <- c("index_b", "dist")
  for (pwfc in intersect(unique(df.a$pwf), unique(df.b$pwf))) {
    pb$tick()
    cur.a <- df.a[pwf == pwfc]
    cur.b <- df.b[df.b$pwf == pwfc]
    
    cur.knn <- RANN::nn2(cur.b[,c("x", "y")], cur.a[,c("x", "y")], k=1)

    nearest.n[cur.a$idx, 1] <- cur.b[cur.knn$nn.idx[,1]]$idx
    nearest.n[cur.a$idx, 2] <- cur.knn$nn.dists[,1]
  }
  pb$terminate()

  nearest.n[nearest.n[,2] > tol,1] <- NA 

  perc <- (sum(!is.na(nearest.n[,1])) / nrow(df.a))*100
  cat("[INFO] ", round(perc, digits=2), "% of objects in object a matched to b at tolerance of ", tol, " distance\n")
  
  b.assays <- names(b@assays)
  
  # Merge assays with the same name
  if (mode == "merge") {
    for (assay in names(a@assays)) {
      if (assay %in% b.assays) {
        
        features.b <- setdiff(colnames(b@assays[[assay]]@data), colnames(a@assays[[assay]]@data))
        bb         <- b@assays[[assay]][nearest.n[,1],features.b]
        
        a@assays[[assay]]@data       <- cbind(a@assays[[assay]]@data, b@assays[[assay]]@data)
        
        # Merge feature DF
        featcols <- unique(c(colnames(a@assays[[assay]]@features), colnames(b@assays[[assay]]@features)))
        features <- data.frame(matrix(NA, nrow=nrow(a@assays[[assay]]@features)+length(features.b), ncol=length(featcols)))
        colnames(features) <- featcols
        rownames(features) <- c(rownames(a@assays[[assay]]@features), features.b)
        features[rownames(a@assays[[assay]]@features), colnames(a@assays[[assay]]@features)] <- a@assays[[assay]]@features
        features[rownames(b@assays[[assay]]@features), colnames(b@assays[[assay]]@features)] <- b@assays[[assay]]@features[features.b,]
        a@assays[[assay]]@features   <- features
        
        
        if (!is.null(a@assays[[assay]]@scale.data)) {
          if (!is.null(b@assays[[assay]]@scale.data)) {
            a@assays[[assay]]@scale.data <- cbind(a@assays[[assay]]@scale.data, b@assays[[assay]]@scale.data)
          } else {
            warning(paste0("@scale.data on object b and assay ", assay, " was null, dropping @scale.data on output"))
            a@assays[[assay]]@scale.data <- NULL
          }
        }
        
        b.assays <- b.assays[!b.assays %in% c(assay)]
      }
      
    } 
  } 
    
  # Add whatever wasn't merged as a new assay
  for (assay in b.assays) {
    if (assay %in% names(a@assays)) {
      new <- paste0(assay.prefix, assay)
    } else {
      new <- assay
    }
    a@assays[[new]] <- b@assays[[assay]][nearest.n[,1],]
    
    # Update the object IDs to the ones in a
    objectIds(a@assays[[new]]) <- objectIds(a)
  }
  
  
  if (!is.null(meta.prefix)) {
    colnames(b@meta) <- paste0(meta.prefix, colnames(b@meta))
  }
  
  meta.names <- colnames(b@meta)[!colnames(b@meta) %in% colnames(a@meta)]
  a@meta <- cbind(a@meta, b@meta[nearest.n[,1], meta.names])
    
  return(a)
    
}
