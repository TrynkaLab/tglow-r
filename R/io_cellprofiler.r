#-------------------------------------------------------------------------------
# Imports
#' @import data.table
NULL

#-------------------------------------------------------------------------------
#' Read a cellprofiler fileset directory tree
#'
#' @description
#' Read a cellprofiler fileset directory tree organized into a <plate>/<well>/<field>.fileset
#' structure
#'
#' @param path path to tglow output dir
#' @param pattern The pattern that uniquely identfies a well. Use '.zip' for type 'B' and the '_experiment.tsv' for type 'A'
#' for type A
#' @param type Must be 'A' or 'B'. See details
#' @param n Read a subset of filesets. If integer, only that fileset is read, otherwise specify indices to read
#' @param skip.orl Skip reading of object relationships, as this can get quite large with many children and is not used
#' @param verbose Should I be chatty?
#' @param col.object The collumn name in the features which contains the per object object identifier. See details
#' @param col.meta.img.id The collumn name in the image level data which contains the image id. See details
#' @param ... Remaining parameters passed to \code{\link{read_cellprofiler_fileset_a}} or \code{\link{read_cellprofiler_fileset_b}}
#'
#' @details
#' 
#' `type`
#' Type A: _cells.tsv, _image.tsv, _experiment.tsv and _objectRelations.tsv
#' See \code{\link{read_cellprofiler_fileset_a}} for detaills
#'
#' Type B: <plate>_<well>.zip with individual files for each child object. Main object is assumed to be _cells
#' See \code{\link{read_cellprofiler_fileset_b}} for detaills
#' 
#' 
#' `col.object` and `col.meta.img.id`
#' 
#' See \code{\link{add_global_ids}} for details on how globally unique Id's are assigned in the case you need to set
#' `col.object` or `col.meta.img.id``.
#' @importFrom progress progress_bar
#' @export
read_cellprofiler_dir <- function(path, pattern, type, n = NULL, skip.orl = TRUE, verbose = F, col.object = "cell_ObjectNumber_Global", col.meta.img.id = "ImageNumber_Global", ...) {
  files <- list.files(path, recursive = T, pattern = paste0("*", pattern), full.names = T)

  if (type == "A") {
    prefixes <- gsub(pattern, "", files)
  } else {
    prefixes <- files
  }

  if (!is.null(n)) {
    prefixes <- prefixes[n]
  }

  # TODO: Replace this global with a function argument
  # Reset global fileset index
  assign("FILESET_ID", 0, envir = .GlobalEnv)

  # Read filesets
  filesets <- list()
  null.filesets <- 0
  pb <- progress::progress_bar$new(format = "[INFO] Reading [:bar] :current/:total (:percent) eta :eta", total = length(prefixes))
  pb$tick(0)
  for (pre in prefixes) {
    pb$tick()
    if (type == "A") {
      cur <- tglowr::read_cellprofiler_fileset_a(pre, return.feature.meta = F, skip.orl = skip.orl, ...)
    } else if (type == "B") {
      cur <- tglowr::read_cellprofiler_fileset_b(pre, return.feature.meta = F, skip.orl = skip.orl, verbose = verbose, ...)
    } else {
      stop(paste0("Invalid type: ", type))
    }

    if (!is.null(cur)) {
      filesets[[pre]] <- cur
      if (verbose) cat("\n[DEBUG] cols:", ncol(cur$cells), " cols meta", ncol(cur$meta), " cols orl:", ncol(cur$orl), "\n")
    } else {
      null.filesets <- null.filesets + 1
      # warning("Fileset was NULL, skipped.")
    }
  }

  if (null.filesets != 0) {
    msg <- paste0("Dectected ", null.filesets, "/", length(filesets), " as NULL (no cells)")
    warning(msg)
    cat(paste0("[WARN] ", msg, "\n"))
  }

  cat("\n[INFO] Merging filesets\n")
  output <- tglowr::merge_filesets(filesets, skip.orl = skip.orl)

  cat("[INFO] names: ", names(output), "\n")

  if (verbose) {
    cat("[DEBUG] colnames:\n", colnames(output$cells))
  }

  features <- tglowr::get_feature_meta_from_names(colnames(output$cells))
  classes <- sapply(output$cells, class)
  features$type <- classes[features$id]

  features <- features[colnames(output$cells), ]
  selector <- !is.na(output$cells[, col.object])

  if (sum(!selector) != 0) {
    warning(paste0("Detected ", sum(selector), " objects with NA in ", col.object, " removing these"))
  }

  output$cells <- output$cells[selector, ]
  rownames(output$cells) <- output$cells[, col.object]
  rownames(output$meta) <- output$meta[, col.meta.img.id]

  return(c(output, list(features = features)))
}

#-------------------------------------------------------------------------------
#' Add a global id to a matrix of image files
#' 
#' @description Take a matrix and extract patterns 'ImageNumber', 'ObjectNumber' and
#' 'Object_Number' and add a globably unique prefix. Will store these globally unique
#' id's in columns <original>_Global.
#' 
#' @param matrix An input matrix or data.frame
#' 
#' @details
#' Matrix with column pattern 'ImageNumber', 'ObjectNumber' and 'Object_Number'
#' which will be duplicated and have a globally unique variable added
#' Output of this is stored in the same column name but with suffix _Global
#'
#' @returns A data.frame with extra columns suffixed by _Global with globally unique ids
#' @export
add_global_ids <- function(matrix) {
  # Fetch global fileset id
  global.prefix <- paste0("FS", FILESET_ID)

  # Image numbers
  cols.i <- grep("ImageNumber", colnames(matrix), value = T)

  for (cur.col in cols.i) {
    selector <- !is.na(matrix[, cur.col])
    matrix[selector, paste0(cur.col, "_Global")] <- paste0(global.prefix, "_I", matrix[selector, cur.col])
  }

  # Object numbers
  cols <- grep("ObjectNumber", colnames(matrix), value = T)
  cols <- c(cols, grep("Object_Number", colnames(matrix), value = T))
  cols <- c(cols, grep("Parent", colnames(matrix), value = T))

  for (cur.col in cols) {
    selector <- !is.na(matrix[, cur.col])
    matrix[selector, paste0(cur.col, "_Global")] <- paste0(global.prefix, "_I", matrix[selector, cols.i[1]], "_O", matrix[selector, cur.col])
  }

  return(matrix)
}



#-------------------------------------------------------------------------------
#' Read a cell level fileset type A
#'
#' @description
#' Reads a CellProfiler fileset into a list
#' Type A: Assumes all features are in a single _cells.tsv / _cells.tsv
#' and all features are matched
#' 
#' @param prefix Path prefix to fileset
#' @param return.feature.meta Should the dataframe with feature metadata be added
#' @param add.global.id Should extra id columns be added that are globally unique
#' @param pat.img The suffix pattern to identify the image level data
#' @param pat.cells The suffix pattern to identify the cell level data
#' @param pat.orl The suffix pattern to identify object relationships
#' @param skip.orl Should object relationships be read (not used, can be quite large)
#'
#' @returns list with data frames:
#'
#' - cells (cell level features)
#'
#' - meta (image level features)
#'
#' - objectRelations
#'
#' - features (optional)
#'
#'
#' Output is NULL if no cells are detected
#'
#' @export
read_cellprofiler_fileset_a <- function(prefix,
                                        return.feature.meta = F,
                                        add.global.id = T,
                                        pat.img = "_image.tsv",
                                        pat.cells = "_cells.tsv",
                                        pat.orl = "_objectRelationships.tsv",
                                        skip.orl = FALSE) {
  if (add.global.id) {
    assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", FILESET_ID)
  }

  # Read header of _cells.tsv
  cells <- data.table::fread(paste0(prefix, pat.cells), data.table = F, nrows = 3, showProgress = FALSE)

  # If the file has no cells empty
  if (nrow(cells) != 3) {
    return(NULL)
  }

  # Clean colnames
  cn <- paste0(colnames(cells), "_", as.character(cells[1, ]))

  if (return.feature.meta) {
    feature.meta <- tglowr::get_feature_meta_from_names(cn)
  }

  # Read content of _cells.tsv
  cells <- data.table::fread(paste0(prefix, pat.cells), data.table = F, skip = 2, showProgress = FALSE)
  colnames(cells) <- cn

  # Read _image.tsv (metadata)
  img <- data.table::fread(paste0(prefix, pat.img), data.table = F, showProgress = FALSE)

  # Read _objectRelation.ships.tsv
  if (skip.orl) {
    orl <- NULL
  } else {
    orl <- data.table::fread(paste0(prefix, pat.orl), data.table = F, showProgress = FALSE)
  }

  # Standardize ID's across filesets into the following format
  # FS#I#O#
  # where FS = file set, I = image within fileset and O = object within image
  # This makes it easier to match across data with a unique ID
  if (add.global.id) {
    # Cell level information
    #-----------
    cells <- add_global_ids(cells)

    # IMG image level information
    #-----------
    img[, "ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])

    # ORL, object relationships
    #-----------
    if (!skip.orl) {
      if (nrow(orl) > 0) {
        orl[, "First Image Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"])
        orl[, "Second Image Number Global"] <- paste0(global.prefix, "_I", orl[, "Second Image Number"])
        orl[, "First Object Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
        orl[, "Second Object Number Global"] <- paste0(global.prefix, "_I", orl[, "Second Image Number"], "_O", orl[, "Second Object Number"])
      }
    }
  }

  out.list <- list(cells = cells, meta = img, orl = orl)

  if (return.feature.meta) {
    return(c(out.list, list(features = feature.meta)))
  } else {
    return(out.list)
  }
}

#-------------------------------------------------------------------------------
#' Read a cell level fileset type B
#'
#'
#' @description Reads a .zip file with cellprofiler features, each file other then
#' _Image, _Experiment and _Object Relationships are assumed to be an object
#' Will match objects on order with an appropriate matching strategy
#'
#' @param prefix Path prefix to fileset
#' @param return.feature.meta Should the dataframe with feature metadata be added
#' @param add.global.id Should extra id columns be added that are globally unique
#' @param merging.strategy How to consolidate 1:many relationships between cell: children. Accepted values: 'mean', 'none'
#' @param pat.exp The suffix pattern to identify the exeperiment file
#' @param pat.img The suffix pattern to identify the image level data
#' @param pat.cells The suffix pattern to identify the cell level data
#' @param pat.orl The suffix pattern to identify object relationships
#' @param pat.others The pattern to use to extract the child object names from the filename. First regex group is used as object name
#' @param na.rm Should NA's be removed when applying merging.strategy
#' @param skip.orl Should object relationships be read (not used, can be quite large)
#' @param verbose Should I be chatty?
#'
#' @returns list with data frames:
#' - cells (cell level features)
#' - meta (image level features)
#' - objectRelations
#' - features (optional)
#' Output is NULL if no cells are detected
#'
#' @export
read_cellprofiler_fileset_b <- function(prefix,
                                        return.feature.meta = F,
                                        add.global.id = T,
                                        merging.strategy = "mean",
                                        parent.col = "Parent_cell",
                                        pat.exp = ".*Experiment.txt",
                                        pat.img = ".*Image.txt",
                                        pat.cells = ".*cell.txt",
                                        pat.orl = ".*Object relationships.txt",
                                        pat.others = "^.*_([a-zA-Z]+\\d*).txt$",
                                        na.rm = F,
                                        skip.orl = F,
                                        verbose = F) {
  if (add.global.id) {
    assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", FILESET_ID)
  }

  index <- unzip(prefix, list = T)
  index$FileName <- basename(index$Name)

  tmpdir <- tempdir()

  # Clean up the tmpdir
  unlink(paste0(tmpdir, "/features"), recursive = T)

  # Unzip into tmp folder
  unzip(prefix, exdir = tmpdir)

  cells <- data.table::fread(paste0(tmpdir, "/", index[grep(pat.cells, index$FileName), "Name"]), data.table = F, showProgress = FALSE)
  img <- data.table::fread(paste0(tmpdir, "/", index[grep(pat.img, index$FileName), "Name"]), data.table = F, showProgress = FALSE)

  if (skip.orl) {
    orl <- NULL
  } else {
    orl <- data.table::fread(paste0(tmpdir, "/", index[grep(pat.orl, index$FileName), "Name"]), data.table = F, showProgress = FALSE)
  }

  if (nrow(cells) == 0) {
    warning("No cells detected for ", index[grep(pat.cells, index$FileName), "Name"], " returning NULL.")
    return(NULL)
  }

  exclude <- c(
    grep(pat.cells, index$FileName),
    grep(pat.img, index$FileName),
    grep(pat.orl, index$FileName),
    grep(pat.exp, index$FileName)
  )

  index <- index[!seq_len(nrow(index)) %in% exclude, ]

  colnames(cells) <- paste0("cell_", colnames(cells))
  rownames(cells) <- paste0(cells$cell_ImageNumber, "_", cells$cell_ObjectNumber)

  children <- list()
  if (nrow(index) > 0) {
    index$object <- gsub(pat.others, "\\1", index$FileName)

    for (i in seq_len(nrow(index))) {
      obj <- index[i, "object"]
      cur <- data.table::fread(paste0(tmpdir, "/", index[i, "Name"]), data.table = T, showProgress = FALSE)

      # Remove these columns from the merging strategy
      exclude.cols <- c("Group.1", grep("ObjectNumber", colnames(cur), value = T), grep("Number_Object_Number", colnames(cur), value = T))

      # If there are no cols, return NA
      if (nrow(cur) == 0) {
        # next()
        # If the file is empty, set these columns to NA
        cur <- as.data.frame(cur)
        colnames(cur) <- paste0(obj, "_", colnames(cur))
        cells[, c(colnames(cur), paste0(obj, "_QC_Object_Count"))] <- NA
        warning(paste0(obj, " assay for ", index[i, "Name"], " is empty. Returning NA for these cols."))
        
      } else if (merging.strategy == "mean") {
        if (verbose) cat("[DEBUG] ", as.character(index[i, ]), "\n")

        cur <- cur[as.logical(cur[[parent.col]] != 0), ]
        colnames(cur) <- paste0(obj, "_", colnames(cur))
        selector <- paste0(cur[[paste0(obj, "_ImageNumber")]], "_", cur[[paste0(obj, "_", parent.col)]])
        counts <- table(selector)
        cur$Group.1 <- selector

        if (verbose) cat("[DEBUG] NA's in selector", sum(is.na(selector)), "\n")
        if (verbose) cat("[DEBUG] Slice: ", head(cur$Group.1), "\n")

        # Calculate the mean per group
        tmp <- as.data.frame(cur[, lapply(.SD, mean, na.rm = na.rm), by = Group.1, .SDcols = colnames(cur)[!colnames(cur) %in% exclude.cols]])
        rownames(tmp) <- tmp$Group.1

        # Add the object count as a sanity check
        tmp[, paste0(obj, "_QC_Object_Count")] <- counts[tmp$Group.1]
        tmp <- tmp[, !colnames(tmp) %in% exclude.cols]

        # Assign the columns to the output matrix
        cells[selector, colnames(tmp)] <- tmp[selector, ]
        
      } else if (merging.strategy == "none") {
        if (add.global.id) {
          cur <- as.data.frame(cur)
          cur <- add_global_ids(cur)
        }
        children[[index[i, "object"]]] <- cur
      } else {
        stop("Only valid merging strategy is 'mean' or 'none'")
      }
    }
  }
  # Clean up the tmpdir
  unlink(paste0(tmpdir, "/features"), recursive = T)


  # Standardize ID's across filesets into the following format
  # FS#I#O#
  # where FS = file set, I = image within fileset and O = object within image
  # This makes it easier to match across data with a unique ID
  if (add.global.id) {
    # Cell level information
    #-----------
    cells <- add_global_ids(cells)

    # IMG image level information
    #-----------
    img[, "ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])
    # ORL, object relationships
    #-----------

    if (!skip.orl) {
      if (nrow(orl) > 0) {
        orl[, "First Image Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"])
        orl[, "First Image Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"])
        orl[, "First Object Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
        orl[, "First Object Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
      }
    }
  }

  if (return.feature.meta) {
    feature.meta <- tglowr::get_feature_meta_from_names(colnames(cells))
  }

  out.list <- list(cells = cells, meta = img, orl = orl)

  if (return.feature.meta) {
    out.list <- c(out.list, list(features = feature.meta))
  }

  if (length(children) > 0) {
    out.list <- c(out.list, list(children = children))
  }

  return(out.list)
}
