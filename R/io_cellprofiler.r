#-------------------------------------------------------------------------------
# Imports
#' @import data.table
#' @include utils.r 
NULL


#-------------------------------------------------------------------------------
#' Read a cellprofiler fileset directory tree
#'
#' Matches images based on <plate>, <well>, <field>
#'
#' @param path path to tglow output dir
#' @param pattern The pattern that uniquely identfies a well. Use .zip for type B and the _experiment.tsv
#' for type A.
#' @param type Must be A or B.
#' Type A: _cells.tsv, _image.tsv, _experiment.tsv and _objectRelations.tsv.
#' See read.cellprofiler.fileset.a for detaills
#' Type B: <plate>_<well>.zip with individual files for each child object. Main object is assumed to be _cells.
#' see
#' See read.cellprofiler.fileset.b for detaills
#' @param n Read a subset of filesets. If integer, only that fileset is read, otherwise specify indices to read.
#' @param verbose Should I be chatty?
#' @param ... Remaining parameters passed to read.cellprofiler.fileset.a/b
#'
#' @export
read.cellprofiler.dir <- function(path, pattern, type, n = NULL, verbose = F, ...) {
  files <- list.files(path, recursive = T, pattern = paste0("*", pattern), full.names = T)

  if (type == "A") {
    prefixes <- gsub(pattern, "", files)
  } else {
    prefixes <- files
  }

  if (!is.null(n)) {
    prefixes <- prefixes[n]
  }

  # Reset global fileset index
  assign("FILESET_ID", 0, envir = .GlobalEnv)

  # Read filesets
  filesets <- list()
  i <- 0
  for (pre in prefixes) {
    i <- i + 1
    cat("\r[INFO] reading fileset ", i, "/", length(prefixes))

    if (type == "A") {
      cur <- read.cellprofiler.fileset.a(pre, return.feature.meta = F, ...)
    } else if (type == "B") {
      cur <- read.cellprofiler.fileset.b(pre, return.feature.meta = F, ...)
    } else {
      stop(paste0("Invalid type: ", type))
    }

    if (!is.null(cur)) {
      filesets[[pre]] <- cur

      if (verbose) {
        cat("\n[DEBUG] cols:", ncol(cur$cells), " cols meta", ncol(cur$meta), " cols orl:", ncol(cur$orl), "\n")
      }
    } else {
      warning("Fileset was NULL, skipped.")
    }
  }

  cat("\n[INFO] Merging filesets\n")
  output <- merge.filesets(filesets)

  cat("[INFO] names: ", names(output), "\n")

  if (verbose) {
    cat("[DEBUG] colnames:\n", colnames(output$cells))
  }

  features <- get.feature.meta.from.names(colnames(output$cells))
  classes <- sapply(output$cells, class)
  features$type <- classes[features$id]

  features <- features[colnames(output$cells), ]
  rownames(output$cells) <- output$cells$cells_ObjectNumber_Global
  rownames(output$meta) <- output$meta$ImageNumber_Global

  return(c(output, list(features = features)))
}

#-------------------------------------------------------------------------------
#' Add a global id to a matrix of image files.
#' Matrix with column pattern 'ImageNumber', 'ObjectNumber' and 'Object_Number'
#' which will be duplicated and have a globally unique variable added.
#' Output of this is stored in the same column name but with suffix _Global
#'
#' @export
add.global.ids <- function(matrix) {
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
#' Reads a CellProfiler fileset into a list.
#' Type A: Assumes all features are in a single _cells.tsv / _cells.tsv
#' and all features are matched
#'
#' @returns list with data frames:
#'
#' - cells (cell level features)
#'
#' - meta (image level features)
#'
#' - objectRelations
#'
#' - features [optional]
#'
#'
#' Output is NULL if no cells are detected.
#'
#' @export
read.cellprofiler.fileset.a <- function(prefix,
                                        return.feature.meta = F,
                                        add.global.id = T,
                                        pat.img = "_image.tsv",
                                        pat.cells = "_cells.tsv",
                                        pat.orl = "_objectRelationships.tsv") {
  if (add.global.id) {
    assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", FILESET_ID)
  }

  # Read header of _cells.tsv
  cells <- fread(paste0(prefix, pat.cells), data.table = F, nrows = 3)

  # If the file has no cells empty
  if (nrow(cells) != 3) {
    return(NULL)
  }

  # Clean colnames
  cn <- paste0(colnames(cells), "_", as.character(cells[1, ]))

  if (return.feature.meta) {
    feature.meta <- tglow.get.feature.meta.from.cells(cn)
  }

  # Read content of _cells.tsv
  cells <- fread(paste0(prefix, pat.cells), data.table = F, skip = 2)
  colnames(cells) <- cn

  # Read _image.tsv (metadata)
  img <- fread(paste0(prefix, pat.img), data.table = F)

  # Read _objectRelation.ships.tsv
  orl <- fread(paste0(prefix, pat.orl), data.table = F)


  # Standardize ID's across filesets into the following format
  # FS#I#O#
  # where FS = file set, I = image within fileset and O = object within image
  # This makes it easier to match across data with a unique ID
  if (add.global.id) {
    # Cell level information
    #-----------
    cells <- add.global.ids(cells)

    # IMG image level information
    #-----------
    img[, "ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])

    # ORL, object relationships
    #-----------
    if (nrow(orl) > 0) {
      orl[, "First Image Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"])
      orl[, "Second Image Number Global"] <- paste0(global.prefix, "_I", orl[, "Second Image Number"])
      orl[, "First Object Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
      orl[, "Second Object Number Global"] <- paste0(global.prefix, "_I", orl[, "Second Image Number"], "_O", orl[, "Second Object Number"])
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
#' Reads a zip file with cellprofiler features, each file other then
#' _Image, _Experiment and _Object Relationships are assumed to be an object.
#' Will match objects on order with an appropriate matching strategy.
#'
#' @returns list with data frames:
#' - cells (cell level features)
#' - meta (image level features)
#' - objectRelations
#' - features [optional]
#' Output is NULL if no cells are detected.
#'
#' @export
read.cellprofiler.fileset.b <- function(prefix,
                                        return.feature.meta = F,
                                        add.global.id = T,
                                        merging.strategy = "mean",
                                        parent.col = "Parent_cell",
                                        pat.exp = ".*Experiment.txt",
                                        pat.img = ".*Image.txt",
                                        pat.cells = ".*cell.txt",
                                        pat.orl = ".*Object relationships.txt",
                                        pat.others = "^.*_([a-z]+\\d*).txt$",
                                        na.rm = F) {
  if (add.global.id) {
    assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", FILESET_ID)
  }

  index <- unzip(prefix, list = T)
  index$FileName <- basename(index$Name)


  # Clean up the tmpdir
  unlink(paste0(tmpdir, "/features"), recursive = T)

  # Unzip into tmp folder
  tmpdir <- tempdir()
  unzip(prefix, exdir = tmpdir)

  cells <- fread(paste0(tmpdir, "/", index[grep(pat.cells, index$FileName), "Name"]), data.table = F)
  img <- fread(paste0(tmpdir, "/", index[grep(pat.img, index$FileName), "Name"]), data.table = F)
  orl <- fread(paste0(tmpdir, "/", index[grep(pat.orl, index$FileName), "Name"]), data.table = F)

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

  index <- index[!1:nrow(index) %in% exclude, ]

  colnames(cells) <- paste0("cell_", colnames(cells))
  rownames(cells) <- paste0(cells$cell_ImageNumber, "_", cells$cell_ObjectNumber)

  children <- list()
  if (nrow(index) > 0) {
    index$object <- gsub(pat.others, "\\1", index$FileName)

    for (i in 1:nrow(index)) {
      obj <- index[i, "object"]
      cur <- fread(paste0(tmpdir, "/", index[i, "Name"]), data.table = T)

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
        cur <- cur[as.logical(cur[[parent.col]] != 0), ]
        colnames(cur) <- paste0(obj, "_", colnames(cur))
        selector <- paste0(cur[[paste0(obj, "_ImageNumber")]], "_", cur[[paste0(obj, "_", parent.col)]])
        counts <- table(selector)
        # exclude.cols      <- c("Group.1", grep("ObjectNumber", colnames(cur), value=T), grep("Number_Object_Number", colnames(cur), value=T))

        cur$Group.1 <- selector

        tmp <- as.data.frame(cur[, lapply(.SD, mean, na.rm = na.rm), by = Group.1, .SDcols = colnames(cur)[!colnames(cur) %in% exclude.cols]])
        rownames(tmp) <- tmp$Group.1

        # The above is the data.table equivalent of aggregate, which is MUCH faster
        # tmp               <- aggregate(cur, by=list(selector), FUN=mean, na.rm=T)

        # Add the object count as a sanity check
        tmp[, paste0(obj, "_QC_Object_Count")] <- counts[tmp$Group.1]
        tmp <- tmp[, !colnames(tmp) %in% exclude.cols]

        # Assign the columns to the output matrix
        cells[selector, colnames(tmp)] <- tmp[selector, ]
      } else if (merging.strategy == "none") {
        if (add.global.id) {
          cur <- add.global.ids(cur)
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
    cells <- add.global.ids(cells)

    # IMG image level information
    #-----------
    img[, "ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])
    # ORL, object relationships
    #-----------
    if (nrow(orl) > 0) {
      orl[, "First Image Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"])
      orl[, "First Image Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"])
      orl[, "First Object Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
      orl[, "First Object Number Global"] <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
    }
  }

  if (return.feature.meta) {
    feature.meta <- tglow.get.feature.meta.from.cells(colnames(cells))
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
