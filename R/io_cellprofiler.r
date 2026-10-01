#-------------------------------------------------------------------------------
#' Read a cellprofiler fileset directory tree
#'
#' @description
#' Read a cellprofiler fileset directory tree organized into a <plate>/<well>/<field>.fileset
#' structure
#'
#' @param path path to tglow output dir
#' @param pattern The pattern that uniquely identifies a fileset. Use '.zip' for type 'B'
#' and the '_Image.txt' or '_Experiment.txt'(or however you exported the image data) for type 'A'
#' @param type Must be 'A' or 'B'. See details
#' @param n Read a subset of filesets. If integer, only that fileset is read, otherwise specify indices to read
#' @param skip.orl Skip reading of object relationships, as this can get quite large with many children and is not used
#' @param verbose Should I be chatty?
#' @param col.object The column name in the features which contains the per object object identifier. See details
#' @param col.meta.img.id The column name in the image level data which contains the image id. See details
#' @param max.per.well Maximumn number of objects in a well to consider. Set to NULL to ignore
#' @param ... Remaining parameters passed to \code{\link{read_cellprofiler_fileset_a}} or \code{\link{read_cellprofiler_fileset_b}}
#'
#' @details
#'
#' `type`
#' Type A: _cell.txt, _Image.txt, _Experiment.txt and _Object relationships.txt
#' See \code{\link{read_cellprofiler_fileset_a}} for details
#'
#' Type B: <plate>_<well>.zip with individual files for each child object. Main object is assumed to be _cells
#' See \code{\link{read_cellprofiler_fileset_b}} for details
#'
#'
#' `col.object` and `col.meta.img.id`
#'
#' See \code{\link{add_global_ids}} for details on how globally unique Id's are assigned in the case you need to set
#' `col.object` or `col.meta.img.id``.
#' @importFrom progress progress_bar
#' @export
read_cellprofiler_dir <- function(path, pattern, type, n = NULL, skip.orl = TRUE, verbose = F, col.object = "cell_ObjectNumber_Global", col.meta.img.id = "ImageNumber_Global", max.per.well=NULL, ...) {
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
  # assign("FILESET_ID", 0, envir = .GlobalEnv)
  fileset.id <- 1
  
  # Read filesets
  filesets <- list()
  null.filesets <- 0
  pb <- progress::progress_bar$new(format = "[INFO] Reading [:bar] :current/:total (:percent) eta :eta, mem :mem (GB)", total = length(prefixes))
  pb$tick(0)
  for (pre in prefixes) {
    if (type == "A") {
      cur <- tglowr::read_cellprofiler_fileset_a(pre, return.feature.meta = F, skip.orl = skip.orl, fileset.id = fileset.id, ...)
    } else if (type == "B") {
      cur <- tglowr::read_cellprofiler_fileset_b(pre, return.feature.meta = F, skip.orl = skip.orl, verbose = verbose, fileset.id = fileset.id, ...)
    } else {
      stop(paste0("Invalid type: ", type))
    }
    
    fileset.id <- fileset.id + 1
    
    if (!is.null(cur)) {
      if (!is.null(max.per.well)) {
        if (nrow(cur[["cells"]]) > max.per.well) {
          cat(paste0("\nWell ", pre, " has ", nrow(cur[["cells"]]) ," (> ", max.per.well, ") objects, skipped. Increase max.per.well to ignore\n"))
          
          warning(paste0("Well ", pre, " has ", nrow(cur[["cells"]]) ," (> ", max.per.well, ") objects, skipped. Increase max.per.well to ignore"))
          next
        }
      }
      
      filesets[[pre]] <- cur
      if (verbose) cat("\n[DEBUG] cols:", ncol(cur$cells), " cols meta", ncol(cur$meta), " cols orl:", ncol(cur$orl), "\n")
    } else {
      null.filesets <- null.filesets + 1
      # warning("Fileset was NULL, skipped.")
    }
    
    memory_used <- format(object.size(filesets), units="GB", digits=2)
    pb$tick(tokens = list(mem = memory_used))      
  }
  
  if (null.filesets != 0) {
    msg <- paste0("Dectected ", null.filesets, "/", length(filesets), " as NULL (no cells)")
    warning(msg)
    cat(paste0("[WARN] ", msg, "\n"))
  }
  
  cat("[INFO] Read filesets into list of ", format(object.size(filesets), units="GB", digits=2), "GB\n")
  
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
    warning(paste0("Detected ", sum(!selector), " objects with NA in ", col.object, " removing these"))
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
#' @param fileset.id The global fileset id to add
#'
#' @details
#' Matrix with column pattern 'ImageNumber', 'ObjectNumber' and 'Object_Number'
#' which will be duplicated and have a globally unique variable added
#' Output of this is stored in the same column name but with suffix _Global
#'
#' @returns A data.frame with extra columns suffixed by _Global with globally unique ids
#' @export
add_global_ids <- function(matrix, fileset.id) {
  # Fetch global fileset id
  global.prefix <- paste0("FS", fileset.id)

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
#' Type A: Assumes all features are in a single _cell.txt (configurable with  pat.cells)
#' and all features are matched
#'
#' @param prefix Path prefix to fileset
#' @param return.feature.meta Should the dataframe with feature metadata be added
#' @param add.global.id Should extra id columns be added that are globally unique
#' @param pat.img The suffix pattern to identify the image level data
#' @param pat.cells The suffix pattern to identify the cell level data
#' @param pat.orl The suffix pattern to identify object relationships
#' @param skip.orl Should object relationships be read (not used, can be quite large)
#' @param fileset.id The global fileset id to add if add.global.id = TRUE
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
#' @importFrom data.table fread
#' @export
read_cellprofiler_fileset_a <- function(prefix,
                                        return.feature.meta = F,
                                        add.global.id = T,
                                        pat.img = "_Image.txt",
                                        pat.cells = "_cell.txt",
                                        pat.orl = "_Object relationships.txt",
                                        skip.orl = FALSE,
                                        fileset.id = NULL) {
  if (add.global.id) {
    if (is.null(fileset.id)) {
      stop("fileset.id cannot be NULL if add.global.id=T")
    }
    # assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", fileset.id)
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
    cells <- add_global_ids(cells, fileset.id)

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
#' @param pat.exp The suffix pattern to identify the experiment file
#' @param pat.img The suffix pattern to identify the image level data
#' @param pat.cells The suffix pattern to identify the cell level data
#' @param pat.orl The suffix pattern to identify object relationships
#' @param pat.others The pattern to use to extract the child object names from the filename. First regex group is used as object name
#' @param na.rm Should NA's be removed when applying merging.strategy
#' @param skip.orl Should object relationships be read (not used, can be quite large)
#' @param skip.children List of child object names to skip during reading.
#' @param verbose Should I be chatty?
#' @param fileset.id The global fileset id to add if add.global.id = TRUE
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
                                        skip.children = NULL,
                                        verbose = F,
                                        fileset.id = NULL) {
  if (add.global.id) {
    if (is.null(fileset.id)) {
      stop("fileset.id cannot be NULL if add.global.id=T")
    }
    # assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", fileset.id)
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
      
      # Optionally skip adding this child object
      if (obj %in% skip.children) {
        next()
      }
      
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

        if (!parent.col %in% colnames(cur)) {
          stop(paste0("Parent col: '", parent.col, "' not found for child object '", obj, "'. Either update the files, or skip reading the child object with skip.children"))
        }

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
    cells <- add_global_ids(cells, fileset.id)

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


#-------------------------------------------------------------------------------
# Parquet readers for the gamma pipeline
#-------------------------------------------------------------------------------

# Feature maps per tglow-pipeline version (manifest.version in nextflow.config)
.PIPELINE_FEATURE_MAPS <- list(
  "0.2.0" = list(x = "centroid_x", y = "centroid_y", z = "centroid_z", plate = "plate", well = "well", field = "field")
)

# Fixed metadata columns in measure_intensity_features_with_debris.py output (+ the ids added on reading)
.PIPELINE_OBJECT_META <- c("plate", "row", "col", "well", "field", "plate_well_field", "cell_label",
                           "centroid_x", "centroid_y", "centroid_z", "plate_id", "image_id", "object_id")
.PIPELINE_IMAGE_META <- c("plate", "row", "col", "well", "field", "plate_well_field", "method",
                          "cellmask_expansion", "n_cells_reg_corr_pass", "plate_id", "image_id")
.PIPELINE_META_PATTERNS <- c("__registration_corr$")

# Per channel debris statistics, used to categorize features
.PIPELINE_DEBRIS_STATS <- c("threshold", "mean_intensity", "threshold_mean_ratio", "debris_percentage",
                            "ori_cell_mask_covered_percent", "cell_mask_covered_percent")


#-------------------------------------------------------------------------------
#' Build a TglowFeatureMap from a named list of meta level feature names
#' @noRd
.feature_map_from_list <- function(features) {
  feature.map <- TglowFeatureMap()
  for (cur.slot in names(features)) {
    slot(feature.map, cur.slot) <- TglowFeatureLocation(features[[cur.slot]])
  }
  return(feature.map)
}


#-------------------------------------------------------------------------------
#' Default feature map for CellProfiler parquet output
#'
#' @description
#' Returns the \linkS4class{TglowFeatureMap} matching the output of \code{\link{read_cellprofiler_parquet}}.
#' x/y/z are taken from `@meta`, plate, well and field from `@image.meta`.
#'
#' @returns A \linkS4class{TglowFeatureMap}
#' @export
tglow_feature_map_cellprofiler <- function() {
  return(.feature_map_from_list(list(
    x = "cell_Location_Center_X",
    y = "cell_Location_Center_Y",
    z = "cell_Location_Center_Z",
    plate = "Metadata_plate",
    well = "Metadata_well",
    field = "Metadata_field"
  )))
}


#-------------------------------------------------------------------------------
#' Feature map for tglow-pipeline intensity output
#'
#' @description
#' Returns the \linkS4class{TglowFeatureMap} matching the output of \code{\link{read_pipeline_parquet}}
#' for a given tglow-pipeline version (manifest.version in the pipeline's nextflow.config).
#'
#' @param version Pipeline version, or "latest" for the most recent supported version
#'
#' @returns A \linkS4class{TglowFeatureMap}
#' @export
tglow_feature_map_pipeline <- function(version = "latest") {
  versions <- names(.PIPELINE_FEATURE_MAPS)

  if (version == "latest") {
    version <- as.character(max(numeric_version(versions)))
  }

  if (!version %in% versions) {
    stop(paste0("Pipeline version '", version, "' is not supported. Supported versions: ", paste(versions, collapse = ", ")))
  }

  return(.feature_map_from_list(.PIPELINE_FEATURE_MAPS[[version]]))
}


#-------------------------------------------------------------------------------
#' Get feature metadata from tglow-pipeline feature names
#'
#' @description
#' Counterpart of \code{\link{get_feature_meta_from_names}} for tglow-pipeline names of the
#' form `ch{N}__{stat}`. Names without '__' (area, n_nuclei) are assigned to object 'cell'.
#'
#' @param feature.names Character vector of feature names
#'
#' @returns A data.frame with columns id, object, measurement, category and name
#' @export
get_feature_meta_from_names_pipeline <- function(feature.names) {
  pos <- regexpr("__", feature.names, fixed = TRUE)
  has.channel <- pos > 0

  object <- ifelse(has.channel, substr(feature.names, 1, pos - 1), "cell")
  measurement <- ifelse(has.channel, substr(feature.names, pos + 2, nchar(feature.names)), feature.names)

  category <- rep("intensity", length(feature.names))
  category[grepl("^background_", measurement)] <- "background"
  category[measurement %in% .PIPELINE_DEBRIS_STATS] <- "debris"
  category[measurement == "registration_corr"] <- "registration"
  category[!has.channel] <- "morphology"

  feature.meta <- data.frame(
    id = feature.names,
    object = object,
    measurement = measurement,
    category = category,
    name = measurement
  )
  rownames(feature.meta) <- feature.meta$id
  return(feature.meta)
}


#-------------------------------------------------------------------------------
#' Read and row bind a set of parquet files
#'
#' @param files Character vector of parquet files
#' @param label Description of the files, used in messages
#' @param idcol If not NULL, add a column with this name holding the source file of each row
#'
#' @returns A data.frame
#' @noRd
.read_parquet_files <- function(files, label = "files", idcol = NULL) {
  # Skip empty files (stub outputs from the pipeline)
  empty <- is.na(file.size(files)) | file.size(files) == 0
  if (any(empty)) {
    warning(paste0("Skipped ", sum(empty), " empty ", label, ": ", paste(files[empty], collapse = ", ")))
    files <- files[!empty]
  }

  if (length(files) == 0) {
    stop(paste0("No non-empty ", label, " to read"))
  }

  pb <- progress::progress_bar$new(format = paste0("[INFO] Reading ", label, " [:bar] :current/:total (:percent) eta :eta"), total = length(files))
  pb$tick(0)
  tables <- lapply(files, function(file) {
    cur <- as.data.frame(arrow::read_parquet(file))
    pb$tick()
    return(cur)
  })
  names(tables) <- files

  # Report columns not present in every file, these are filled with NA
  col.counts <- table(unlist(lapply(tables, colnames)))
  missing <- col.counts[col.counts != length(tables)]
  if (length(missing) > 0) {
    msg <- paste0("Not all ", label, " have the same columns, filling with NA. Missing in n files: ",
                  paste0(names(missing), " (", length(tables) - missing, ")", collapse = ", "))
    warning(msg)
  }

  return(as.data.frame(data.table::rbindlist(tables, use.names = TRUE, fill = TRUE, idcol = idcol)))
}


#-------------------------------------------------------------------------------
#' Split a data.frame into metadata and a numeric feature matrix
#'
#' @param df Input data.frame
#' @param meta.patterns Regex patterns, matching columns are considered metadata
#' @param meta.cols Column names considered metadata
#'
#' @details Non numeric columns are always considered metadata.
#'
#' @returns list with meta (data.frame), features (numeric matrix) and types (original class per feature)
#' @noRd
.split_meta_features <- function(df, meta.patterns = NULL, meta.cols = NULL) {
  # integer64 does not behave as numeric, convert to double
  for (col in colnames(df)[sapply(df, inherits, "integer64")]) {
    df[[col]] <- as.numeric(df[[col]])
  }

  is.meta <- !sapply(df, is.numeric) | colnames(df) %in% meta.cols
  for (pattern in meta.patterns) {
    is.meta <- is.meta | grepl(pattern, colnames(df))
  }

  features <- as.matrix(df[, !is.meta, drop = F])
  storage.mode(features) <- "double"

  return(list(
    meta = df[, is.meta, drop = F],
    features = features,
    types = sapply(df[, !is.meta, drop = F], function(x) class(x)[1])
  ))
}


#-------------------------------------------------------------------------------
#' Remove columns matching any of a set of perl regex patterns
#'
#' @param df Input data.frame
#' @param patterns Perl regex patterns, NULL drops nothing
#' @param keep Columns that are never dropped
#' @param label Description of the columns, used in messages
#' @param verbose Should I be chatty?
#' @noRd
.drop_columns <- function(df, patterns, keep = NULL, label = "", verbose = FALSE) {
  drop <- rep(FALSE, ncol(df))
  for (pattern in patterns) {
    drop <- drop | grepl(pattern, colnames(df), perl = TRUE)
  }
  drop <- drop & !colnames(df) %in% keep

  if (verbose) {
    cat("[DEBUG] dropping ", sum(drop), " ", label, " columns: ", paste(colnames(df)[drop], collapse = ", "), "\n", sep = "")
  }

  return(df[, !drop, drop = F])
}


#-------------------------------------------------------------------------------
#' Warn for feature map features not present on a dataset
#' @noRd
.check_feature_map <- function(dataset) {
  for (cur.slot in slotNames(dataset@feature.map)) {
    loc <- slot(dataset@feature.map, cur.slot)

    # Slot not set
    if (length(loc@feature) == 0) {
      next
    }

    if (is.null(loc@assay)) {
      available <- loc@feature %in% c(colnames(dataset@meta), colnames(dataset@image.meta))
    } else {
      available <- loc@feature %in% colnames(slot(dataset@assays[[loc@assay]], loc@slot))
    }

    if (!available) {
      warning(paste0("Feature map ", cur.slot, " feature '", loc@feature, "' not found on dataset"))
    }
  }
}


#-------------------------------------------------------------------------------
#' Set the plate column from the plate folder each row was read from
#'
#' @param df Output of .read_parquet_files with idcol ".source_file"
#' @param file.plate Named vector, names are files, values the plate folder
#' @param label Description of the rows, used in messages
#' @noRd
.set_folder_plate <- function(df, file.plate, label) {
  folder.plate <- unname(file.plate[df$.source_file])

  mismatch <- df$plate != folder.plate
  if (any(mismatch, na.rm = T)) {
    warning(paste0("plate column does not match the plate folder for ", sum(mismatch, na.rm = T), " ", label, ", using the folder name"))
  }

  df$plate <- folder.plate
  df$.source_file <- NULL
  return(df)
}


#-------------------------------------------------------------------------------
#' Warn for non numeric columns that are not in the expected metadata columns
#' @noRd
.warn_unexpected_meta <- function(df, expected, label) {
  unexpected <- setdiff(colnames(df)[!sapply(df, is.numeric)], expected)
  if (length(unexpected) > 0) {
    warning(paste0("Unexpected non numeric columns in ", label, ", placing on meta: ", paste(unexpected, collapse = ", ")))
  }
}


#-------------------------------------------------------------------------------
#' Build a TglowDataset from split object and image level data
#'
#' @param obj Output of .split_meta_features for the objects, with rownames set
#' @param img Output of .split_meta_features for the images, with rownames set
#' @param image.ids Image id for each object
#' @param feature.meta.fun Function to generate the feature metadata from feature names
#' @param feature.map TglowFeatureMap or NULL
#' @param assay.out The assay name to store objects under
#' @noRd
.build_parquet_dataset <- function(obj, img, image.ids, feature.meta.fun, feature.map, assay.out = "raw") {
  # No numeric image features, use a dummy as TglowDatasetFromList does
  if (ncol(img$features) == 0) {
    img$features <- matrix(0, nrow = nrow(img$meta), ncol = 1, dimnames = list(rownames(img$meta), "dummy"))
    img$types <- c(dummy = "numeric")
  }

  overlap <- intersect(colnames(obj$meta), colnames(img$meta))
  if (length(overlap) > 0) {
    warning(paste0("Columns present in both @meta and @image.meta: ", paste(overlap, collapse = ", ")))
  }

  dataset <- tglowr::TglowDatasetFromMatrices(obj$features, img$features, image.ids,
                                              object.meta = obj$meta, image.meta = img$meta, assay.out = assay.out)

  # Feature level metadata
  features <- feature.meta.fun(colnames(obj$features))
  features$type <- obj$types[features$id]
  features$analyze <- TRUE
  dataset@assays[[assay.out]]@features <- features

  features <- feature.meta.fun(colnames(img$features))
  features$type <- img$types[features$id]
  features$analyze <- TRUE
  dataset@image.data@features <- features

  dataset@active.assay <- assay.out
  dataset@feature.map <- feature.map

  if (!is.null(feature.map)) {
    .check_feature_map(dataset)
  }

  cat("[INFO] Read ", nrow(obj$features), " objects with ", ncol(obj$features), " features and ",
      nrow(img$features), " images with ", ncol(img$features), " features\n", sep = "")

  return(dataset)
}


#-------------------------------------------------------------------------------
#' Read per-plate CellProfiler parquet files
#'
#' @description
#' Read the per-plate `<plate>_cells.parquet` and `<plate>_image.parquet` files produced by
#' concat_cellprofiler in the gamma version of tglow-pipeline into a \linkS4class{TglowDataset}.
#'
#' @param path Directory to search recursively for parquet files, or a character vector of parquet files
#' @param pattern.cells Pattern identifying the object level files. Removing it from the filename gives the plate name
#' @param pattern.image Pattern identifying the image level files. Removing it from the filename gives the plate name
#' @param plates Character vector of plates to read. NULL reads all plates found
#' @param meta.patterns Regex patterns for object level columns that are put on `@meta` instead of the assay
#' @param img.meta.patterns Regex patterns for image level columns that are put on `@image.meta` instead of `@image.data`
#' @param drop.patterns Perl regex patterns for object level columns that are removed entirely. NULL keeps all columns
#' @param img.drop.patterns Perl regex patterns for image level columns that are removed entirely. NULL keeps all columns
#' @param col.object Column with the globally unique object id
#' @param col.img.id Column in the object level data with the globally unique image id
#' @param col.meta.img.id Column in the image level data with the globally unique image id
#' @param feature.map \linkS4class{TglowFeatureMap} to set on the dataset, or NULL. See \code{\link{tglow_feature_map_cellprofiler}}
#' @param assay.out The assay name to store objects under
#' @param verbose Should I be chatty?
#'
#' @details
#' Child objects are expected to be merged onto the parent objects and object/image ids to be globally
#' unique already, as concat_cellprofiler does.
#'
#' Columns matching `drop.patterns` / `img.drop.patterns` are removed first. By default these are:
#' - plate and well, which duplicate Metadata_plate and Metadata_well on `@image.meta`
#' - `<child>_Parent_*` (and `_Global`) for all non cell objects. Parent_cell equals the cell's own object number after
#'   merging, other child parent columns are averaged over the children during merging and no longer refer to an object
#' - `<child>_ImageNumber` and `<child>_ObjectNumber`/`Number_Object_Number` (and `_Global`) for all non cell objects
#'
#' `col.object` and `col.img.id` are never dropped.
#'
#' Columns that are not numeric are always placed on `@meta` (object level) or `@image.meta` (image level).
#' Numeric columns matching `meta.patterns` / `img.meta.patterns` are placed there as well, the
#' remaining numeric columns form the assay and `@image.data`.
#'
#' @returns A \linkS4class{TglowDataset}
#' @importFrom arrow read_parquet
#' @export
read_cellprofiler_parquet <- function(path,
                                      pattern.cells = "_cells.parquet$",
                                      pattern.image = "_image.parquet$",
                                      plates = NULL,
                                      meta.patterns = c("ImageNumber", "ObjectNumber", "Object_Number", "Parent",
                                                        "_Location_", "BoundingBox", "^global_"),
                                      img.meta.patterns = c("ImageNumber", "^Metadata_", "^Group_", "^ExecutionTime_",
                                                            "^ModuleError_", "^Frame_", "^Series_", "^Height_", "^Width_",
                                                            "^global_"),
                                      drop.patterns = c("^plate$", "^well$", "^(?!cell_)[^_]+_Parent_",
                                                        "^(?!cell_)[^_]+_(ImageNumber|ObjectNumber|Number_Object_Number)(_Global)?$"),
                                      img.drop.patterns = c("^plate$", "^well$"),
                                      col.object = "cell_ObjectNumber_Global",
                                      col.img.id = "cell_ImageNumber_Global",
                                      col.meta.img.id = "ImageNumber_Global",
                                      feature.map = tglow_feature_map_cellprofiler(),
                                      assay.out = "raw",
                                      verbose = FALSE) {
  if (length(path) == 1 && dir.exists(path)) {
    files <- list.files(path, recursive = T, full.names = T)
  } else {
    files <- path
  }

  # Pair cell and image files by plate
  cell.files <- grep(pattern.cells, files, value = T)
  img.files <- grep(pattern.image, files, value = T)
  names(cell.files) <- sub(pattern.cells, "", basename(cell.files))
  names(img.files) <- sub(pattern.image, "", basename(img.files))

  unpaired <- setdiff(union(names(cell.files), names(img.files)), intersect(names(cell.files), names(img.files)))
  if (length(unpaired) > 0) {
    stop(paste0("Plates without both a cells and image file: ", paste(unpaired, collapse = ", ")))
  }

  if (is.null(plates)) {
    plates <- sort(names(cell.files))
  }

  if (length(plates) == 0) {
    stop(paste0("No files matching '", pattern.cells, "' found"))
  }

  missing <- setdiff(plates, names(cell.files))
  if (length(missing) > 0) {
    stop(paste0("Plates not found: ", paste(missing, collapse = ", ")))
  }

  cells <- .read_parquet_files(cell.files[plates], label = "cell files")
  img <- .read_parquet_files(img.files[plates], label = "image files")

  # Check ids
  for (col in c(col.object, col.img.id)) {
    if (!col %in% colnames(cells)) {
      stop(paste0("Column ", col, " not found in cell files"))
    }
  }

  if (!col.meta.img.id %in% colnames(img)) {
    stop(paste0("Column ", col.meta.img.id, " not found in image files"))
  }

  selector <- !is.na(cells[[col.object]])
  if (sum(!selector) != 0) {
    warning(paste0("Detected ", sum(!selector), " objects with NA in ", col.object, " removing these"))
    cells <- cells[selector, , drop = F]
  }

  if (any(duplicated(cells[[col.object]]))) {
    stop(paste0("Duplicated ids in ", col.object, ". Check plates were not given the same plate_id in the pipeline"))
  }

  if (any(duplicated(img[[col.meta.img.id]]))) {
    stop(paste0("Duplicated ids in ", col.meta.img.id, ". Check plates were not given the same plate_id in the pipeline"))
  }

  missing <- setdiff(cells[[col.img.id]], img[[col.meta.img.id]])
  if (length(missing) > 0) {
    stop(paste0(length(missing), " image ids in ", col.img.id, " not found in the image files"))
  }

  # Remove redundant columns
  cells <- .drop_columns(cells, drop.patterns, keep = c(col.object, col.img.id), label = "object", verbose = verbose)
  img <- .drop_columns(img, img.drop.patterns, keep = col.meta.img.id, label = "image", verbose = verbose)

  # Split into meta and features
  obj <- .split_meta_features(cells, meta.patterns = meta.patterns)
  rownames(obj$meta) <- cells[[col.object]]
  rownames(obj$features) <- cells[[col.object]]

  im <- .split_meta_features(img, meta.patterns = img.meta.patterns)
  rownames(im$meta) <- img[[col.meta.img.id]]
  rownames(im$features) <- img[[col.meta.img.id]]

  if (verbose) {
    cat("[DEBUG] object meta cols: ", ncol(obj$meta), " image meta cols: ", ncol(im$meta), "\n")
  }

  return(.build_parquet_dataset(obj, im,
    image.ids = cells[[col.img.id]],
    feature.meta.fun = tglowr::get_feature_meta_from_names,
    feature.map = feature.map,
    assay.out = assay.out
  ))
}


#-------------------------------------------------------------------------------
#' Read tglow-pipeline intensity parquet files
#'
#' @description
#' Read the per-well object_features.parquet and image_features.parquet files produced by
#' measure_intensity_features_with_debris in the gamma version of tglow-pipeline into a \linkS4class{TglowDataset}.
#' Files are expected as `<path>/<plate>/<row>/<col>/object_features.parquet` and image_features.parquet.
#'
#' @param path The output directory containing one folder per plate
#' @param plates Character vector of plate folder names to read. NULL reads all plate folders in alphabetical order
#' @param feature.map \linkS4class{TglowFeatureMap} to set on the dataset, or NULL. See \code{\link{tglow_feature_map_pipeline}}
#'
#' @details
#' Object ids are constructed as `<plate_id>_<well>_I<field>_L<cell_label>`, image ids as `<plate_id>_<well>_I<field>`.
#' The L indicates the cell mask label, which is not the same as the CellProfiler ObjectNumber.
#'
#' Plate ids (P1, P2, ...) follow the order of `plates`. They match the plate ids in the CellProfiler parquet
#' output only if `plates` is given in the same order as the pipeline manifest, which is how the pipeline
#' assigns them. To align with a CellProfiler dataset, either provide `plates` in manifest order, or use
#' \code{\link{match_objects_xy_nn}}, which matches objects on plate name, well, field and position and
#' does not depend on the ids.
#'
#' Metadata columns are fixed. Registration correlations are placed on `@meta` and n_cells_reg_corr_pass on
#' `@image.meta`. Image level metadata (plate, row, col, well, field, plate_well_field, plate_id) is placed on
#' `@image.meta` only, `@meta` holds the object level metadata (cell_label, centroids, registration correlations)
#' and image_id. Image level columns are still available per object through \code{\link{getDataByObject}}.
#'
#' @returns A \linkS4class{TglowDataset}
#' @importFrom arrow read_parquet
#' @export
read_pipeline_parquet <- function(path, plates = NULL, feature.map = tglow_feature_map_pipeline()) {
  if (!dir.exists(path)) {
    stop(paste0("Directory not found: ", path))
  }

  if (is.null(plates)) {
    plates <- sort(list.dirs(path, recursive = F, full.names = F))
  }

  if (any(duplicated(plates))) {
    stop("plates contains duplicates")
  }

  missing <- plates[!dir.exists(file.path(path, plates))]
  if (length(missing) > 0) {
    stop(paste0("Plate folders not found in ", path, ": ", paste(missing, collapse = ", ")))
  }

  # Find files per plate
  obj.files <- c()
  img.files <- c()
  for (plate in plates) {
    cur.obj <- Sys.glob(file.path(path, plate, "*", "*", "object_features.parquet"))
    cur.img <- Sys.glob(file.path(path, plate, "*", "*", "image_features.parquet"))

    if (length(cur.obj) == 0 || length(cur.img) == 0) {
      warning(paste0("No object_features.parquet or image_features.parquet found for plate ", plate, ", skipping"))
      next
    }

    obj.files <- c(obj.files, setNames(cur.obj, rep(plate, length(cur.obj))))
    img.files <- c(img.files, setNames(cur.img, rep(plate, length(cur.img))))
  }

  if (length(obj.files) == 0) {
    stop(paste0("No parquet files found in ", path))
  }

  # Plate ids follow the order of plates
  plates <- plates[plates %in% names(obj.files)]
  plate.ids <- setNames(paste0("P", seq_along(plates)), plates)
  cat("[INFO] Plate ids:\n")
  cat(paste0("  ", plate.ids, " = ", plates, "\n"), sep = "")

  # Read, the source file is used to set the plate from the folder name
  objects <- .read_parquet_files(obj.files, label = "object files", idcol = ".source_file")
  objects <- .set_folder_plate(objects, setNames(names(obj.files), obj.files), "objects")
  images <- .read_parquet_files(img.files, label = "image files", idcol = ".source_file")
  images <- .set_folder_plate(images, setNames(names(img.files), img.files), "images")

  # Older 2D output has no centroid_z, use 0 as the pipeline does for 2D input
  if (!"centroid_z" %in% colnames(objects)) {
    objects$centroid_z <- 0
  } else {
    objects$centroid_z[is.na(objects$centroid_z)] <- 0
  }

  # Construct ids
  objects$plate_id <- unname(plate.ids[objects$plate])
  images$plate_id <- unname(plate.ids[images$plate])
  images$image_id <- paste0(images$plate_id, "_", images$well, "_I", images$field)
  objects$image_id <- paste0(objects$plate_id, "_", objects$well, "_I", objects$field)
  objects$object_id <- paste0(objects$image_id, "_L", objects$cell_label)

  if (any(duplicated(objects$object_id))) {
    stop("Duplicated object ids, check for duplicate wells or fields in the input")
  }

  if (any(duplicated(images$image_id))) {
    stop("Duplicated image ids, check for duplicate wells or fields in the input")
  }

  missing <- setdiff(objects$image_id, images$image_id)
  if (length(missing) > 0) {
    stop(paste0(length(missing), " image ids of objects not found in the image files"))
  }

  # Non numeric columns that are not part of the expected output still go to meta
  .warn_unexpected_meta(objects, .PIPELINE_OBJECT_META, "objects")
  .warn_unexpected_meta(images, .PIPELINE_IMAGE_META, "images")

  # Split into meta and features
  obj <- .split_meta_features(objects, meta.patterns = .PIPELINE_META_PATTERNS, meta.cols = .PIPELINE_OBJECT_META)
  rownames(obj$meta) <- objects$object_id
  rownames(obj$features) <- objects$object_id

  im <- .split_meta_features(images, meta.patterns = .PIPELINE_META_PATTERNS, meta.cols = .PIPELINE_IMAGE_META)
  rownames(im$meta) <- images$image_id
  rownames(im$features) <- images$image_id

  # Image level metadata (plate, well, field, ...) is kept on @image.meta only, image_id on @meta only
  obj$meta <- obj$meta[, !colnames(obj$meta) %in% colnames(im$meta) | colnames(obj$meta) == "image_id", drop = F]
  im$meta <- im$meta[, colnames(im$meta) != "image_id", drop = F]

  return(.build_parquet_dataset(obj, im,
    image.ids = objects$image_id,
    feature.meta.fun = get_feature_meta_from_names_pipeline,
    feature.map = feature.map
  ))
}
