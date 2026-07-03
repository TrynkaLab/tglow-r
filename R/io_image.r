#-------------------------------------------------------------------------------
#' Read a binary matrix, used for reading the registration matrices
#'
#' @description
#' Adapted from: https://www.r-bloggers.com/2012/06/getting-numpy-data-into-r/
#'
#' @param path Path to the binary matrix file
#'
#' @export
tglow_read_binmat <- function(path) {
    con <- file(path, "rb")
    on.exit(close(con))
    dims <- readBin(con, "integer", 2)
    raw <- readBin(con, "numeric", prod(dims))
    if (length(raw) != prod(dims)) {
        stop("File truncated or corrupt: expected ", prod(dims), " values, got ", length(raw))
    }
    mat <- matrix(raw, dims[1], dims[2])
    return(mat)
}

#-------------------------------------------------------------------------------
#' Read images for objects in a TglowDataset from the tglow-pipeline .h5 format
#' 
#' @description
#' Read image data from a .h5 file that is orgnaized as such: 
#' /plate/row/col/field.h5 where each h5 group represents a cells image
#' 
#' @param dataset A \linkS4class{TglowDataset}
#' @param objects A selection vector for which objects to read images for
#' @param index.col Column on dataset@meta which has the h5 groups (object ids)
#' @param out.size Pad the images to a constant square size setting 0. NULL = no padding.
#' @param path The root path for the .h5 files. 'root path'/plate/row/col/field.h5
#' @param path.col Column in dataset@meta with the paths to .h5 file. Overrides path if not NULL
#' 
#' @returns A list with EBImages and channel names
#' @export
tglow_read_imgs <- function(dataset, objects, index.col, out.size=NULL, path=NULL, path.col=NULL) {
  
  check_package("EBImage")
  check_package("hdf5r")

  if (is.null(path) && is.null(path.col)) {
    stop("Need to set either the root 'path' or specify a @meta column with the path through 'path.col'")
  }
  
  if (is.null(dataset@feature.map)) {
    stop("@feature.map must be set on object when calling read_images")
  }
  
  if (!is.null(path.col)) {
    if (!isAvailable(dataset, c(path.col), assay=NULL, slot=NULL)) {
      stop(paste0("path.col: ", path.col, " not available on dataset"))
    }
  }
  
  if (!index.col %in% colnames(dataset@meta)) {
    stop(paste0("index.col: ", index.col, " not available on @meta of this dataset"))
  }
  
  
  object.names <- dataset@object.ids[objects]
  dataset      <- dataset[objects,]
  object.index <- getDataByObject(dataset, index.col)
  
  # Check if object index have NAs in there (the h5 group names)
  if (any(is.na(object.index))) {
    warning("NA's found in index.col. Removing these")
    dataset      <- dataset[!is.na(object.index),]
    object.names <- object.names[!is.na(object.index)]
    object.index <- object.index[!is.na(object.index)]
  }
  
  if (!is.character(object.index)) {
    stop("h5 group index (provided through index.col) must be a chracter")
  }
  
  # Construct the path to the h5 file
  plate <- getDataByObject(dataset, dataset@feature.map@plate@feature, assay=dataset@feature.map@plate@assay, slot=dataset@feature.map@plate@slot)
  well  <- getDataByObject(dataset, dataset@feature.map@well@feature, assay=dataset@feature.map@well@assay, slot=dataset@feature.map@well@slot)
  field <- getDataByObject(dataset, dataset@feature.map@field@feature, assay=dataset@feature.map@field@assay, slot=dataset@feature.map@field@slot)
  
  pwf   <- data.frame(plate=plate, well=well, field=field, cell_index=object.index, object_id=object.names)
  pwf   <- cbind(pwf, t(sapply(pwf[,2], well_to_index)))
  rownames(pwf) <- object.names
  
  if (is.null(path.col)) {
    # col is unpadded numeric (see well_to_index()); directory layout uses plain numbers, not zero-padded
    pwf$h5_path <- paste0(path, "/", pwf[,1], "/", pwf[,"row_letter"], "/", pwf[,"col"], "/", pwf[,"field"], ".h5")
  } else {
    pwf$h5_path <- getDataByObject(dataset, path.col)
  }
  
  # Aggregate the unique files, this avoids unneded IO
  files <- unique(pwf$h5_path)
  
  images        <- list()
  channel.names <- NULL
  
  pb <- progress::progress_bar$new(format = "[INFO] Reading images [:bar] :current/:total (:percent) eta :eta", total = length(object.names))
  for (h5file in files) {
    
    cur.pwf <- pwf[pwf$h5_path == h5file,]
    file    <- hdf5r::H5File$new(h5file, mode = "r")

    tryCatch({
      for (obj in 1:nrow(cur.pwf)) {
        cur.pwf.obj <- cur.pwf[obj,,drop=F]
        m           <- file[[cur.pwf.obj$cell_index]]$read()
        if (!is.null(out.size)) {
          m <- img_pad_center(m, out.size, out.size)
        }
        images[[cur.pwf.obj$object_id]] <- EBImage::Image(m)
        pb$tick()
      }

      # Try to read channel names group if available
      if (is.null(channel.names)) {
        if ("channel_names" %in% names(file)) {
          channel.names <- file[["channel_names"]]$read()
        }
      } else if ("channel_names" %in% names(file)) {
        cur.channel.names <- file[["channel_names"]]$read()
        if (!identical(cur.channel.names, channel.names)) {
          warning("channel_names in ", h5file, " differ from previously read channel_names. Using the first file's channel_names for all objects.")
        }
      }
    }, finally = {
      file$close_all()
    })
  }

  if (!all(object.names %in% names(images))) {
    stop("Missing images for some objects — internal inconsistency in tglow_read_imgs")
  }

  return(list(img=images[object.names], channels=channel.names, channel.order="XY[Z]C"))
  
}

