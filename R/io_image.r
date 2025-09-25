#-------------------------------------------------------------------------------
#' Build an index of where the source images are stored
#'
#' @description
#' This assumes images are organized as follows:
#' <plate>/<row>/<col>/<field>.ome.tiff
#'
#' Matches images based on <plate>/<row>/<col>/<field>.ome.tiff
#'
#' @param path path to tglow image dir
#' @param plate_filter plate names to run
#' @param pattern the pattern to match image names and then extract the fields (gsubbed away)
#'
#' @returns A data frame with the image paths
#' @export
depr_tglow_build_img_index <- function(path, plate_filter = NULL, pattern = ".ome.tiff$") {
    warning("[WARN] Deprecation warning. This will not be actively maintained. Please use the new way of pre-caching imagecrops.")

    plates <- list.files(path, recursive = F)

    if (!is.null(plate_filter)) {
        plates <- plates[plates %in% plate_filter]
    }

    # cat("[INFO] Indexing plates: ", plates, "\n")

    index <- data.frame(matrix(NA, nrow = 0, ncol = 6))

    nfiles <- 0
    for (plate in plates) {
        cat("[INFO] Indexing plate: ", plate, "\n")

        rows <- list.files(paste0(path, "/", plate), pattern = "^[A-Z]$")
        for (row in rows) {
            cols <- list.files(paste0(path, "/", plate, "/", row), pattern = "^\\d+$")
            for (col in cols) {
                files <- list.files(paste0(path, "/", plate, "/", row, "/", col), pattern = pattern)
                nfiles <- nfiles + length(files)
                for (file in files) {
                    field <- gsub(pattern, "", file)
                    tmp <- c(plate, paste0(row, sprintf("%02d", as.numeric(col))), row, col, field, paste0(path, "/", plate, "/", row, "/", col, "/", file))
                    index <- rbind(index, tmp)
                }
            }
        }
    }
    colnames(index) <- c("plate", "well", "row", "col", "field", "path")

    cat("[INFO] Indexed ", nfiles, " image files\n")
    if (nfiles == 0) {
        warning("No files detected, suggest you check pattern is correct.")
    } else {
        rownames(index) <- paste0(index$plate, ":", index$well, ":", index$field)
    }


    return(index)
}

#-------------------------------------------------------------------------------
#' Read a binary matrix, used for reading the registration matrices
#'
#' @description
#' Adapted from: https://www.r-bloggers.com/2012/06/getting-numpy-data-into-r/
#'
#' @param path to binary matrix
#'
#' @export
tglow_read_binmat <- function(path) {
    con <- file(path, "rb")
    dim <- readBin(con, "integer", 2)
    mat <- matrix(readBin(con, "numeric", prod(dim)), dim[1], dim[2])
    close(con)
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
#' @param path The root path for the .h5 files. 'root path'/plate/row/field.h5
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
  
  if (!isAvailable(dataset, c(index.col), assay=NULL, slot=NULL)) {
    stop(paste0("index.col: ", index.col, " not available on dataset"))
  }
  
  
  object.names <- dataset@object.ids[objects]
  dataset <- dataset[objects,]
  
  plate <- getDataByObject(dataset, dataset@feature.map@plate@feature, assay=dataset@feature.map@plate@assay, slot=dataset@feature.map@plate@slot)
  well  <- getDataByObject(dataset, dataset@feature.map@well@feature, assay=dataset@feature.map@well@assay, slot=dataset@feature.map@well@slot)
  field <- getDataByObject(dataset, dataset@feature.map@field@feature, assay=dataset@feature.map@field@assay, slot=dataset@feature.map@field@slot)
  
  pwf <- data.frame(plate=plate, well=well, field=field, cell_index=getDataByObject(dataset, index.col), object_id=object.names)
  pwf <- cbind(pwf, t(sapply(pwf[,2], well_to_index)))
  rownames(pwf) <- object.names
  
  if (is.null(path.col)) {
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
      } else {
        channel.names <- NULL
      }
    }
    
    file$close_all()
  }

  return(list(img=images[object.names], channels=channel.names, channel.order="XY[Z]C"))
  
}


#-------------------------------------------------------------------------------
#' Load images around the center of a cell object
#'
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param cell.subset The subset of cells to retrieve
#' @param group.col The feature to match cells to img.index. Usually 'Metadata_group'
#' @param img.index A matrix with image paths generated by tglow.build.img.index
#' @param reg.index A matrix with registration paths generated by tglow.build.img.index
#' Must have an additional collumn named 'plate2' to indicate the registraiton plate
#' @param window Window in px around the cell center to retrieve
#' @param feature.x The feature that describes the object pos in px in x
#' @param feature.y The feature that describes the object pos in px in y
#' @param channels The channels to read. A list of Vector of indices per cycle. (default NULL = all channels)
#' @param planes The planes to read. A list of Vector of indices per cycle. (default NULL = all planes)
#' @param max.project Should the stack be max projected per channel? (default TRUE)
#' @param assay The assay to grab feature.x and feature.y from. Not used for anything else. Defaults to 'raw'. 
#' If features are on meta slot, leave this default.
#' @param slot The slot to use for feature.x and feature.y, defaults to "data". Can be "data" or "scale.data"
#' 
#' importFrom EBImage colorMode affine combine Grayscale colorMode<-
#' importFrom RBioFormats read.image
#' @returns A list of EBImage objects by object id
#' @export
depr_tglow_read_imgs <- function(dataset,
                            assay = "raw",
                            slot = "data",
                            cell.subset,
                            group.col,
                            img.index,
                            reg.index = NULL,
                            window = 75,
                            feature.x = "cell_AreaShape_Center_X",
                            feature.y = "cell_AreaShape_Center_Y",
                            channels = NULL,
                            planes = NULL,
                            max.project = T) {
    check_dataset_assay_slot(dataset, assay, slot)
    check_package("EBImage")
    check_package("RBioFormats")
    warning("[WARN] Deprecation warning. This will not be actively maintained. Please use the new way of pre-caching imagecrops.")
    
    #cur.cells <- slot(dataset@assays[[assay]], slot)
       
    # Change this to getDataByObject instead
    cur.cells <- getDataByObject(dataset, j=c(feature.x, feature.y), assay=assay, slot=slot, drop=F)
     
    if ((class(cell.subset) != "character") && (class(cell.subset) != "integer")) {
        stop("Cell.subset must be character or integer")
    }

    if (class(cell.subset) == "integer") {
        warning("cell.subset is integer, make sure the indices are properly matched to assay.")
        if (max(cell.subset) > nrow(cur.cells)) {
            stop("Index in cell.subset larger then number of rows in assay. Are you using the correct assay?")
        }
    }

    pb <- progress::progress_bar$new(format = "[INFO] Reading images [:bar] :current/:total (:percent) eta :eta", total = length(cell.subset))
    out <- lapply(cell.subset, function(i) {
        pb$tick()
        cur.group <- dataset@image.meta[dataset@image.ids[i], group.col]
        cur.img <- img.index[cur.group, ]
        x.pos <- cur.cells[i, feature.x]
        y.pos <- cur.cells[i, feature.y]

        # cat("[INFO] Reading ", cur.img, "\n")
        cx <- round((x.pos - window):(x.pos + window))
        cy <- round((y.pos - window):(y.pos + window))

        subset <- list(X = cx, Y = cy)

        if (!is.null(channels)) {
            if (class(channels) != "list") {
                stop("Channels must be a list with one numeric vector per cycle containing channel indices to load.")
            }
            subset[["C"]] <- channels[[1]]
        }

        if (!is.null(planes)) {
            if (class(planes) != "list") {
                stop("Planes must be a list with one numeric vector per cycle containing plane indices to load.")
            }
            subset[["Z"]] <- planes[[1]]
        }

        img <- RBioFormats::read.image(cur.img$path, normalize = F, subset = subset)
        EBImage::colorMode(img) <- EBImage::Grayscale

        # Append mutliple cycles
        if (!is.null(reg.index)) {
            cur.reg.index <- reg.index[reg.index$plate == cur.img$plate & reg.index$well == cur.img$well & reg.index$field == cur.img$field, ]

            for (p2 in 1:nrow(cur.reg.index)) {
                new.img <- img.index[paste0(
                    cur.reg.index[p2, "plate2"], ":",
                    cur.reg.index[p2, "well"], ":",
                    cur.reg.index[p2, "field"]
                ), ]

                # Possible fix, use base::aperm() and a flexible order argument, to reorder the matrices after reading
                # for now this is too much effort. 
               # order <- c()
                if (!is.null(channels)) {
                    stop("Currently subset + reg.index is broken, please subset manually after reading and set channels=NULL when registering")
                    if (length(channels) != nrow(cur.reg.index) + 1) {
                        stop("Length of channel list does not match the number of cycles found in reg.index")
                    }
                    subset[["C"]] <- as.numeric(channels[[p2 + 1]])
                }

                if (!is.null(planes)) {
                    stop("Currently subset + reg.index is broken, please subset manually after reading and set planes=NULL when registering")
                    if (length(planes) != nrow(cur.reg.index) + 1) {
                        stop("Length of plane list does not match the number of cycles found in reg.index")
                    }
                    subset[["Z"]] <- planes[[p2 + 1]]
                }

                # Read the registration matrix
                #cat("[DEBUG] reading cycle x: ", cur.reg.index[p2, "path"], "\n")
                reg <- tglowr::tglow_read_binmat(cur.reg.index[p2, "path"])
                #cat("[DEBUG] read registration: \n", reg, "\n")

                #cat("[DEBUG] reading image: ", new.img$path, "\n")
                #cat("[DEBUG] susbet image: \n")
                #print(subset)
                tmp.img <- RBioFormats::read.image(new.img$path, subset = subset, normalize = F)
                #tmp.img <- base::aperm(tmp.img, order)
                EBImage::colorMode(tmp.img) <- EBImage::Grayscale
                #cat("[DEBUG] read image: ", dim(tmp.img), "\n")

                # Recode the matrix so it works with EBImage
                m <- t(reg)[, 1:2]
                m[3, ] <- -1 * m[3, ]
                #cat("[DEBUG] prepared registration mat \n", m, "\n")

                # Apply affine transform
                tmp.img <- EBImage::affine(tmp.img, m, filter = "none", antialias = F)
                #cat("[DEBUG] Applied affine to cycle x\n")


                #cat("[DEBUG] ", dim(img), " - ", dim(tmp.img), "\n")
                # Combine with the image
                #img <-c(img, tmp.img)
                # YXCZ
                img <- EBImage::abind(img, tmp.img, along = 3)
                #cat("[DEBUG] ", dim(img), "\n")

            }
        }

        #cat("[DEBUG] Done reading aditional cycles\n")
        #cat("[DEBUG] ", dim(img), "\n")

        crop <- img

        if (length(dim(img)) == 4) {
            # Max project
            if (max.project) {
                crop <- apply(img, c(1, 2, 3), max)
            }
        } else if (length(dim(img)) == 3) {
            # Max project
            if (max.project) {
                crop <- apply(img, c(1, 2), max)
            }
        } else if (length(dim(img)) == 2) {
            crop <- img
        }
        return(crop)
    })

    cat("\n")
    names(out) <- dataset@object.ids[cell.subset]

    return(out)
}



#-------------------------------------------------------------------------------
#' Load images around the center of a cell object using AICSImageIO
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param cell.subset The subset of cells to retrieve
#' @param python python install with AICSImageIO
#' @param group.col The feature to match cells to img.index. Usually 'Metadata_group'
#' @param img.index A matrix with image paths generated by tglow.build.img.index
#' @param reg.index A matrix with registration paths generated by tglow.build.img.index
#' Must have an additional collumn named 'plate2' to indicate the registraiton plate
#' @param window Window in px around the cell center to retrieve
#' @param feature.x The feature that describes the object pos in px in x
#' @param feature.y The feature that describes the object pos in px in y
#' @param channels The channels to read. A list of Vector of indices per cycle. (default NULL = all channels)
#' @param planes The planes to read. A list of Vector of indices per cycle. (default NULL = all planes)
#' @param max.project Should the stack be max projected per channel? (default TRUE)
#' @param assay The assay to grab feature.x and feature.y from. Not used for anything else. Defaults to 'raw'. 
#' If features are on meta slot, leave this default.
#' @param slot The slot to use for feature.x and feature.y, defaults to "data". Can be "data" or "scale.data"
#' 
#' importFrom EBImage colorMode affine combine Grayscale colorMode<-
#' importFrom RBioFormats read.image
#' @description 
#' 
#' Same as tglow_read_imgs but uses AICSImageIO as the backend. This is slower then tglow_read_imgs.
#' 
#' @returns A list of EBImage objects by object id
#' @export
depr_tglow_read_imgs_aics <- function(dataset,
                            cell.subset,
                            python,
                            group.col,
                            img.index,
                            reg.index = NULL,
                            window = 75,
                            feature.x = "cell_AreaShape_Center_X",
                            feature.y = "cell_AreaShape_Center_Y",
                            channels = NULL,
                            planes = NULL,
                            max.project = T,
                            assay = "raw",
                            slot = "data") {
    check_dataset_assay_slot(dataset, assay, slot)
    check_package("reticulate")
    warning("[WARN] Deprecation warning. This will not be actively maintained. Please use the new way of pre-caching imagecrops.")

    reticulate::use_python(python)
    aics <- reticulate::import("aicsimageio")

    # Fetch x and y positions
    cur.cells <- getDataByObject(dataset, j=c(feature.x, feature.y), assay=assay, slot=slot, drop=F)
     
    if ((class(cell.subset) != "character") && (class(cell.subset) != "integer")) {
        stop("Cell.subset must be character or integer")
    }

    if (class(cell.subset) == "integer") {
        warning("cell.subset is integer, make sure the indices are properly matched to assay.")
        if (max(cell.subset) > nrow(cur.cells)) {
            stop("Index in cell.subset larger then number of rows in assay. Are you using the correct assay?")
        }
    }

    pb <- progress::progress_bar$new(format = "[INFO] Reading images [:bar] :current/:total (:percent) eta :eta", total = length(cell.subset))
    out <- lapply(cell.subset, function(i) {
        pb$tick()
        cur.group <- dataset@image.meta[dataset@image.ids[i], group.col]
        cur.img <- img.index[cur.group, ]
        x.pos <- cur.cells[i, feature.x]
        y.pos <- cur.cells[i, feature.y]

        # cat("[INFO] Reading ", cur.img, "\n")
        cx <- round((x.pos - window):(x.pos + window))
        cy <- round((y.pos - window):(y.pos + window))

        subset <- list(X = cx - 1, Y = cy - 1)

        aics.img <- aics$AICSImage(cur.img$path)

        if (!is.null(channels)) {
            if (class(channels) != "list") {
                stop("Channels must be a list with one numeric vector per cycle containing channel indices to load.")
            }
            subset[["C"]] <- array(channels[[1]] - 1)
        } else {
            subset[["C"]] <- 0:(aics.img$dims[["C"]][[1]] - 1)
        }

        if (!is.null(planes)) {
            if (class(planes) != "list") {
                stop("Planes must be a list with one numeric vector per cycle containing plane indices to load.")
            }
            subset[["Z"]] <- array(planes[[1]] - 1)
        } else {
            subset[["Z"]] <- 0:(aics.img$dims[["Z"]][[1]] - 1)
        }

        d <- aics.img$get_image_dask_data("XYCZ", T=0, X=subset[["X"]], Y=subset[["Y"]],  C=subset[["C"]], Z=subset[["Z"]])
        img <- EBImage::Image(d$compute())
        #tmp.img <- aics.img$get_image_data("XYCZ", T=0, Y=subset[["Y"]], X=subset[["X"]], C=subset[["C"]], Z=subset[["Z"]])

        #img <- RBioFormats::read.image(cur.img$path, normalize = F, subset = subset)
        EBImage::colorMode(img) <- EBImage::Grayscale

        # Append mutliple cycles
        if (!is.null(reg.index)) {
            cur.reg.index <- reg.index[reg.index$plate == cur.img$plate & reg.index$well == cur.img$well & reg.index$field == cur.img$field, ]

            for (p2 in 1:nrow(cur.reg.index)) {
                new.img <- img.index[paste0(
                    cur.reg.index[p2, "plate2"], ":",
                    cur.reg.index[p2, "well"], ":",
                    cur.reg.index[p2, "field"]
                ), ]

                
                aics.img2 <- aics$AICSImage(new.img$path)

                if (!is.null(channels)) {
                    if (length(channels) != nrow(cur.reg.index) + 1) {
                        stop("Length of channel list does not match the number of cycles found in reg.index")
                    }
                    subset[["C"]] <- array(channels[[p2 + 1]] - 1)
                } else {
                    subset[["C"]] <- 0:(aics.img2$dims[["C"]][[1]] - 1)
                }

                if (!is.null(planes)) {
                    if (length(planes) != nrow(cur.reg.index) + 1) {
                        stop("Length of plane list does not match the number of cycles found in reg.index")
                    }
                    subset[["Z"]] <- array(planes[[p2 + 1]] - 1)
                } else {
                    subset[["Z"]] <- 0:(aics.img2$dims[["Z"]][[1]] - 1)
                }

                # Read the registration matrix
                reg <- tglowr::tglow_read_binmat(cur.reg.index[p2, "path"])
                
                # Read the image
                d2 <- aics.img2$get_image_dask_data("XYCZ", T=0, X=subset[["X"]], Y=subset[["Y"]], C=subset[["C"]], Z=subset[["Z"]])
                tmp.img <- EBImage::Image(d2$compute())
                #tmp.img <- aics.img2$get_image_data("XYCZ", T=0, Y=subset[["Y"]], X=subset[["X"]], C=subset[["C"]], Z=subset[["Z"]])
    
                EBImage::colorMode(tmp.img) <- EBImage::Grayscale

                # Recode the matrix so it works with EBImage
                m <- t(reg)[, 1:2]
                m[3, ] <- -1 * m[3, ]

                # Apply affine transform
                tmp.img <- EBImage::affine(tmp.img, m, filter = "none", antialias = F)

                # YXCZ
                img <- EBImage::abind(img, tmp.img, along = 3)
            }
        }

        #cat("[DEBUG] Done reading aditional cycles\n")
        #cat("[DEBUG] ", dim(img), "\n")

        crop <- img

        if (length(dim(img)) == 4) {
            # Max project
            if (max.project) {
                crop <- apply(img, c(1, 2, 3), max)
            }
        } else if (length(dim(img)) == 3) {
            # Max project
            if (max.project) {
                crop <- apply(img, c(1, 2), max)
            }
        } else if (length(dim(img)) == 2) {
            crop <- img
        }
        return(crop)
    })

    cat("\n")
    names(out) <- dataset@object.ids[cell.subset]

    return(out)
}


