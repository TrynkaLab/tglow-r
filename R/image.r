

#-------------------------------------------------------------------------------
#' Convert a hex code to RGB
#' @param hex String with hex code
#'
#' @returns a vector of length 3 with relative scaling between RGB channels between 0 and 1
#' @export
hex_to_rgb <- function(hex) {
    if (startsWith(hex, "#")) {
        hex <- gsub("^#", "", hex)
    }

    rgb.hex <- strsplit(hex, "(?<=..)", perl = TRUE)[[1]]
    rgb <- sapply(rgb.hex, strtoi, base = 16)
    rgb <- rgb / 255

    names(rgb) <- c("r", "g", "b")
    return(rgb)
}

#-------------------------------------------------------------------------------
#' Apply a relative scaling of RGB values to a 3d array
#'
#' @description
#' 3rd dimension is assumed to contain the color channels.
#' If input is grey, returns the image in the color
#'
#' @param image 3d array with 3rd dimension as rgb
#' @param rgb vector of length 3 with scaling factors for rgb channels
#'
#' @returns The scaled image
#' @export
img_apply_color <- function(image, rgb) {
    check_package("EBImage")
    image[, , 1] <- image[, , 1] * rgb[1]
    image[, , 2] <- image[, , 2] * rgb[2]
    image[, , 3] <- image[, , 3] * rgb[3]

    return(image)
}

#-------------------------------------------------------------------------------
#' Combine a series of RGB images into a composite image
#'
#' @description
#' Apply a color shift and average their RGB values at each pixel to create a composite image
#'
#' @param images list of RGB EBImages
#' @param colors list of hex codes to color shift, must be same length as images
#'
#' @returns A single RGB images representing the average color values
#'
#' importFrom EBImage channel
#' @importFrom gplots col2hex
#' @export
img_composite <- function(images, colors) {
    check_package("EBImage")
    if (length(images) != length(colors)) {
        stop("Images and suplied colors don't match")
    }

    i <- 1
    comp <- NULL

    for (image in images) {
        if (length(dim(image)) != 2) {
            stop(paste0("Dimensions of image not equal to 2, is it greyscale. At index ", i))
        }

        rgb <- hex_to_rgb(gplots::col2hex(colors[i]))
        img <- EBImage::channel(image, mode = "rgb")

        if (is.null(comp)) {
            # Create an RGB image from greyscale
            comp <- img_apply_color(img, rgb)
        } else {
            comp <- comp + img_apply_color(img, rgb)
        }

        i <- i + 1
    }

    comp <- comp / length(images)

    return(comp)
}


#-------------------------------------------------------------------------------
#' Find the maximum value in each channel in a collection of images
#' This is used for normalizing images
#'
#' @param images A list of EBImage objects, with 3 or 4 dimensions. 3d dimension
#' is assumed to be channel
#' @param channel.dim Dimension in image stack that is channel
#' @param q The quantile to normalize to
#'
#' @returns A vector with the max values in each channel
#'
#' importFrom EBImage imageData
#' @export
img_max_per_channel <- function(images, channel.dim = 3, q = 1) {
    check_package("EBImage")
    channel.dim <- 3
    max <- c(rep(0, dim(images[[1]])[channel.dim]))

    for (img in images) {
        if (length(dim(img)) == 3) {
            for (ch in seq_len(dim(img)[channel.dim])) {
                if (quantile(img[, , ch], probs = q) > max[ch]) {
                    max[ch] <- quantile(img[, , ch], probs = q)
                }
            }
        } else if (length(dim(img)) == 4) {
            for (ch in seq_len(dim(img)[channel.dim])) {
                if (quantile(img[, , ch, ], probs = q) > max[ch]) {
                    max[ch] <- quantile(img[, , ch, ], probs = q)
                }
            }
        } else {
            stop(paste0(dim(img), " is not a valid number of dimensions for an image. "))
        }
    }

    return(max)
}

#-------------------------------------------------------------------------------
#' Normalize a list of EBImages to a common norm factor
#'
#' @param images A list of EBImage objects, with 3 or 4 dimensions. 3d dimension
#' is assumed to be channel.
#' @param norm.factors Value to divide the images by
#' Defaults to NULL in which case tglow.max.per.channel is called
#' @param q Quantile to normalize for, passed to \code{\link{img_max_per_channel}}
#'
#' importFrom EBImage imageData
#' @export
img_norm <- function(images, norm.factors = NULL, q = 1) {
    check_package("EBImage")
    if (is.null(norm.factors)) {
        cat("[INFO] Calculating norm factors\n")
        norm.factors <- img_max_per_channel(images, q = q)
    }

    if (length(dim(images[[1]])) == 3) {
        imgs <- lapply(seq_along(images), function(x) {
            clls <- images[[x]]

            for (ch in seq_len(dim(clls)[3])) {
                EBImage::imageData(clls)[, , ch] <- clls[, , ch] / norm.factors[ch]
            }
            return(clls)
        })
    } else if (length(dim(images[[1]])) == 4) {
        imgs <- lapply(seq_along(images), function(x) {
            clls <- images[[x]]

            for (ch in seq_len(dim(clls)[3])) {
                EBImage::imageData(clls)[, , ch, ] <- clls[, , ch, ] / norm.factors[ch]
            }
            return(clls)
        })
    } else {
        stop(paste0(dim(img), " is not a valid number of dimensions for an image. "))
    }


    return(imgs)
}

#-------------------------------------------------------------------------------
#' Max project a EBImage along the 4th axis
#'
#' @param img EBImage
#' @export
img_max_project <- function(img) {
    check_package("EBImage")
    apply(img, c(1, 2, 3), max)
}


#-------------------------------------------------------------------------------
#' Pad an array by zero-padding the first two dimensions with centering
#'
#' This function takes an input array or matrix of any number of dimensions and
#' pads or crops it along the first two dimensions to a specified target size.
#' The input is centered in the output array, and zeros are used for padding.
#' Dimensions beyond the first two are preserved.
#' If the input is larger than the target size, it is cropped accordingly.
#'
#' @param input An array or matrix of any dimensionality.
#' @param target.rows Integer. The desired number of rows in the output (first dimension).
#' @param target.cols Integer. The desired number of columns in the output (second dimension).
#'
#' @return An array with the first two dimensions padded or cropped to the target size,
#'   with the input array centered. Other dimensions remain the same as the input.
#'
#' @examples
#' input <- array(1, dim = c(5, 5, 3))
#' padded <- pad_center_nd(input, 10, 10)
#' dim(padded)  # 10 10 3
#' print(padded[,,1])  # Shows the first 2D slice padded and centered
#'
#' @export
img_pad_center <- function(input, target.rows, target.cols) {
  input_dims <- dim(input)
  other_dims <- if(length(input_dims) > 2) input_dims[3:length(input_dims)] else NULL
  
  if (is.null(other_dims)) {
    output_array <- array(0, dim = c(target.rows, target.cols))
  } else {
    output_array <- array(0, dim = c(target.rows, target.cols, other_dims))
  }
  
  start_row <- max(1, floor((target.rows - input_dims[1]) / 2) + 1)
  start_col <- max(1, floor((target.cols - input_dims[2]) / 2) + 1)
  end_row <- min(target.rows, start_row + input_dims[1] - 1)
  end_col <- min(target.cols, start_col + input_dims[2] - 1)
  
  input_start_row <- max(1, 1 + (1 - start_row))
  input_start_col <- max(1, 1 + (1 - start_col))
  input_end_row <- input_start_row + (end_row - start_row)
  input_end_col <- input_start_col + (end_col - start_col)
  
  idx_out <- list(start_row:end_row, start_col:end_col)
  idx_in <- list(input_start_row:input_end_row, input_start_col:input_end_col)
  if (!is.null(other_dims)) {
    for (i in seq_along(other_dims)) {
      idx_out[[2 + i]] <- 1:other_dims[i]
      idx_in[[2 + i]] <- 1:other_dims[i]
    }
  }
  
  output_array <- do.call(`[<-`, c(list(output_array), idx_out,
                                   list(value = do.call(`[`, c(list(input), idx_in)))))
  
  return(output_array)
}


#' Convert a matrix to a base64 text PNG
#' Takes a 2d (greyscale) or 3d (RGB) matrix and converts it to a base64 encoded png
#' 
#' @param mat A 2d or 3d matrix
#' @param add.html Should a HTML tag be added
#' @param width If add.html, what should the width be of the <img> tag
#' 
#' @returns a base64 encoded version of the image, or HTML <img> tag with the data
#' @export
img_to_base64png <- function(mat, add.html=F, width=250) {
  check_package("png")
  check_package("base64enc")
  
  raw_conn <- base::rawConnection(raw(0), "wb")
  png::writePNG(mat, raw_conn)
  
  img_data <- base::rawConnectionValue(raw_conn)
  base::close(raw_conn)
  
  # Convert the binary data to base64
  txt <- base64enc::base64encode(img_data)
  
  if (add.html) {
    txt <- paste0("<img src='data:image/png;base64,", txt, "' width='",width,"/>")
  }
  
  return(txt)
}

