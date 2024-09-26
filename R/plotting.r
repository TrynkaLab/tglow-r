#------------------------------------------------------------------------------
#' Simple plotting theme wrapper for ggplot objects
#'
#' @param p A ggplot object
#' @param base_size Text size
#' @param base_family Text font
#' @param legend Keep the legend
#'
#' @importFrom ggplot2 theme theme_grey %+replace%
#' @export
theme_plain <- function(p, base_size = 11, base_family = "ArialMT", legend = TRUE) {
    p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", size = 0.75),
            axis.ticks = element_line(size = 0.75),
            axis.text = element_text(size = base_size, family = base_family, face = "plain"),
            strip.background = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size = base_size, family = base_family, face = "plain"),
            complete = TRUE,
            plot.title = element_text(hjust = 0.5)
        )

    if (!legend) {
        p <- p + theme(legend.position = "none")
    }
    return(p)
}


#-------------------------------------------------------------------------------
#' Chart of average execution times
#'
#' @description
#' Plot the execution time of CellProfiler modules
#'
#' @param object A TglowDataset
#' @param object As percentage, logical. Should values be returned as a percentage

#' @returns ggplot object containing the plot
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme element_text element_blank element_line ylab xlab
#' @export
tglow_plot_execution_time <- function(object, as.percentage = FALSE) {
    if (!is(object, "TglowDataset")) {
        stop("Object must be a TglowDataset")
    }
    meta <- object@image.meta

    # Calculate average execution time
    df.plot <- meta[, grep("ExecutionTime", colnames(meta))]

    if (ncol(df.plot) <= 0) {
        stop("No columns of pattern 'ExecutionTime' found on object@image.meta")
    }

    df.plot <- data.frame(name = gsub("ExecutionTime_", "", colnames(df.plot)), time = colMeans(df.plot, na.rm = T))

    if (as.percentage) {
        df.plot$time <- (df.plot$time / sum(df.plot$time)) * 100
    }

    df.plot <- df.plot[order(df.plot$time, decreasing = T), ]
    df.plot$name <- factor(df.plot$name, levels = unique(df.plot$name))

    plot <- ggplot(df.plot, aes(x = name, y = time, fill = time)) +
        geom_bar(width = 0.8, stat = "identity") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey"),
            panel.background = element_blank()
        ) +
        ylab("Average time (seconds)") +
        xlab("")

    return(plot)
}


#-------------------------------------------------------------------------------
#' Plot a reduction
#'
#' @description Plot a reduction on a \linkS4class{TglowDataset}
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param reduction A name of a reduction on dataset
#' @param ident The item to use for coloring points, passed to \code{\link{getDataByObject}}
#' @param assay The assay to use for coloring, passed to \code{\link{getDataByObject}}
#' @param slot The slot to use for coloring, passed to \code{\link{getDataByObject}} Can be "data" or "scale.data"
#' @param downsample Should downsampling be applied prior to plot. NA's are removed first
#' @param log.ident Should color vector be log2 transformed
#' @param axis.x Which dimension from reduction to plot in x
#' @param axis.y Which dimension from reduction to plot in y
#' @param xlab Override xlab
#' @param ylab Override ylab
#' @param no.colscale Skip adding colorscale
#' @param ... Remaining arguments passed to \code{\link{plot_xy}}
#'
#' @export
tglow_dimplot <- function(object, reduction, ident = NULL, assay = NULL, slot = NULL, downsample = NULL, log.ident = FALSE, axis.x = 1, axis.y = 2, xlab = NULL, ylab = NULL, no.colscale = FALSE, ...) {
    if (!reduction %in% names(object@reduction)) {
        stop("Reduction not found on object")
    }

    if (!is.numeric(downsample) && !is.integer(downsample) && !is.logical(downsample)) {
        stop("Downsample must be a number indicating samples to select, or a numeric or logical selection vector")
    }

    dim <- object@reduction[[reduction]]@x[, c(axis.x, axis.y)]

    if (is.null(ident)) {
        col <- "blue"
    } else {
        col <- getDataByObject(object, ident, assay = assay, slot = slot)
        col <- col[rowSums(is.na(dim)) != ncol(dim)]
    }

    dim <- dim[rowSums(is.na(dim)) != ncol(dim), ]

    if (!is.null(downsample)) {
        if (length(downsample) == 1) {
            if (downsample > nrow(dim)) {
                stop("Downsample larger then number of non NA rows")
            }

            downsample <- sample.int(nrow(dim), downsample)
        }

        dim <- dim[downsample, , drop = F]
        col <- col[downsample]
    }


    if (is.null(xlab)) {
        xlab <- paste0(reduction, " - ", axis.x)
    }

    if (is.null(ylab)) {
        ylab <- paste0(reduction, " - ", axis.y)
    }

    if (log.ident) {
        col <- log2(log.ident)
    }

    p <- theme_plain(plot_xy(
        dim[, axis.x],
        dim[, axis.y],
        col = col,
        do.lm = F,
        xlab = xlab,
        ylab = ylab,
        ...
    ))

    if (is.numeric(col)) {
        p <- p + scale_color_viridis_c(name = ident)
    } else {
        p <- p + scale_color_viridis_d(name = ident)
    }

    return(p)
}



#-------------------------------------------------------------------------------
#' Plot a hexbin colored on a 3rd variable
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param xlab xlab
#' @param ylab ylab
#' @param zlab Scale title
#' @param main title
#' @param bins bins
#' @param scale.to.well.mean Divide z by the mean (fold change compared to mean)
#' @param scale.z Scale z to mean 0 variance 1
#' @param trim.outliers.z Remove outliers in z prior to plotting using modified z-score
#' @param trim.outliers.z.thresh Modified z-score hreshold to consider an outlier
#' @param bins.mincount Minimum number of objects in a bin to render it
#'
#' @returns A ggplot2 object
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab stat_summary_hex scale_fill_viridis_c
#' @export
tglow_plot_location_hex <- function(dataset,
                                    assay,
                                    slot,
                                    feature.z,
                                    feature.x = "nucl_AreaShape_Center_X",
                                    feature.y = "nucl_AreaShape_Center_Y",
                                    xlab = "X centroid",
                                    ylab = "Y centroid",
                                    zlab = "Z",
                                    main = NULL,
                                    bins = 30,
                                    scale.to.well.mean = F,
                                    scale.z = F,
                                    trim.outliers.z = F,
                                    trim.outliers.z.thresh = 3.5,
                                    bins.mincount = 0) {
    check_dataset_assay_slot(dataset, assay, slot)


    # Build the plot
    df.plot <- getDataByObject(dataset, c(feature.z, feature.x, feature.y), assay = assay, slot = slot)

    colnames(df.plot) <- c("z", "x", "y")

    if (trim.outliers.z) {
        df.plot <- df.plot[abs(mod_zscore(df.plot$z)) < trim.outliers.z.thresh, ]
    }

    if (scale.to.well.mean) {
        df.plot$z <- df.plot$z / mean(df.plot$z)
    }

    if (scale.z) {
        df.plot$z <- scale(df.plot$z)
    }

    if (is.null(main) && is.character(feature.z)) {
        main <- feature.z
    }


    p1 <- ggplot(
        data = df.plot,
        mapping = aes(x = x, y = y, z = z)
    ) +
        stat_summary_hex(
            fun = function(x) {
                if (length(x) > bins.mincount) {
                    mean(x)
                } else {
                    return(NA)
                }
            },
            bins = bins
        ) +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(main) +
        scale_fill_viridis_c(name = zlab)




    return(theme_plain(p1))
}

#-------------------------------------------------------------------------------
#' Plot a hexbin scatterplot colored on density
#'
#' @param x x
#' @param y y
#' @param bins bins
#' @param do.lm Should a linear fit be rendered
#' @param lm.col Color of the line
#' @param xlab xlab
#' @param ylab ylab
#' @param main title
#' @param facet Vector to facet on
#' @param ... Arguments to facet_wrap
#'
#' @returns A ggplot2 object
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab geom_smooth geom_hex facet_wrap
#' @export
plot_hex <- function(x, y, bins = 250, do.lm = T, lm.col = "lightgrey", xlab = "x", ylab = "y", main = NULL, facet = NULL, ...) {
    df.plot <- data.frame(x = x, y = y)

    if (!is.null(facet)) {
        df.plot$facet <- facet
    }

    p1 <- ggplot(data = df.plot, mapping = aes(x = x, y = y)) +
        geom_hex(bins = bins) +
        xlab(xlab) +
        ylab(ylab)

    if (do.lm) {
        p1 <- p1 + geom_smooth(method = "lm", col = lm.col)
    }

    # Add title
    if (is.null(main) & do.lm) {
        main <- paste0(
            main.prefix, "R: ", format(ct$estimate, digits = 2),
            " p-value: ", format(ct$`p.value`, digits = 2, scientific = T)
        )
    } else {
        main <- NA
    }
    p1 <- p1 + ggtitle(p1)

    if (!is.null(facet)) {
        p1 <- p1 + facet_wrap(~facet, ...)
    }

    return(p1)
}

#-------------------------------------------------------------------------------
#' Plot an individual EBImage image with a title
#'
#' @param img The image
#' @param main title
#' @param marker.add Should a marker be added at the center of the image
#' @param marker.x x position of marker. If NULL = center
#' @param marker.y y position of marker. If NULL = center
#' @param marker.col color of the marker
#' @param marker.size size of marker
#' @param marker.shape shape of marker
#'
#' @returns A ggplot2 object
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab annotation_raster geom_point coord_fixed theme_void
#' @export
plot_img <- function(img, main = "img", marker.add = T, marker.x = NULL, marker.y = NULL, marker.col = "white", marker.size = 8, marker.shape = 18) {
    p <- ggplot() +
        ggtitle(main) +
        annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        coord_fixed()

    if (marker.add) {
        if (is.null(marker.x)) {
            marker.x <- round(dim(img)[1] / 2)
        }

        if (is.null(marker.y)) {
            marker.y <- round(dim(img)[2] / 2)
        }

        marker.x

        p <- p + geom_point(
            mapping = aes(x = marker.x, y = marker.y),
            col = marker.col,
            shape = marker.shape,
            size = marker.size
        )
    }
    return(p + theme_void())
}


#-------------------------------------------------------------------------------
#' Plot a series of EBImage images with a title
#'
#' @param imgs List of images
#' @param ncol Numer of columns in grid
#' @param main Main plot title
#' @param main.sub Vector of length(imgs) indicating the subplot titles
#' @param text.col Color of the text
#' @param background.col Color of the background
#' @param byrow Should plots be filled per row or per column
#' @param ... Remaining arguments passed to \code{\link{plot_img}}
#'
#' @returns A ggplot2 object
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom ggplot2 element_rect theme
#' @export
plot_img_set <- function(imgs, ncol, main = "", main.sub = NULL, text.col = "white", background.col = "black", byrow = F, ...) {
    plots <- list()
    i <- 0
    for (img in seq_along(imgs)) {
        i <- i + 1
        if (!is.null(main.sub)) {
            p <- plot_img(imgs[[img]], main = main.sub[i], ...)
        } else {
            p <- plot_img(imgs[[img]], main = "", ...)
        }

        p <- p + theme(plot.title = element_text(colour = text.col)) +
            theme(
                plot.background = element_rect(fill = background.col),
                plot.title = element_text(hjust = 0.5)
            )
        plots[[img]] <- p
    }

    p.tmp <- plot_grid(plotlist = plots, ncol = ncol, byrow = byrow) +
        theme(plot.background = element_rect(fill = background.col, colour = NA))

    title <- ggdraw() + draw_label(main, fontface = "bold")

    p.final <- plot_grid(title, p.tmp, ncol = 1, rel_heights = c(0.05, 1))

    return(p.final)
}


#------------------------------------------------------------------------------
#' Scatterplot using ggplot
#'
#' @param x x
#' @param y y
#' @param xlab xlab
#' @param ylab ylab
#' @param main Main plot title, overrides default correlation
#' @param main.prefix Text to add to plot title before correlation
#' @param size size
#' @param col col
#' @param fixed Should axes be fixed to the same x and y
#' @param alpha alpha
#' @param shape shape
#' @param lm.col Color of the line
#' @param do.lm Should a linear fit be rendered
#' @param method Method for fitting the line. Anything supported by \code{\link{ggplot2::geom_smooth()}}
#' @param lm.group Group to fit model in
#' @param raster Should points be rasterized (usefull for large numbers of points)
#' @param dpi When raster=TRUE the DPI to render at
#' @param facet Vector with facet variable
#' @param facet.ncol Number of columns in facet
#'
#' @returns A ggplot2 object
#' @importFrom stats cor.test
#' @importFrom ggrastr rasterize
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab annotation_raster geom_point coord_fixed geom_smooth geom_abline facet_wrap xlim ylim
#' @export
plot_xy <- function(x, y, xlab = "X", ylab = "Y", main = NA, main.prefix = "", size = 1, col = "black", fixed = F, alpha = 0.75, shape = 16, lm.col = "blue", do.lm = T, method = "lm", lm.group = NULL, raster = F, dpi = 300, facet = NULL, facet.ncol = NULL) {
    df.plot <- data.frame(
        x = x,
        y = y
    )

    # Fix axis limits between x and y
    if (fixed) {
        lims <- c(min(c(x, y), na.rm = T), max(c(x, y), na.rm = T))
    }

    # Facet y or n
    if (!is.null(facet)) {
        df.plot$facet <- facet
    }

    if (nrow(df.plot) > 3) {
        ct <- cor.test(df.plot$x, df.plot$y)
    } else {
        ct <- list(`p.value` = NA, `estimate` = NA)
    }

    if (length(shape) > 1 && length(col > 1)) {
        if (raster) {
            p <- ggplot(data = df.plot, mapping = aes(x = x, y = y, col = col, shape = shape)) +
                rasterize(geom_point(alpha = alpha, size = size), dpi = dpi)
        } else {
            p <- ggplot(data = df.plot, mapping = aes(x = x, y = y, col = col, shape = shape)) +
                geom_point(alpha = alpha, size = size)
        }
    } else if (length(col) > 1) {
        if (raster) {
            p <- ggplot(data = df.plot, mapping = aes(x = x, y = y, col = col)) +
                rasterize(geom_point(alpha = alpha, shape = shape, size = size), dpi = dpi)
        } else {
            p <- ggplot(data = df.plot, mapping = aes(x = x, y = y, col = col)) +
                geom_point(alpha = alpha, shape = shape, size = size)
        }
    } else {
        if (raster) {
            p <- ggplot(data = df.plot, mapping = aes(x = x, y = y)) +
                rasterize(geom_point(alpha = alpha, color = col, shape = shape, size = size), dpi = dpi)
        } else {
            p <- ggplot(data = df.plot, mapping = aes(x = x, y = y)) +
                geom_point(alpha = alpha, color = col, shape = shape, size = size)
        }
    }

    p <- p +
        xlab(xlab) +
        ylab(ylab)

    # Add linear line
    if (do.lm) {
        p <- p + geom_smooth(method = method, col = lm.col, inherit.aes = F, mapping = aes(x = x, y = y, group = lm.group))
    }

    # If axis are fixed, add an abline and set limits
    if (fixed) {
        p <- p +
            xlim(lims) +
            geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
            coord_fixed() +
            ylim(lims)
    }

    # Add title
    if (is.na(main) & do.lm) {
        main <- paste0(
            main.prefix, "R: ", format(ct$estimate, digits = 2),
            " p-value: ", format(ct$`p.value`, digits = 2, scientific = T)
        )
    }

    if (!is.na(main)) {
        p <- p + ggtitle(main)
    }

    if (!is.null(facet)) {
        p <- p + facet_wrap(~facet, ncol = facet.ncol)
    }

    return(p)
}

#------------------------------------------------------------------------------
#' Plot a grouped histogram
#'
#' @param x x
#' @param z grouping variable
#' @param xlab xlab
#' @param main title
#' @param bins bins
#' @param facet Vector with facet variable
#' @param facet.ncol Number of columns in facet
#' @param density Plot as a density plot instead
#'
#' @returns A ggplot2 object
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab geom_vline geom_histogram scale_fill_viridis_d facet_wrap geom_density scale_color_viridis_d
#' @export
plot_hist_dens_grouped <- function(x,
                                   z = NULL,
                                   xlab = "Value",
                                   main = NULL,
                                   bins = 100,
                                   facet = NULL,
                                   facet.ncol = NULL,
                                   density = FALSE) {
    if (is.null(z)) {
        z <- "1"
    }

    # Build the plot
    df.plot <- data.frame(
        x = x,
        z = z
    )

    # Facet y or n
    if (!is.null(facet)) {
        df.plot$facet <- facet
    }

    med <- median(df.plot$x, na.rm = TRUE)

    if (!density) {
        p1 <- ggplot(
            data = df.plot,
            mapping = aes(x = x, fill = z), alpha = 0.5, color = "black"
        ) +
            geom_histogram() +
            geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) + # Add vertical line for the median
            xlab(xlab) +
            ggtitle(main) +
            scale_fill_viridis_d()
    } else {
        p1 <- ggplot(
            data = df.plot,
            mapping = aes(x = x, col = z), alpha = 0.5, color = "black"
        ) +
            geom_density() +
            geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) + # Add vertical line for the median
            xlab(xlab) +
            ggtitle(main) +
            scale_color_viridis_d()
    }

    if (!is.null(facet)) {
        p <- p + facet_wrap(~facet, ncol = facet.ncol)
    }
    return(p1)
}

#-------------------------------------------------------------------------------
#' Plot boxplot with a numerically ordered x and a loess line
#'
#' @param x x
#' @param y y
#' @param levels The order of the levels on the x-axis
#' @param xlab xlab
#' @param ylab ylab
#' @param main title
#' @param line.col Color of the line
#' @param facet Vector with facet variable
#' @param facet.ncol Number of columns in facet
#' @param method Method for fitting the line. Anything supported by \code{\link{ggplot2::geom_smooth()}}
#'
#' @returns A ggplot2 object
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab scale_x_continuous facet_wrap geom_boxplot geom_smooth
#' @export
plot_boxline <- function(x, y, levels, xlab = "x", ylab = "y", main = "", line.col = "#3d403d", facet = NULL, facet.ncol = NULL, method = NULL) {
    df.plot <- data.frame(x = as.numeric(factor(x, levels = levels)), y = y)
    if (!is.null(facet)) {
        df.plot$facet <- facet
    }

    p1 <- ggplot(data = df.plot, mapping = aes(x = x, y = y, group = x, facet = facet)) +
        geom_boxplot() +
        geom_smooth(se = F, aes(group = 1), col = line.col, method = method) +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(main) +
        scale_x_continuous(breaks = 1:length(levels), labels = unique(levels))

    if (!is.null(facet)) {
        p1 <- p1 + facet_wrap(~facet, ncol = facet.ncol)
    }

    return(p1)
}


#------------------------------------------------------------------------------
#' Simple heatmap with auto labels
#'
#' @param data A matrix
#' @param cellsize The size of the cells in x and y
#' @param cellwidth The size of the cells in x
#' @param cellheight THe size of the cells in y
#' @param limits Clip the data to upper and lower, must be a numeric of length 2
#' @param cluster Should clustering be applied
#' @param range Automatic setting for color range: 'symmetric', 'absolute', 'auto'
#' @param palette Palette to use in \code{\link{RColorBrewer::brewer.pal()}}
#' @param border Border color
#' @param ... Remaining arguments passed to  \code{\link{pheatmap::pheatmap()}}
#'
#' @returns The pheatmap plot
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
plot_simple_hm <- function(data, cellsize = -1, cellwidth = 12, cellheight = 12, limits = NULL, cluster = T, range = "symmetric", palette = NULL, border = NA, ...) {
    if (cellsize > 0) {
        cellwidth <- cellsize
        cellheight <- cellsize
    }

    if (is.null(limits)) {
        max <- max(data)
        min <- min(data)
        max.abs <- max(abs(data))
        min.abs <- min(abs(data))
    } else {
        if (length(limits) == 2 && is(limits, "numeric")) {
            max <- limits[2]
            min <- limits[1]
            max.abs <- max(abs(limits))
            min.abs <- min(abs(limits))
            range <- "auto"
        } else {
            stop("Limits must be a vector of length 2")
        }
    }

    if (range == "symmetric") {
        break.list <- seq(-max.abs, max.abs, by = max.abs / 50)
        if (is.null(palette)) {
            palette <- "RdBu"
        }
        cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = palette)))(length(break.list))
    } else if (range == "absolute") {
        if (is.null(palette)) {
            palette <- "Reds"
        }
        break.list <- seq(min, max.abs, by = (max.abs - min.abs) / 100)
        cols <- grDevices::colorRampPalette(c("#FFFFFF", RColorBrewer::brewer.pal(n = 7, name = palette)))(length(break.list))
    } else if (range == "auto") {
        break.list <- seq(min, max, by = (max - min) / 100)
        if (is.null(palette)) {
            palette <- "RdBu"
        }
        cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = palette)))(length(break.list))
    } else {
        stop("Range must be 'symmetric', 'auto', or 'absolute'\n")
    }

    if (!cluster) {
        pheatmap::pheatmap(data,
            breaks = break.list,
            col = cols,
            cellwidth = cellwidth,
            cellheight = cellheight,
            border = border,
            cluster_rows = F,
            cluster_cols = F,
            ...
        )
    } else {
        pheatmap::pheatmap(data,
            breaks = break.list,
            col = cols,
            cellwidth = cellwidth,
            cellheight = cellheight,
            border = border,
            ...
        )
    }
}
