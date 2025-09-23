#-------------------------------------------------------------------------------
#' Caclulate the effective dimensionality of a dataset
#'
#' @description Calculates the effective dimensionality of a multidimensional dataset.
#'
#' @param data A matrix with data
#' @param method 'LiJi', 'pc_var', 'Chevrud', 'Galwey' . See details
#' @param var.test Percentage of variance cutoff
#' @details
#' `method`
#'
#' `eigenval` are the eigenvalues of the correlation matrix.
#'
#' `m` are the number of eigenvalues
#'
#' - `LiJi` Use Li and Ji's method - https://www.nature.com/articles/jhg201134#Sec2
#'
#'    `eff.tests <- sum(I( eigenval > 1) + ( eigenval - floor( eigenval)))`
#'
#' - `pc_var` Number of PC's needed to pass var.thresh
#'
#'     `eff.tests <- which(cumsum(eigenval / sum(eigenval)) >= var.thresh)`
#'
#' - `Chevrud`:
#'
#'     `eff.tests <- m +(m - 1) * (1-(var(eigenval) / m))`
#'
#' - `Galwey`
#'
#' @export
effective_dimensionality <- function(data, method = "LiJi", var.thresh = 0.95) {
  if (!is.matrix(matrix)) {
    stop("data should be of class matrix")
  }
  c <- cor(matrix, use = "pairwise.complete.obs")
  e <- eigen(c, only.values = T)

  eigenval <- e$values
  m <- length(eigenval)

  if (method == "LiJi") {
    eff.tests <- sum(I(eigenval > 1) + (eigenval - floor(eigenval)))
  } else if (method == "pc_var") {
    eff.tests <- which(cumsum(eigenval / sum(eigenval)) >= var.thresh)[0]
  } else if (method == "Cheverud") {
    eff.tests <- m + (m - 1) * (1 - (var(eigenval) / m))
  } else if (method == "Galwey") {
    eff.tests <- sum(eigenval)^2 / sum(eigenval^2)
  } else {
    stop(paste0(method, " is not a valid method"))
  }

  return(eff.tests)
}

#-------------------------------------------------------------------
#' Calculate the skewness, kurtosis, geometric mean, mode
#' 
#' @param x Numeric vector
#' @param na.rm Should NA's be removed
#' @rdname moments
#' @returns estimates for skewness, kurtosis respecitvely
#' @export
skewness <- function(x, na.rm=F) {
  if (na.rm){
    x <- x[!is.na(x)]
  }
  
  sum((x-mean(x))^3)/((length(x)-1)*sd(x)^3)
}

#' @rdname moments
#' @export
kurtosis <- function(x, na.rm=F) {
  if (na.rm){
    x <- x[!is.na(x)]
  }
  sum((x-mean(x))^4)/((length(x)-1)*sd(x)^4)
}


#' @rdname moments
#' @export
geometric_mean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' @rdname moments
#' @export
mode <- function(x, na.rm=TRUE) {
  
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  xu <- unique(x)
  xu[which.max(tabulate(match(x, xu)))]
}


#-------------------------------------------------------------------------------
#' Calculate summary stats on features
#'
#' @description
#' Calculate the mean, median, mode, skewness, kurtosus, geometric mean, min, max
#' Removes NA's for each feature individually.
#' Stats are put on the assay@features slot under the name <slot>_<stat>
#' 
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#'
#' @returns The \linkS4class{TglowDataset}
#' @export
calculate_sumstats <- function(dataset, assay, slot) {
  
  check_dataset_assay_slot(dataset, assay, slot)
  
  data <- slot(dataset[[assay]], slot)

  dataset@assays[[assay]]@features[,paste0(slot, "_", "mean")] <- base::colMeans(data, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "median")] <- matrixStats::colMedians(data, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "mode")] <- apply(data, 2, mode, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "skew")] <- apply(data, 2, skewness, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "kurt")] <- apply(data, 2, kurtosis, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "geom.mean")] <- apply(data, 2, geometric_mean, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "min")] <- apply(data, 2, min, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "max")] <- apply(data, 2, max, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "var")] <- matrixStats::colVars(data, na.rm=T)
  dataset@assays[[assay]]@features[,paste0(slot, "_", "unique")] <- apply(data, 2, function(x){length(unique(x[!is.na(x)]))})

  return(dataset)
  
}
