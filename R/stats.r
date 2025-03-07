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
#' Calculate the skewness, kurtosis
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
