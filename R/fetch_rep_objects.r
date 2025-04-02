#-------------------------------------------------------------------------------
#' Retrieve a cell (and its neighbours) based on a feature sumstat
#'
#' @description
#' Gets a objects and its closes neighbours based on a single feature and a
#' sumstat (mean, median, upper.q, lower.q). To customise the quantile used
#' specify q.
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param feature The feature to find a representative objects for
#' @param metric Can be 'mean', 'median', 'upper.q', 'lower.q'
#' @param na.rm Should NA's be removed
#' @param n How many objects either side of the representitative objects should be returend
#' @param subset Look only in a subset of objects. Must be a selection vector
#' @param q Override the quantile when using upper.q (0.75) or lower.q (0.25)
#'
#' @returns vector of indices in objects matrix
#' @export
fetch_representative_object <- function(dataset, assay, slot, feature, metric = "mean", na.rm = F, n = 0, subset = NULL, q = NULL) {
  check_dataset_assay_slot(dataset, assay, slot)

  x <- getDataByObject(dataset, feature, assay = assay, slot = slot, drop = T)
  # x <- slot(dataset@assays[[assay]], slot)[,feature]
  i <- seq_len(length(x))

  if (!is.null(subset)) {
    x <- x[subset]
    i <- i[subset]
  }

  i <- i[order(x)]
  x <- sort(x)

  if (metric == "mean") {
    m <- mean(x, na.rm = na.rm)
    out <- nearest_index(x, m)
  } else if (metric == "median") {
    m <- median(x, na.rm = na.rm)
    out <- nearest_index(x, m)
  } else if (metric == "upper.q") {
    if (is.null(q)) {
      q <- 0.75
    }
    m <- quantile(x, probs = q)
    out <- nearest_index(x, m)
  } else if (metric == "lower.q") {
    if (is.null(q)) {
      q <- 0.25
    }
    m <- quantile(x, probs = q)
    out <- nearest_index(x, m)
  }

  if (n > 0) {
    out <- (out - n):(out + n)
  }

  return(i[out])
}


#-------------------------------------------------------------------------------
#' Retrieve a set of objects (and its neighbours) based on a feature sumstat
#'
#' @description
#' Fetches (n*2)+1 objects arround the 0th, 10th, 25th, 50th, 75th, 90th and 100th quantiles
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param feature The feature to find a representative cell for
#' @param name Prefix to add to the names of the names vector
#' @param n How many objects either side of the representitative cell should be returend
#'
#' @returns A list with row ids and labels for the objects
#' @export
fetch_representative_object_quantiles <- function(dataset, assay, slot, feature, name = NULL, n = 1) {
  check_dataset_assay_slot(dataset, assay, slot)

  if (is.null(name)) {
    name <- feature
  }
  #cur.assay <- slot(dataset@assays[[assay]], slot)@.Data
  #cur.assay <- slot(dataset@assays[[assay]], slot)
  
  #f <- cur.assay[, feature]
  
  f <- getDataByObject(dataset, feature, assay=assay, slot=slot)
  
  idx <- dataset@object.ids[!is.na(f)]
  f <- f[!is.na(f)]
  fs <- order(f)

  lower <- which(dataset@object.ids %in% idx[head(fs, n = (n * 2) + 1)])
  upper <- which(dataset@object.ids %in% idx[tail(fs, n = (n * 2) + 1)])

  l10 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.10,
    n = n
  )

  l25 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.25,
    n = n
  )

  l50 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.5,
    n = n
  )

  l75 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.75,
    n = n
  )

  l90 <- fetch_representative_object(dataset, assay, slot, feature,
    metric = "lower.q",
    q = 0.90,
    n = n
  )

  cn <- c(
    rep(paste0("q0 - ", name), (n * 2) + 1),
    rep(paste0("q10 - ", name), (n * 2) + 1),
    rep(paste0("q25 - ", name), (n * 2) + 1),
    rep(paste0("q50 - ", name), (n * 2) + 1),
    rep(paste0("q75 - ", name), (n * 2) + 1),
    rep(paste0("q90 - ", name), (n * 2) + 1),
    rep(paste0("q100 - ", name), (n * 2) + 1)
  )

  return(list(ids = c(lower, l10, l25, l50, l75, l90, upper), names = cn))
}


#-------------------------------------------------------------------------------
#' Retrieve a cell (and its neighbours) based on a feature sumstat
#'
#' @description
#' Gets a objects and its closes neighbours based on a single feature and a
#' sumstat (mean, median, upper.q, lower.q). To customise the quantile used
#' specify q.
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay to use
#' @param slot The slot to use for calculating filters, defaults to "data". Can be "data" or "scale.data"
#' @param feature The feature to find a representative objects for, used as a grouping variable
#' @param reduction The reduction to use to calculate distance. Reccomend a PCA.
#' @param assay if feature is on an assay, the assay to use
#' @param slot if feature is on an assay, which slot to use
#' @param n How many nearest neighbours to the central most object should be fetched
#'
#' @returns vector of indices in objects matrix
#' @export
fetch_representative_object_nn <- function(dataset, feature, reduction, assay=NULL, slot=NULL, n=0) {
  check_dataset_assay_slot(dataset, assay, slot)
  
  x <- getDataByObject(dataset, feature, assay = assay, slot = slot, drop = F)
  
  if(!reduction %in% names(dataset@reduction) ){
    stop(paste0(reduction, " is not available on dataset@reduction"))
  }
  
  red <- dataset@reduction[[reduction]]

  cells  <- c()
  groups <- c()
  
  # Find for each group, the closest cell to the mean in reduction space  
  for (cur.group in unique(x[,1])) {
    
    # Find the cell in the group closest to the mean
    tmp         <- red@x[x==cur.group,]
    tmp         <- tmp - colMeans(tmp)
    dist        <- apply(tmp, 1, function(x){ sum(x^2)})
    
    cur.cell <- names(which(dist==min(dist)))

    # Find its nearest neighbours
    if (n > 0) {
      cur.knn <- RANN::nn2(tmp, tmp[cur.cell,, drop=F],  k=n)
      cells <- c(cells, rownames(tmp)[cur.knn$nn.idx])
      groups <- c(groups, rep(cur.group, n))
    } else {
      cells <- c(cells,cur.cell)
      groups <- c(groups, cur.group)
    }
  }
  
  return(list(ids=cells, names=groups))
  
}