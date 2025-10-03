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
    
    if (sum(out <= 0) > 0 ) {
      out[out <=0] <- max(out):((max(out)+sum(out <=0))-1)
      warning("Ties detected leading to negative indices. Returning the closet as per base::sort() tie breaking")
    }
    #out <- out[out > 0 & out <= length(out)]
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
#' @param q A vector of quantiles [0-100]
#' @param add.extremes Should the top and bottom n values be added
#'
#' @returns A list with row ids and labels for the objects
#' @export
fetch_representative_object_quantiles <- function(dataset, assay, slot, feature, name = NULL, n = 1, q=c(10, 25, 50, 75, 90), add.extremes=TRUE) {
  check_dataset_assay_slot(dataset, assay, slot)

  if (length(feature) > 1) {
    stop("Can only supply one feature at the time")
  }

  if (is.null(name)) {
    name <- feature
  }
  
  f   <- getDataByObject(dataset, feature, assay=assay, slot=slot)
  idx <- dataset@object.ids[!is.na(f)]
  f   <- f[!is.na(f)]
  fs  <- order(f)

  # Store the output
  ids <- c()
  cn  <- c()

  if (add.extremes) {
      lower <- which(dataset@object.ids %in% idx[head(fs, n = (n * 2) + 1)])
      ids   <- c(ids, lower)
      cn    <- c(cn, rep(paste0("q0 - ", name), (n * 2) + 1)) 
  }

  for (cur.q in q) {
    l <- fetch_representative_object(dataset, assay, slot, feature,
      metric = "lower.q",
      q = cur.q/100,
      n = n
    )
    ids  <- c(l, ids)
    cn   <- c(cn, rep(paste0("q", cur.q, " - ", name), (n * 2) + 1)) 
  }
  
  if (add.extremes) {
      upper <- which(dataset@object.ids %in% idx[tail(fs, n = (n * 2) + 1)])
      ids   <- c(ids, upper)
      cn    <- c(cn, rep(paste0("q100 - ", name), (n * 2) + 1)) 
  }
  
  return(list(ids = ids, names = cn))
}


#-------------------------------------------------------------------------------
#' Retrieve a cell (and its neighbours) based on a feature sumstat
#'
#' @description
#' Gets a objects and its closes neighbours based on a reduction and a feature.
#' Feature should be catagorical, can be used to find the object that is closest to the centroid
#' of a group in reduction space. To get more than one object per group, specify n.
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
  
  if (class(x) == "numeric") {
    stop("Currently does not support numeric features, if it is categorical numeric, convert to a string/factor first")
  }
  
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