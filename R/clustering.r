#-------------------------------------------------------------------------------
#' Cluster using a K nearest neighbor graph in the PCA space
#'
#' @description
#' Builds a k-nearest neighbor graph from the supplied PCA reduction and clusters
#' it using Louvain or Leiden community detection. Clusters are stored as new
#' columns in the meta slot of the returned dataset.
#'
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param reduction The reduction to use for calculating clustering. If NULL re-calculated
#' @param pc.n How many PC's to use
#' @param k How many NN to calculate
#' @param method Clustering method to use 'louvain' or 'leiden'
#' @param resolution Resolution for method 'louvain' or 'leiden'. Can be a vector of resolutions
#' @param exact.nn Instead of using Seurat ANNOY for kNN use \code{\link[=nn2]{RANN::nn2()}}
#' @param col.out Prefix for column name to store the clustering under added to the meta slot on the output object. <col.out>_res_<res>
#'
#' @details
#' Here I use the Seurat implementation of the kNN, which is NOT and exact
#' kNN graph, but it is very fast. Practically it should perform pretty well
#' The parameter k controls how many nearest neighbors to find
#' The distance matrix is built from the principal components
#'
#'
#' Alternatively:
#' Build exact knn graph, better but MUCH slower
#' Set eps to something not 0 to allow error tolerance and
#' get an approximate nn
#'
#'
#' Cluster the graph using louvain (default) or leiden clustering
#' The resolution parameter controls how many clusters get generated
#' you will need to play with this to find a reasonable number of clusters
#' lower resolution tends to result in fewer clusters
#'
#' @returns The \linkS4class{TglowDataset} with clusters added to meta slot and the graph
#' @importFrom igraph graph_from_edgelist membership cluster_louvain cluster_leiden
#' @importFrom RANN nn2
#' @export
calculate_clustering <- function(dataset, reduction, pc.n = NULL, k = 10, method = "louvain", resolution = 0.1, exact.nn = FALSE, col.out = "clusters_") {
    if (!reduction %in% names(dataset@reduction)) {
        stop(paste0("Reduction ", reduction, " not found in dataset"))
    }

    if (!is.null(pc.n)) {
        if (ncol(dataset@reduction[[reduction]]@x) < pc.n) {
            stop(paste0("Reduction ", reduction, " does not have ", pc.n, " PC's available"))
        }
    } else {
        pc.n <- ncol(dataset@reduction[[reduction]]@x)
    }

    pcs <- dataset@reduction[[reduction]]@x[, 1:pc.n, drop = FALSE]

    # Re-use graph if available
    if (!is.null(dataset@graph) && dataset@graph$k == k && dataset@graph$reduction == reduction && dataset@graph$exact.nn == exact.nn) {
        cat("[INFO] Re-using existing graph\n")
        graph <- dataset@graph$graph
    } else {
        # Build  knn graph
        if (exact.nn) {
            # RANN KNN impelemntation
            knn <- nn2(pcs, k = k, eps = 0)
        } else {
            # Seurat KNN approximate implementation
            knn <- AnnoyNN(pcs, k = k)
        }

        # Next convert the kNN to an edge list for clustering
        # Create edge list
        edges <- matrix(0, nrow = nrow(pcs) * k, ncol = 2)
        for (i in seq_len(nrow(pcs))) {
            edges[((i - 1) * k + 1):(i * k), 1] <- i
            edges[((i - 1) * k + 1):(i * k), 2] <- knn$nn.idx[i, ]
        }

        graph <- graph_from_edgelist(edges, directed = FALSE)

        dataset@graph <- list(graph = graph, k = k, reduction = reduction, exact.nn = exact.nn)
    }


    for (res in resolution) {
        cat("[INFO] Calculating clustering for resolution: ", res, "\n")
        cur.out <- paste0(col.out, "res_", res)
        if (method == "louvain") {
            cl <- cluster_louvain(graph, resolution = res)
        } else if (method == "leiden") {
            cl <- cluster_leiden(graph, resolution = res)
        } else {
            stop("No valid cluster method")
        }

        # Make a vector of the cluster memberships for each cell
        dataset@meta[, cur.out] <- as.character(membership(cl))
    }


    return(dataset)
}

#-------------------------------------------------------------------------------
#' Calculate feature clustering using correlation-based graph or hierarchical clustering
#'
#' @description
#' Clusters features based on their correlation patterns using either graph-based 
#' community detection or hierarchical clustering methods.
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay on dataset to use
#' @param slot The assay slot to use ("data", "scale.data")
#' @param feature.group Column name in features slot or vector of length ncol(assay) 
#'        defining feature groups to cluster separately
#' @param resolution Resolution parameter for clustering. Controls cluster granularity.
#'        For graph methods, lower values give fewer clusters. For hierarchical methods,
#'        represents height to cut the tree.
#' @param threshold Correlation threshold. Only correlations above this value are kept
#' @param correlation Pre-calculated correlation matrix (optional). If feature.group is 
#'        provided, should be a named list of correlation matrices matching group levels.
#'        Otherwise, a single correlation matrix.
#' @param col.out Prefix for output column names in features slot
#' @param method Clustering method to use: "louvain", "leiden", "complete", "ward.D", or "ward.D2"
#' @param return.object If TRUE returns updated dataset, if FALSE returns list with 
#'        clusters and correlation matrices
#'
#' @details
#' For graph-based methods (louvain/leiden):
#' Creates a graph where nodes are features and edges are correlations above the threshold.
#' The graph is then clustered using community detection.
#'
#' For hierarchical methods (complete/ward):
#' Uses 1-correlation as distance metric for hierarchical clustering.
#' The resolution parameter determines where to cut the dendrogram.
#' When method is complete, the resolution can be seen as the maximum
#' difference in correlation within a cluster between the highest and lowest correlated features.
#' So a value of 0.1 would mean all features in a cluster have at least 0.9 correlation with each other.
#'
#' Features can be clustered in groups by providing feature.group. This is useful
#' when different feature types should be clustered separately.
#'
#' @return If return.object=TRUE, returns the \linkS4class{TglowDataset} with 
#'         feature clusters added to features slot. If FALSE, returns list with
#'         cluster assignments and correlation matrices.
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain cluster_leiden membership
#' @importFrom stats hclust cutree cor dist as.dist
#' @export
calculate_feature_clustering <- function(
    dataset,
    assay,
    slot,
    feature.group=NULL,
    resolution = 0.1,
    threshold = 0,
    correlation = NULL,
    col.out = "fcl_",
    method = "complete",
    return.object=TRUE
) {
        
    check_dataset_assay_slot(dataset, assay, slot)
    
    # Validate method
    valid_methods <- c("louvain", "leiden", "complete", "ward.D", "ward.D2")
    if (!method %in% valid_methods) {
        stop(paste0("Invalid method. Must be one of: ", paste(valid_methods, collapse=", ")))
    }

    data         <- slot(dataset@assays[[assay]], slot)
    corr.list    <- list()
    cluster.cols <- c()
    graph.list   <- list()

    if (is.null(feature.group)) {
        grouping <- rep("all", ncol(dataset@assays[[assay]]))
        
        if (!is.null(correlation)) {
            corr.list[["all"]] <- correlation
        }
    } else {
        if (length(feature.group) == 1) {
            if (!feature.group %in% colnames(dataset@assays[[assay]]@features)) {
                stop("feature.group must be a column on assay@features or a vector of ncol(assay)")
            }
            grouping <- dataset@assays[[assay]]@features[,feature.group]
        } else {
            if (length(feature.group) != ncol(dataset@assays[[assay]])) {
                stop("feature.group must be a column on assay@features or a vector of ncol(assay)")
            }
            grouping <- feature.group
        }
        
        if (!is.null(correlation) && (!is.list(correlation) || length(correlation) != length(unique(grouping)))) {
            stop("When feature.group is set, correlation must be a list matching the levels of feature.group")
        }
        
        if (!is.null(correlation)) {
            corr.list <- correlation
            
            if (!all(unique(grouping) %in% names(corr.list))) {
                stop(paste0("Not all grouping levels found in supplied correlation list"))
            }
        }
    }

    for (group in unique(grouping)) {
        cat("[INFO] Running feature clustering for group: ", group, "\n")
        selector <- grouping == group
        
        # Calculate correlation matrix if not provided
        if (is.null(corr.list[[group]])) {
            cat("[INFO] Running correlation for group: ", group, "\n")
            correlation <- cor(data[,selector, drop=F], use = "pairwise.complete.obs")
            
            if (!return.object) {
                corr.list[[group]] <- correlation
            }
        } else {
            correlation <- corr.list[[group]]
        }
        
        # Threshold correlation matrix
        correlation[is.na(correlation)] <- 0
        correlation[correlation < threshold] <- 0
    
        if (method %in% c("louvain", "leiden")) {

            if (threshold < 0) {
                stop("Threshold must be non-negative for graph-based methods.")
            }

            # Create igraph object from adjacency matrix
            graph <- igraph::graph_from_adjacency_matrix(
                correlation,
                mode = "undirected",
                weighted = TRUE,
                diag = FALSE
            )

            graph.list[[group]] <- graph
        } else if (method %in% c("complete","ward.D", "ward.D2")) {
            # A use 1 - correlation as a distance measure
            dist_mat <- as.dist(1 - correlation)
            
            # Perform hierarchical clustering
            hc <- hclust(dist_mat, method = method)
        }
        
        # Calculate clustering for each resolution
        for (res in resolution) {
            cur.out            <- paste0(col.out, method,    "_res_", res)
            cluster.cols <- c(cluster.cols,    paste0(col.out, method, "_group"), cur.out)
            cat("[INFO] Calculating clustering for: ", cur.out, "\n")
            
            if (method == "louvain") {
                cl <- igraph::cluster_louvain(graph, resolution = res)
                clusters <- as.character(igraph::membership(cl))
            } else if (method == "leiden") {
                cl <- igraph::cluster_leiden(graph, resolution = res)
                clusters <- as.character(igraph::membership(cl))
            } else if (method %in% c("complete", "ward.D", "ward.D2")) {
                clusters <- cutree(hc, h=res)
            } else {
                stop("No valid cluster method")
            }
            
            # Make a vector of the cluster memberships for each cell
            dataset@assays[[assay]]@features[selector, cur.out] <- paste0(group, "__", clusters)
            dataset@assays[[assay]]@features[selector, paste0(col.out, method, "_group")] <- group
        }
    }
    
    if (!return.object) {
        return(list(
            clusters=dataset@assays[[assay]]@features[, cluster.cols, drop=F],
            correlation=corr.list,
            graph=graph.list
        ))
    } else {
        return(dataset)
    }
}

#-------------------------------------------------------------------------------
#' Calculate the eigenfeature (first PC) for a cluster of features
#'
#' @description
#' Calculates eigenfeatures for a set of features defined by a column (see calculate_feature_clustering).
#' The eigenfeature is the first principal component of the features in the cluster.
#' This assay can be used for more fair testing of how well a featureset is associated to a trait.
#' 
#' @param dataset A \linkS4class{TglowDataset}
#' @param assay The assay on dataset to use
#' @param slot The assay slot to use ("data", "scale.data")
#' @param cluster.col Column in features slot containing cluster assignments to calculate eigenfeatures for
#' 
#' @details
#' NOTE: Eigenfeatures have arbitrary sign, so the direction of the eigenfeature may not be consistent across clusters.
#' NOTE: No scaling is applied prior to eigenfeature calculation, and it is assumed you use the scale-data slot.
#'
#' @return Returns the \linkS4class{TglowDataset} with a new assay named "<assay>_eigenfeatures"
#'  containing the eigenfeature for each cluster. The features slot contains the variance explained by the eigenfeature.
#' @export 
calculate_eigenfeatures <- function(dataset, assay, slot="scale.data", cluster.col) {
    check_dataset_assay_slot(dataset, assay, slot)
    
    data <- slot(dataset@assays[[assay]], slot)
    
    if (!cluster.col %in% colnames(dataset@assays[[assay]]@features)) {
        stop(paste0("cluster.col: ", cluster.col, " not found in features slot"))
    }
    
    clusters <- unique(dataset@assays[[assay]]@features[, cluster.col])
    
    if (any(is.na(clusters))) {
        stop(paste0("NA values found in cluster.col: ", cluster.col, ". All features must be assigned to a cluster."))
    }
    
    nclust <- length(clusters)
    
    result.mat <- matrix(NA, nrow = nrow(dataset@assays[[assay]]), ncol = nclust)
    colnames(result.mat) <- clusters
    rownames(result.mat) <- rownames(dataset@assays[[assay]])
    
    meta <- data.frame(id=clusters, var=NA, var_tot=NA, var_exp=NA)
    rownames(meta) <- clusters
    
    for (cluster in clusters) {
        selector <- dataset@assays[[assay]]@features[, cluster.col] == cluster
        meta[cluster, "n_feature"] <- sum(selector)
        meta[cluster, "features"]  <- paste0(rownames(dataset@assays[[assay]]@features)[selector], collapse=",")
        
        if (sum(selector) == 1) {
            result.mat[, cluster] <- data[, selector]
            meta[cluster, "var"] <- 1
            meta[cluster, "var_tot"] <- 1
            meta[cluster, "var_exp"] <- 1
            next
        }
        
        # Calculate the first principal component for the features in this cluster
        pca    <- irlba::prcomp_irlba(data[,selector,drop=F], n=1, center=T, scale=F)
    
        # Store the first principal component as the eigenfeature for this cluster
        result.mat[, cluster] <- pca$x[, 1]
        
        meta[cluster, "var"] <- pca$sdev[1]^2
        meta[cluster, "var_tot"] <- pca$totalvar
        
        if (pca$totalvar == 0) {
            meta[cluster, "var_exp"] <- NA
        } else {
            meta[cluster, "var_exp"] <- pca$sdev[1]^2 / pca$totalvar
        }
    }
    
    dataset@assays[[paste0(assay, "_eigenfeatures")]] <- TglowAssayFromMatrix(result.mat, features = meta) 
    
    return(dataset)
}

