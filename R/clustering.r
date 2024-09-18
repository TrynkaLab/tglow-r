#-------------------------------------------------------------------------------
#' Cluster using a a K nearest neighbor graph in the PCA space
#'
#' @description
#' Cluster using a a K nearest neighbor graph in the PCA space
#'
#'
#' @param dataset A tglow dataset
#' @param reduction The reduction to use for calculating UMAPs. If NULL re-calculated
#' @param pc.n How many PC's to use
#' @param k How many NN to calculate
#' @param method Clustering method to use 'louvain' or 'leiden'
#' @param resolution Resolution for method 'louvain' or 'leiden'
#' @param exact.nn Instead of using Seurat ANNOY for kNN use \code{\link[=nn2]{RANN::nn2()}}
#' @param col.out Column name to store the clustering under added to the meta slot on the output object
#'
#' @details
#' Here I use the Seurat implementation of the kNN, which is NOT and exact
#' kNN graph, but it is very fast. Practically it should perform pretty well.
#' The parameter k controls how many nearest neighbors to find.
#' The distance matrix is built from the prinicpal components
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
apply_clustering <- function(dataset, reduction, pc.n = NULL, k = 10, method = "louvain", resolution = 0.1, exact.nn = FALSE, col.out="clusters") {
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

    if (method == "louvain") {
        cl <- cluster_louvain(graph, resolution = resolution)
    } else if (method == "leiden") {
        cl <- cluster_leiden(graph, resolution = resolution)
    } else {
        stop("No valid cluster method")
    }

    # Make a vector of the cluster memberships for each cell
    dataset@meta[,col.out] <- membership(cl)

    return(dataset)
}
