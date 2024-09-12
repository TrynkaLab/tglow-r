setClassUnion(name = "MatrixOrNull", members = c("matrix", "NULL"))

#-------------------------------------------------------------------------------
#' TglowAssay
#'
#' @slot data Slot to store unscaled HCI features. Rows are objects, columns are features
#' @slot scale.data Slot to store a scaled version of data
#' @slot features Data frame with feature level metadata
#'
setClass("TglowAssay",
    slots = list(
        data = "matrix",
        scale.data = "MatrixOrNull",
        features = "data.frame"
    ),
    prototype = list(scale.data = NULL)
)

#-------------------------------------------------------------------------------
#' TglowDataset
#'
#' Object to store HCI features along side object (cell) and image level metadata.
#'
#' @slot assays List of TglowAssays containing the data
#' @slot meta Data frame containing object (cell) level metadata
#' @slot object.ids Vector containing object (cell) level identifiers
#' @slot image.meta Image level metadata
#' @slot image.data TglowAssay object with image level features
#' @slot image.ids Vector of length objects with corresponding image ids
#' @slot reduction List to store PCA/UMAP embeddings
#' @slot graph Slot to store kNN/AkNN graph
#' @slot active.assay The default assay to analyze. Not all functions may respect this
#'
setClass("TglowDataset",
    slots = list(
        assays = "list",
        meta = "data.frame",
        object.ids = "character",
        image.meta = "data.frame",
        image.data = "TglowAssay",
        image.ids = "character",
        reduction = "list",
        graph = "ANY",
        active.assay = "character"
    ),
    prototype = list(
        reduction = list(),
        graph = NULL
    )
)
