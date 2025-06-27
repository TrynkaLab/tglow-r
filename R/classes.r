setClassUnion(name = "MatrixOrNull", members = c("matrix", "NULL"))
setClassUnion(name = "NumericOrNull", members = c("numeric", "NULL"))
setClassUnion(name = "CharacterOrNull", members = c("character", "NULL"))

#-------------------------------------------------------------------------------
#' TglowMatrix
#'
#' Wrapper arround matrix to enable the use of $
setClass("TglowMatrix", contains = "matrix")
setClassUnion(name = "TglowMatrixOrNull", members = c("matrix", "NULL"))

#-------------------------------------------------------------------------------
#' TglowAssay
#'
#' @slot data Slot to store unscaled HCI features. Rows are objects, columns are features
#' @slot scale.data Slot to store a scaled version of data
#' @slot features Data frame with feature level metadata
#'
#' @exportMethod colnames
#' @exportMethod rownames
#' @exportMethod ncol
#' @exportMethod nrow
setClass("TglowAssay",
    slots = list(
        data = "TglowMatrix",
        scale.data = "TglowMatrixOrNull",
        features = "data.frame"
    ),
    prototype = list(scale.data = NULL)
)
setClassUnion(name = "TglowAssayOrNull", members = c("TglowAssay", "NULL"))


#-------------------------------------------------------------------------------
#' TglowDataset
#'
#' @description
#' Object to store HCI features along side object (cell) and image level metadata
#' Matrices are assumed to have row and column names set in the constructor function
#' This is important for downstream functionality
#'
#' @slot assays List of TglowAssays containing the data
#' @slot meta Data frame containing object (cell) level metadata
#' @slot object.ids Vector containing object (cell) level identifiers
#' @slot image.meta Image level metadata
#' @slot image.data TglowAssay object with image level features
#' @slot image.data.trans TglowAssay object with image level features after transformation
#' @slot image.data.norm TglowAssay object with image level features after transformation and normalization
#' @slot image.ids Vector of length objects with corresponding image ids
#' @slot reduction List to store PCA/UMAP embeddings
#' @slot graph Slot to store kNN/AkNN graph
#' @slot active.assay The default assay to analyze. Not all functions may respect this
#'
#' @exportMethod nrow
setClass("TglowDataset",
    slots = list(
        assays = "list",
        meta = "data.frame",
        object.ids = "character",
        image.meta = "data.frame",
        image.data = "TglowAssay",
        image.data.trans = "TglowAssayOrNull",
        image.data.norm = "TglowAssayOrNull",
        image.ids = "character",
        reduction = "list",
        graph = "ANY",
        active.assay = "character",
        feature.map="TglowFeatureMapOrNull"
    ),
    prototype = list(
        reduction = list(),
        graph = NULL,
        image.data.trans = NULL,
        image.data.norm = NULL,
        feature.map=NULL
    )
)

#-------------------------------------------------------------------------------
#' TglowFilter
#'
#' @description An object to describe a filter for either objects (cells) or
#' features
#'
#' @slot name Filter name
#' @slot column_pattern Collumn pattern to apply filters on. "all" represents all
#' @slot func Name of filter function
#' @slot threshold Filter threshold value
#' @slot active Logical indicating if the filter should be applied
#' @slot transpose Should data matrix be transposed prior to running func
setClass("TglowFilter",
    slots = list(
        name = "character",
        column_pattern = "character",
        func = "character",
        threshold = "numeric",
        transpose = "logical",
        active = "logical"
    ),
    prototype = list(
        transpose = FALSE,
        active = TRUE
    )
)


#-------------------------------------------------------------------------------
#' TglowReduction
#'
#' @slot x Slot to store PCA / UMAP coordinates
#' @slot sdev optional standard deviations
#' @slot sdev_total optional sum of sdev
#' @slot object Optional PCA / UMAP object
#'
setClass("TglowReduction",
    slots = list(
        x = "matrix",
        var = "NumericOrNull",
        var_total = "NumericOrNull",
        object = "ANY"
    ),
    prototype = list(var = NULL, var_total = NULL, object = NULL)
)


#-------------------------------------------------------------------------------
#' TglowFeatureMap
#'
#' @slot feature.x Feature describing object X
#' @slot feature.y Feature describing object Y
#' @slot feature.z Feature describing object Z
#' @slot feature.well Feature describing well position
#' @slot feature.plate Feature describing plate
#' @slot feature.field Feature describing field
setClass("TglowFeatureMap",
    slots = list(
        x = "TglowFeatureLocation",
        y = "TglowFeatureLocation",
        z = "TglowFeatureLocation",
        well = "TglowFeatureLocation",
        plate = "TglowFeatureLocation",
        field = "TglowFeatureLocation"
    )
)
setClassUnion(name = "TglowFeatureMapOrNull", members = c("TglowFeatureMap", "NULL"))


#-------------------------------------------------------------------------------
#' TglowFeatureMap
#'
#' @slot feature Character with the feature position
#' @slot assay to grab feature from, or NULL if metadata feature
#' @slot slot to grab feature from, or NULL if metadata feature
setClass("TglowFeatureLocation",
    slots = list (
        feature = "character",
        assay = "CharacterOrNull",
        slot = "CharacterOrNull"
    )
)
setClassUnion(name = "TglowFeatureLocationOrNull", members = c("TglowFeatureLocation", "NULL"))
