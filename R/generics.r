#-------------------------------------------------------------------------------
#' Get and set object Ids
#'
#' @description Get and set object Id's on TglowDataset, TglowAssay and TglowReduction
#' @param object The object to get or set
#' @param value If setting using <- the new object ids to assign
#'
#' @returns The object with the new Ids
#'
#' @rdname objectIds
#' @export
setGeneric("objectIds", function(object) {
    standardGeneric("objectIds")
})

#' @rdname objectIds
#' @export
setGeneric("objectIds<-", function(object, value) {
    standardGeneric("objectIds<-")
})

#-------------------------------------------------------------------------------
#' Get and set image Ids
#'
#' @description Get and set image Id's on TglowDataset
#' @param object The object to get or set
#' @param value If setting using <- the new object ids to assign
#'
#' @returns The object with the new Ids
#'
#' @rdname imageIds
#' @export
setGeneric("imageIds", function(object) {
    standardGeneric("imageIds")
})

#' @rdname imageIds
#' @export
setGeneric("imageIds<-", function(object, value) {
    standardGeneric("imageIds<-")
})



#-------------------------------------------------------------------------------
#' Fetch image data or meta data from a tglow object
#'
#' @description Select columns from assay.image, image.meta and return them
#' as a data.frame per image
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param j character with column names from image.meta or assay.image to select
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetching image.data columns, leave at NULL
#' @param slot The assay slot to use ("data", "scale.data")
#' @param drop should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setGeneric("getImageData", function(dataset, j, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getImageData")
})


#-------------------------------------------------------------------------------
#' Fetch image data or meta data from a tglow object per object (cell)
#'
#' @description Select columns from assay.image, image.meta and return them
#' as a data.frame per object (cell)
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param j Character with column names from image.meta or assay.image to select
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetching image.data columns, leave at NULL
#' @param slot The assay slot to use ("data", "scale.data")
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
#-------------------------------------------------------------------------------
setGeneric("getImageDataByObject", function(dataset, j, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getImageDataByObject")
})

#-------------------------------------------------------------------------------
#' Get image data and features per object (cell)
#'
#' @description Select columns from assay, assay.image, image.meta from 'data' or 'scale.data' slots
#' and return them as a data.frame
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param j Character with column names from assay, assay.image, image.meta to select
#' @param assay The assay to select from. If not fetching assay columns, leave at NULL
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetching image.data columns, leave at NULL
#' @param slot The assay slot to use ("data", "scale.data")
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setGeneric("getDataByObject", function(dataset, j, assay = NULL, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getDataByObject")
})

#-------------------------------------------------------------------------------
#' Check if column names are available on a dataset
#'
#' @description Check whether column names are available on the dataset
#'
#' @param dataset A \linkS4class{TglowDataset}
#' @param j Character with column names from assay, assay.image, image.meta to select
#' @param assay The assay to select from. If not fetching assay columns, leave at NULL
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetching image.data columns, leave at NULL
#' @param slot The assay slot to use ("data", "scale.data")
#' @param drop Should cols be dropped or not
#' @param return.names If TRUE, return a character vector of available names instead of a logical
#' @returns A logical (or named logical vector if return.names=TRUE) indicating column availability
#' @export
setGeneric("isAvailable", function(dataset, j, assay, assay.image = NULL, slot, return.names = FALSE) {
    standardGeneric("isAvailable")
})


#-------------------------------------------------------------------------------
#' Check if a Tglow family of objects is valid
#'
#' @description Greedily check if a Tglow family of objects is valid.
#' Works on TglowDataset, TglowAssay, TglowReduction, TglowMatrix
#'
#' @details
#' Checks are peformed greedily, so as soon as an issue if found, FALSE is returned
#' and a warning for that issue is raised with more details. I.e. finding one issue
#' does not mean there are not other issues!
#'
#' For TglowDataset, object.names can be left at NULL, as it uses @object.ids slot to compare against.
#' If not NULL it will use the vector of names you supply to check
#'
#' TglowAssay, TglowReduction and TglowMatrix must all have object.names supplied, otherwise an error
#' is thrown.
#'
#' @param object A \linkS4class{TglowDataset}, \linkS4class{TglowAssay}, \linkS4class{TglowReduction} or \linkS4class{TglowMatrix}
#' @param object.names A character vector of object names to validate against. If NULL, uses object@object.ids
#' @returns A logical indicating validitiy
#' @export
setGeneric("isValid", function(object, object.names = NULL) {
    standardGeneric("isValid")
})
