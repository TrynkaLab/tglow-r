#-------------------------------------------------------------------------------
#' Get and set object Ids
#'
#' @description Get and set object Id's on TglowDataset, TglowAssay and TglowReduction
#' @param object The object to get or set
#' @param value If setting using <- the new object ids to assing
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
#' Fetch image data or meta data from a tglow object
#'
#' @description Select columns from assay.image, image.meta and return them
#' as a data.frame per image
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j character with column names from image.meta or assay.image to select
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetcing image.data columns, leave at NULL
#' @param slot slot to fetch features fromdata or scale.data
#' @param drop should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setGeneric("getImageData", function(object, j, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getImageData")
})


#-------------------------------------------------------------------------------
#' Fetch image data or meta data from a tglow object per object (cell)
#'
#' @description Select columns from assay.image, image.meta and return them
#' as a data.frame per object (cell)
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j Character with column names from image.meta or assay.image to select
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetcing image.data columns, leave at NULL
#' @param slot Slot to fetch features from data or scale.data
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
#-------------------------------------------------------------------------------
setGeneric("getImageDataByObject", function(object, j, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getImageDataByObject")
})

#-------------------------------------------------------------------------------
#' Get image data and features per object (cell)
#'
#' @description Select columns from assay, assay.image, image.meta from 'data' or 'scale.data' slots
#' and return them as a data.frame
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j Character with column names from assay, assay.image, image.meta to select
#' @param assay The assay to select from. If not fetcing assay columns, leave at NULL
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetcing image.data columns, leave at NULL
#' @param slot Slot to fetch features from: "data" or "scale.data"
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setGeneric("getDataByObject", function(object, j, assay = NULL, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getDataByObject")
})

#-------------------------------------------------------------------------------
#' Check if column names are available on a dataset
#'
#' @description Select columns from assay, assay.image, image.meta from 'data' or 'scale.data' slots
#' and return them as a data.frame
#'
#' @param object A \linkS4class{TglowDataset}
#' @param j Character with column names from assay, assay.image, image.meta to select
#' @param assay The assay to select from. If not fetcing assay columns, leave at NULL
#' @param assay.image Which image assay to use, "image.data", "image.data.trans" or "image.data.norm". If not fetcing image.data columns, leave at NULL
#' @param slot Slot to fetch features from: "data" or "scale.data"
#' @param drop Should cols be dropped or not
#' @returns A data frame with the corresponding columns
#' @export
setGeneric("isAvailable", function(object, j, assay, assay.image = NULL, slot, return.names = FALSE) {
    standardGeneric("isAvailable")
})
