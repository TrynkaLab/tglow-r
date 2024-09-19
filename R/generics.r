#-------------------------------------------------------------------------------
setGeneric("getImageData", function(object, j, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getImageData")
})

#-------------------------------------------------------------------------------
setGeneric("getImageDataByObject", function(object, j, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getImageDataByObject")
})

#-------------------------------------------------------------------------------
setGeneric("getDataByObject", function(object, j, assay = NULL, assay.image = NULL, slot = "data", drop = TRUE) {
    standardGeneric("getDataByObject")
})
