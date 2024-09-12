#-------------------------------------------------------------------------------
setGeneric("getImageData", function(object, j, slot = "data", drop = T) {
    standardGeneric("getImageData")
})

#-------------------------------------------------------------------------------
setGeneric("getImageDataByCell", function(object, j, slot = "data", drop = T) {
    standardGeneric("getImageDataByCell")
})

#-------------------------------------------------------------------------------
setGeneric("getImageDataAndFeatures", function(object, j, assay, slot = "data", drop = F) {
    standardGeneric("getImageDataAndFeatures")
})
