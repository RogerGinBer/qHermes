#'@export
setGeneric("annotateChromPeaks", function(object, param) {
    standardGeneric("annotateChromPeaks")
})

setMethod("annotateChromPeaks",
          signature(object = "XCMSnExp", param = "ChemFormulaParam"),
          function(object, param){
    pks <- chromPeaks(object) %>% as.data.frame
    ionF <- param@ionFormulas
    ppm <- param@ppm
    anot <- lapply(pks$mz, function(mz){
        match <- RHermes:::binarySearch(matrix(ncol = 1, ionF$m), mz, ppm)
        paste(ionF$f[match], collapse = "_")
    })
    cpdata <- chromPeakData(object)
    cpdata$annotation <- anot
    chromPeakData(object) <- cpdata
    return(object)
})

