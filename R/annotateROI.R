annotateROI <- function(object, DBfile = NA, adfile = NA, ppm = 3){
    struct <- RHermesExp()
    ppm <- struct@metadata@ExpParam@ppm

    #Selecting formulas and adducts
    if(!is.na(DBfile)){
        if(!is.na(adfile)){
            struct <- setDB(struct, db = "custom", filename = DBfile,
                            adductfile = adfile)
        } else {
            struct <- setDB(struct, db = "custom", filename = DBfile)
        }
    } else {
        struct <- setDB(struct)
    }
    
    prepro <- RHermes:::preprocessing(struct)
    IF_DB <- prepro[[1]]
    IC <- prepro[[2]]
    
    pks <- chromPeaks(object) %>% as.data.frame
    IF_DB[[1]] <- IF_DB[[1]][order(IF_DB[[1]]$m),]
    anot <- lapply(pks$mz, function(mz){
        match <- RHermes:::binarySearch(matrix(ncol=1, IF_DB[[1]]$m), mz, ppm)
        paste(IF_DB[[1]]$f[match], collapse = "_")
    })
    cpdata <- chromPeakData(object)
    cpdata$anot <- anot
    chromPeakData(object) <- cpdata
    
    return(object)
}
