#'@import MSnbase
#'@import RHermes
#'@export
findSOIpeaks <- function(MSnExp, DBfile, adfile){
    struct <- RHermesExp()
    struct <- setCluster(struct, BiocParallel::SerialParam())
    ppm <- struct@metadata@ExpParam@ppm

    #This would be what the end user should do:
    # struct <- setDB(struct, db = DBfile, adductfile = adfile)
    
    #For testing let's just use the "demo DB" defaults 
    #(just a few formulas and adducts)
    struct <- setDB(struct)
    
    prepro <- RHermes:::preprocessing(struct)
    IF_DB <- prepro[[1]]
    IC <- prepro[[2]]
    
    
    files <- fileNames(MSnExp)
    toAdd <- lapply(seq_along(files), function(i) {
        cur_MSnExp <- filterFile(MSnExp, i)
        imported <- extractToHermes(cur_MSnExp)
        ss <- RHermes:::OptScanSearch(DB = IF_DB[[1]],
                            raw = imported[[3]],
                            ppm = ppm,
                            labelled = FALSE,
                            IsoList = IC,
                            BiocParallelParam = struct@metadata@cluster)

        #Construction of S4 Object output
        RHermesPL(peaklist = ss, header = imported[[2]], raw = imported[[1]],
                    labelled = FALSE, filename = files[i])
    })
    struct@data@PL <- c(struct@data@PL, toAdd)
    struct@metadata@ExpParam@ionF <- IF_DB
    struct@metadata@ExpParam@isoList <- IC
    struct@metadata@filenames <- c(struct@metadata@filenames, files)
    
    struct <- findSOI(struct, getSOIpar(), fileID = seq_along(files))
    
    SOIs <- do.call("rbind", lapply(seq_along(files), function(soi){
        soi_list <- struct@data@SOI[[soi]]@SOIList
        soi_list$file <- files[[soi]]
        return(soi_list)
    }))
    return(SOIs)
}

#'@importFrom MSnbase extractSpectraData
#'@importFrom data.table rbindlist
extractToHermes <- function(MSnExp){
    df <- extractSpectraData(MSnExp)
    h <- as.data.frame(df[, 1:35])
    names(h)[names(h) == "rtime"] <- "retentionTime"
    
    #This is a bit slow, any ideas to speed it up?
    pks <- rbindlist(lapply(1:nrow(df), function(x){
        data.table(mz = df$mz[[x]],
                    rtiv = df$intensity[[x]],
                    rt = df$rtime[x])
    }))
    
    #To match with the expected format of RHermes: list(raw, header, filtered)
    return(list(pks, h, pks))
}
