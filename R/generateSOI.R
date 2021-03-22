#'@import MSnbase
#'@import RHermes
#'@export
findSOIpeaks <- function(MSnExp, DBfile = NA,
                         adfile = NA,
                         BPPARAM = bpparam(),
                         SOIParam = getSOIpar()){
    struct <- RHermesExp()
    struct <- setCluster(struct, BPPARAM)
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
    
    
    struct <- setCluster(struct, SerialParam())
    files <- fileNames(MSnExp)
    toAdd <- bplapply(seq_along(files), function(i, MSnExp, IF_DB, IC, struct, ppm, files) {
        cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
        imported <- XCHermes:::extractToHermes(cur_MSnExp)
        ss <- RHermes:::OptScanSearch(DB = IF_DB[[1]],
                            raw = imported[[3]],
                            ppm = ppm,
                            labelled = FALSE,
                            IsoList = IC,
                            BiocParallelParam = struct@metadata@cluster)

        #Construction of S4 Object output
        RHermes:::RHermesPL(peaklist = ss, header = imported[[2]], raw = imported[[1]],
                    labelled = FALSE, filename = files[i])
    }, BPPARAM = BPPARAM, MSnExp = MSnExp, IF_DB = IF_DB, IC = IC,
        struct = struct, ppm = ppm, files = files)
    struct@data@PL <- c(struct@data@PL, toAdd)
    struct@metadata@ExpParam@ionF <- IF_DB
    struct@metadata@ExpParam@isoList <- IC
    struct@metadata@filenames <- c(struct@metadata@filenames, files)
    
    struct <- setCluster(struct, BPPARAM)
    struct <- findSOI(struct, SOIParam, fileID = seq_along(files))
    SOIs <- do.call("rbind", lapply(seq_along(files), function(soi){
        soi_list <- struct@data@SOI[[soi]]@SOIList
        
        mzs <- lapply(soi_list$peaks, calculateMz)
        soi_list$mzmin <- sapply(mzs, function(x){x[[1]]})
        soi_list$mzmax <- sapply(mzs, function(x){x[[2]]})
        soi_list$mz <- sapply(mzs, function(x){x[[3]]})
        
        rts <- lapply(soi_list$peaks, calculateRt)
        soi_list$rtmin <- sapply(rts, function(x){x[[1]]})
        soi_list$rtmax <- sapply(rts, function(x){x[[2]]})
        soi_list$rt <- sapply(rts, function(x){x[[3]]})
        
        soi_list$into <- lapply(soi_list$peaks, function(x){sum(x[,2])})
        soi_list$intb <- soi_list$into
        soi_list$maxo <- soi_list$MaxInt

        soi_list$sample <- soi
        return(soi_list)
    })) %>% as.data.frame()
    
    MSnExp <- as(MSnExp, "XCMSnExp")
    if (nrow(SOIs) > 0) {
        names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
                    "intb", "maxo", "sample")
        chromPeaks(MSnExp) <- SOIs[, names] %>% as.matrix() %>%
            apply(., 2, as.numeric)
        chromPeakData(MSnExp) <- S4Vectors::DataFrame(SOIs[, 1:10],
                                                      ms_level = as.integer(1),
                                                      is_filled = FALSE)
    }
    return(MSnExp)
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


calculateMz <- function(peaks){
    mzmin <- peaks$mz[which.min(peaks$rt)]
    mzmax <- peaks$mz[which.max(peaks$rt)]
    mz <- peaks$mz[which.max(peaks$rtiv)]
    return(list(mzmin, mzmax, mz))
}

calculateRt <- function(peaks){
    rtmin <- min(peaks$rt)
    rtmax <- max(peaks$rt)
    rt <- peaks$rt[which.max(peaks$rtiv)]
    return(list(rtmin, rtmax, rt))
}
