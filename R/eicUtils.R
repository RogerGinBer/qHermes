extractEICDataFromPeaks <- function(XCMSnExp, PeakMatrix){
    # if(is(XCMSnExp, "XCMSnExp")) XCMSnExp <- as(XCMSnExp, "MSnExp")
    
    all_eics <- lapply(unique(PeakMatrix$sample), function(i){
        curPeakMatrix <- filter(PeakMatrix, sample == i)
        cur_MSnExp <- xcms::filterFile(XCMSnExp, i)
        raw_data <- qHermes:::extractToHermes(cur_MSnExp)
        mzs <- unlist(raw_data[[3]][, 1])
        ints <- unlist(raw_data[[3]][, 2])
        h <- raw_data[[2]]
        scantime <- h$retentionTime
        valsPerSpect <- h$originalPeaksCount
        scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
        
        apply(curPeakMatrix, 1, function(peak){
            peak <- as.numeric(peak[1:12])
            sr <- which(scantime >= (peak[5]) & scantime <= (peak[6]))
            sr <- c(min(sr), max(sr))
            if(any(is.infinite(sr))) return(data.frame())
            eic <- .Call("getEIC", mzs, ints, scanindex, 
                         as.double(as.numeric(c(peak[2], peak[3]))), 
                         as.integer(sr),
                         as.integer(length(scanindex)), PACKAGE = "xcms")
            eic <- as.data.frame(eic)
            eic$rt <- scantime[eic$scan]
            return(eic)
        })
    })
    all_eics <- unlist(all_eics, recursive = FALSE)
    PeakMatrix$eics <- all_eics
    return(PeakMatrix)
}

eic_sim <- function(eic1, eic2){
    eic1 <- bin_eic(eic1)
    eic2 <- bin_eic(eic2)
    RHermes:::cosineSim(eic1,eic2)
}

bin_eic <- function(eic, bwidth = 1){
    bin_breaks <- seq(floor(min(eic$rt)) - 5,
                      ceiling(max(eic$rt)) + 5,
                      bwidth)
    eic$bins <- cut(eic$rt, bin_breaks, labels = FALSE)
    binned_eic <- eic %>%
        group_by(bins) %>%
        summarise(rtiv = sum(intensity)) %>% 
        mutate(rt = bin_breaks[bins])
    return(binned_eic[,c("rt", "rtiv")])
}
