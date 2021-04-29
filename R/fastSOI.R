soi_from_eic <- function(eic, thr = 1000, tol = 5){
    above <- which(eic > thr)
    if(length(above) == 0){return(matrix(nrow=0,ncol=2))}
    m <- diff(above) < tol
    in_soi <- F
    st <- c()
    end <- c()
    for(i in seq_along(m)){
        if(m[i] & !in_soi){
            st <- c(st, above[i])
            in_soi <- T
        } else if(!m[i] & in_soi){
            end <- c(end, above[i])
            in_soi <- F
        }
    }
    if(in_soi){end <- c(end, above[i])}
    if(length(st) == 0 | length(end) == 0){return(matrix(nrow=0,ncol=2))}
    return(cbind(scstart = st, scend = end))
}

#'@importFrom xcms chromPeaks chromPeakData
#'@importFrom dplyr select rename
#'@importFrom BiocParallel SerialParam
#'@export
fastSOI <- function(MSnExp, minint = 1000, ChemFormulaParam){
    ppm <- ChemFormulaParam@ppm
    ionf <- ChemFormulaParam@ionFormulas
    message("Calculating fast SOIs")
    ionf$s <- ionf$m*(1-ppm*1e-6)
    ionf$e <- ionf$m*(1+ppm*1e-6)
    
    ##Processing each file
    files <- fileNames(MSnExp)
    SOIs <- bplapply(seq_along(files),
                      function(i) {
        cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
        pks <- XCHermes:::extractToHermes(cur_MSnExp)
        mzs <- unlist(pks[[3]][,1])
        ints <- unlist(pks[[3]][,2])

        ## Define the values per spectrum:
        h <- pks[[2]]
        valsPerSpect <- h$originalPeaksCount
        scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
        scantime <- h$retentionTime
        
        ## Calculate SOIs from ionic formula list
        SL <- apply(ionf, 1, function(x){
            eic <- .Call("getEIC", mzs, ints, scanindex,
                            as.double(as.numeric(c(x[5], x[6]))),
                            as.integer(c(1,length(scanindex))),
                            as.integer(length(scanindex)),
                            PACKAGE = "xcms")
            if(all(eic$intensity < minint)){return()}
            sois <- XCHermes:::soi_from_eic(eic$intensity)
            if(nrow(sois) == 0){return()}
            sois <- as.data.frame(sois)
            
            sois$start <- scantime[sois$scstart]
            sois$end <- scantime[sois$scend]
            sois$length <- sois$end - sois$start
            sois$formula <- x[1]
            sois$peaks <- apply(sois, 1, function(row){
                data.frame(rt = scantime[row[1]:row[2]],
                           rtiv = eic$intensity[row[1]:row[2]])
            })
            sois$mass <- as.numeric(x[2]) 
            sois <- select(sois, start, end, length, formula, peaks, mass)
            return(sois)
        })
        if(is.null(SL)){return()}
        SL <- do.call("rbind", SL)
        
        suppressMessages({
            SL <- RHermes:::groupShort(SL, maxlen = 30,
                                    BPPARAM = BiocParallel::SerialParam())
        })
        SL <- rename(SL, rtmin = start, rtmax = end, mz = mass)
        
        SL$mzmin <- as.numeric(SL$mz) * (1 - ppm * 1e-6)
        SL$mzmax <- as.numeric(SL$mz) * (1 + ppm * 1e-6)
        SL$into <- sapply(SL$peaks, function(x){sum(x[,2])})
        SL$intb <- SL$into
        SL$maxo <- sapply(SL$peaks, function(x){max(x[,2])})
        SL$rt <- sapply(SL$peaks, function(x){x[which.max(x[,2]),1]})
        SL <- SL[SL$length > 5,] #Filter by min time
        SL$sample <- i
        
        ## Do some filtering here
        
        return(SL)
    },  BPPARAM = bpparam())
    SOIs <- do.call("rbind", SOIs)
    row.names(SOIs) <- NULL
    
    MSnExp <- as(MSnExp, "XCMSnExp")
    if (nrow(SOIs) > 0) {
        names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
                    "intb", "maxo", "sn", "sample")
        SOIs$sn <- 1e3
        cp_mat <- SOIs[, names] %>% as.matrix() %>% apply(., 2, as.numeric)
        row.names(cp_mat) <- seq(nrow(SOIs))
        chromPeaks(MSnExp) <- cp_mat
        chromPeakData(MSnExp) <- S4Vectors::DataFrame(SOIs[, 3:11],
                                                    ms_level = as.integer(1),
                                                    is_filled = FALSE,
                                                    row.names = seq(nrow(SOIs)))
    }
    return(MSnExp)
}

#'@export
filterSOIFromTemplate <- function(XCMSnExp, RHermesExp, SOI_id){
    cp <- chromPeaks(XCMSnExp)
    cpdata <- chromPeakData(XCMSnExp)
    SL <- SOI(RHermesExp, SOI_id)@SOIList
    unique_ions <- unique(SL$formula)
    toRemove <- rep(TRUE, nrow(cpdata))
    FALSE -> toRemove[sapply(cpdata$formula, function(x){x %in% unique_ions})]
    for(i in which(!toRemove)){
        matching_sois <- filter(SL, formula == cpdata$formula[[i]])
        good <- (cp[i, "rtmin"] > matching_sois$start - 20) &
                (cp[i, "rtmax"] < matching_sois$end + 20)
        if(!any(good)) toRemove[i] <- TRUE
    }
    chromPeaks(XCMSnExp) <- cp[!toRemove,]
    chromPeakData(XCMSnExp) <- cpdata[!toRemove,]
    return(XCMSnExp)
}

#'@export
plotFeature <- function(XCMSnExp, feature){
    cp <- chromPeaks(XCMSnExp)
    cpdata <- chromPeakData(XCMSnExp)%>% as.data.frame
    ft <- featureDefinitions(XCMSnExp)[feature,] %>% as.data.frame
    pks <- cpdata[ft$peakidx[[1]], ] 
    pknames <- fileNames(XCMSnExp)[cp[ft$peakidx[[1]], "sample"]] %>%
        basename
    plot_data <- lapply(1:nrow(pks), function(x){
        pk <- pks$peaks[[x]]
        pk <- pk[pk$rtiv != 0,]
        pk$samp <- pknames[x]
        pk <- rbind(data.frame(rt = min(pk$rt), rtiv = 0, samp = pknames[x]),
                    pk,
                    data.frame(rt = max(pk$rt), rtiv = 0, samp = pknames[x]))
    })
    plot_data <- do.call(rbind, plot_data)
    ggplot(plot_data) + geom_polygon(aes(x=rt, y=rtiv, fill = samp),
                                     alpha = 0.5)
} 
