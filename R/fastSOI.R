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

#'@export
fastSOIfromList <- function (MSnExp, struct, SOI_id =1, rtwin=10) 
{
  ppm <- struct@metadata@ExpParam@ppm
  target_list <- RHermes::SOI(struct,SOI_id)@SOIList
  target_list$mmin <- target_list$mass * (1 - ppm * 1e-06)
  target_list$mmax <- target_list$mass * (1 + ppm * 1e-06)
  files <- fileNames(MSnExp)
  SOIs <- bplapply(seq_along(files), function(i,MSnExp,target_list,ppm) {
    cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
    pks <- extractToHermes(cur_MSnExp)
    mzs <- unlist(pks[[3]][, 1])
    ints <- unlist(pks[[3]][, 2])
    h <- pks[[2]]
    scantime <- h$retentionTime
    valsPerSpect <- h$originalPeaksCount
    scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
    SL <- lapply(seq(nrow(target_list)), function(x) {
      x <- target_list[x, ]
      sr <- which(scantime >= (x$start - rtwin) & scantime <= (x$end + rtwin))
      sr <- c(min(sr),max(sr))
      eic <- .Call("getEIC", mzs, ints, scanindex, 
                   as.double(as.numeric(c(x$mmin, x$mmax))), 
                   as.integer(sr),
                   as.integer(length(scanindex)), PACKAGE = "xcms")
      sois <- XCHermes:::soi_from_eic(eic$intensity)
      sl <- apply(sois, 1, function(x) x[2] - x[1])
      sois <- sois[which(sl > (x$nscans * 0.5)),,drop=F]
      if (nrow(sois) == 0) {return()}
      # Filters that could be possible applied or deleted (Jordi)
      # smaxint <- apply(sois, 1, function(x) max(eic$intensity[x[1]:x[2]]))
      # sois <- sois[which(smaxint > 1e4),,drop=F]
      # if (nrow(sois) == 0) {return()}
      # if (nrow(sois) >1) {
      #   xpeaks <- x$peaks[[1]]
      #   # xpeaks$rt <- floor(xpeaks$rt-min(xpeaks$rt)+1)
      #   xpeaks$rt <- 1:nrow(xpeaks)
      #   cossim <- apply(sois,1,function(s){
      #     eicpeaks <- data.table::data.table("rt"=1:(s[2]-s[1]+1),
      #                                        "rtiv"=eic$intensity[s[1]:s[2]])
      #     # eicpeaks$rt <- floor(eicpeaks$rt-min(eicpeaks$rt)+1)
      #     # plot(eicpeaks,type="l", main="EIC")
      #     # plot(xpeaks,type="l", main="SOIpeaks")
      #     RHermes:::cosineSim(xpeaks,eicpeaks)
      #   })
      #   # soiSIM per quedarme la soi_eic més semblant?
      #   sois <- sois[which(cossim > 0.5),, drop = FALSE]
      #   if (nrow(sois) == 0) {return()}
      # }
      sois <- as.data.frame(sois)
      sois$start <- scantime[sr[1]:sr[2]][sois$scstart]
      sois$end <- scantime[sr[1]:sr[2]][sois$scend]
      sois$rt <- sapply(seq(nrow(sois)),function(j) median(scantime[sr[1]:sr[2]][sois$scstart[j]:sois$scend[j]]))
      sois$length <- sois$end - sois$start
      sois$formula <- x$formula
      sois$peaks <- apply(sois, 1, function(row) {
        data.frame(rt = scantime[row[1]:row[2]],
                   rtiv = eic$intensity[row[1]:row[2]])
      })
      sois$mass <- x$mass
      sois <- dplyr::select(sois, start, rt, end, length, formula, 
                            peaks, mass)
      return(sois)
    })
    if (is.null(SL)) {
      return()
    }
    SL <- do.call("rbind", SL)
    # suppressMessages({
    #   SL <- RHermes:::groupShort(SL, maxlen = 30, BPPARAM = BiocParallel::SerialParam())
    # })
    SL <- dplyr::rename(SL, rtmin = start, rtmax = end, mz = mass)
    SL$mzmin <- as.numeric(SL$mz) * (1 - ppm * 1e-06)
    SL$mzmax <- as.numeric(SL$mz) * (1 + ppm * 1e-06)
    SL$into <- sapply(SL$peaks, function(x) {
      sum(x[, 2])
    })
    SL$intb <- SL$into
    SL$maxo <- sapply(SL$peaks, function(x) {
      max(x[, 2])
    })
    SL <- SL[SL$length > 5, ,drop=F]
    SL$sample <- i
    return(SL)
  },
  MSnExp = MSnExp, target_list = target_list, ppm = ppm,
  BPPARAM = bpparam()
  )
  SOIs <- do.call("rbind", SOIs)
  row.names(SOIs) <- NULL
  MSnExp <- as(MSnExp, "XCMSnExp")
  if (nrow(SOIs) > 0) {
    names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", 
               "into", "intb", "maxo", "sn", "sample")
    SOIs$sn <- 1000
    cp_mat <- SOIs[, names] %>% as.matrix() %>% apply(., 2, as.numeric)
    row.names(cp_mat) <- seq(nrow(SOIs))
    chromPeaks(MSnExp) <- cp_mat
    chromPeakData(MSnExp) <- S4Vectors::DataFrame(
      SOIs[,3:11],
      ms_level = as.integer(1),
      is_filled = FALSE,
      row.names = seq(nrow(SOIs))
    )
  }
  return(MSnExp)
}

