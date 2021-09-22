#' @title fastSOI
#' @description A new implementation of SOI detection algorithm that avoids the
#'   datapoint abstraction.
#' @details The algorithm is faster than the regular SOI detection used in
#'   RHermes, but on the other hand loses the computational commodity of having
#'   all experimental datapoints annotated.
#' @param MSnExp An MSnExp object containing one or multiple MS1 files
#' @param ChemFormulaParam A ChemFormulaParam object, containing all ionic
#'   formulas to be annotated in the data.
#' @param minint Minimum intensity of the scans used in the calculated XICs.
#'   Defaults to 1000.
#' @importFrom xcms chromPeaks chromPeakData
#' @importFrom dplyr select rename
#' @importFrom BiocParallel SerialParam
#' @export
fastSOI <- function(MSnExp, ChemFormulaParam, minint = 1000){
    ppm <- ChemFormulaParam@ppm
    ionf <- ChemFormulaParam@ionFormulas
    message("Calculating fast SOIs")
    ionf$s <- ionf$m * (1 - ppm * 1e-6)
    ionf$e <- ionf$m * (1 + ppm * 1e-6)
    
    ##Processing each file
    files <- fileNames(MSnExp)
    SOIs <- bplapply(seq_along(files),
                      function(i) {
        cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
        pks <- qHermes:::extractToHermes(cur_MSnExp)
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
            sois <- qHermes:::soi_from_eic(eic$intensity)
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

## Internal function, similar to ROI detection in XCMS, but with a scan gap
## tolerance implemented.

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


#' @title fastSOIfromList
#' @description Quantify SOIs in multiple files using an RHermes SOI list as a
#'   template.
#' @details The function uses the fastSOI approach to quickly detect SOIs in
#'   multiple files using a target list of all SOIs in a given RHermesExp
#'   SOIList.
#' @param MSnExp
#' @param struct
#' @param SOI_id
#' @param rtwin
#' @param tol
#' @param thr
#'@export
fastSOIfromList <- function (MSnExp, struct, SOI_id = 1, rtwin = 3, tol = 3, 
                                thr = 1000) {
    ppm <- struct@metadata@ExpParam@ppm
    target_list <- RHermes::SOI(struct,SOI_id)@SOIList
    target_list$mmin <- target_list$mass * (1 - ppm * 1e-06)
    target_list$mmax <- target_list$mass * (1 + ppm * 1e-06)
    files <- fileNames(MSnExp)
    SOIs <- bplapply(seq_along(files), single_fastSOI_from_list,
                     MSnExp = MSnExp, target_list = target_list, ppm = ppm,
                     rtwin = rtwin, tol = tol, thr = thr,
                     BPPARAM = bpparam()
    )
    SOIs <- do.call("rbind", SOIs)
    row.names(SOIs) <- NULL
    MSnExp <- as(MSnExp, "XCMSnExp")
    if (nrow(SOIs) > 0) {
        names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", 
                   "into", "intb", "maxo", "sn", "sample",
                   "SOIidx", "nscans")
        SOIs$sn <- 1000
        cp_mat <- SOIs[, names] %>% as.matrix() %>% apply(., 2, as.numeric)
        row.names(cp_mat) <- seq(nrow(SOIs))
        chromPeaks(MSnExp) <- cp_mat
        chromPeakData(MSnExp) <- S4Vectors::DataFrame(
            SOIs[, c("rtmin", "peaks", "sample", "SOIidx", "formula")],
            ms_level = rep(as.integer(1), times=nrow(SOIs)),
            is_filled = rep(FALSE, times=nrow(SOIs)),
            row.names = seq(nrow(SOIs))
        )
    }
    return(MSnExp)
}

single_fastSOI_from_list <- function(i, MSnExp, target_list, ppm, rtwin, tol,
                                      thr){
    cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
    raw_data <- extractToHermes(cur_MSnExp)
    mzs <- unlist(raw_data[[3]][, 1])
    ints <- unlist(raw_data[[3]][, 2])
    h <- raw_data[[2]]
    scantime <- h$retentionTime
    valsPerSpect <- h$originalPeaksCount
    scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
    SL <- lapply(seq(nrow(target_list)), function(x) {
        tgt <- target_list[x, ]
        sr <- which(scantime >= (tgt$start - rtwin) &
                      scantime <= (tgt$end + rtwin))
        sr <- c(min(sr), max(sr))
        eic <- .Call("getEIC", mzs, ints, scanindex, 
                     as.double(as.numeric(c(tgt$mmin, tgt$mmax))), 
                     as.integer(sr),
                     as.integer(length(scanindex)), PACKAGE = "xcms")
        sois <- soi_from_eic(eic = eic$intensity, tol = tol)
        sl <- apply(sois, 1, function(x) x[2] - x[1])
        sois <- sois[which(sl > 0), , drop = FALSE]
        sl <- sl[which(sl > 0)]
        if (nrow(sois) == 0) {return()}
        sois <- as.data.frame(sois)
        sois$start <- scantime[sr[1]:sr[2]][sois$scstart]
        sois$end <- scantime[sr[1]:sr[2]][sois$scend]
        sois$rt <- sapply(seq(nrow(sois)),
                          function(j) median(scantime[sr[1]:sr[2]][sois$scstart[j]:sois$scend[j]]))
        sois$length <- sois$end - sois$start
        sois$formula <- tgt$formula
        sois$peaks <- apply(sois, 1, function(row) {
          data.frame(rt = scantime[row[1]:row[2]],
                        rtiv = eic$intensity[row[1]:row[2]])
        })
        sois$mass <- tgt$mass
        sois <- dplyr::select(sois, start, rt, end, length, formula, 
                              peaks, mass)
        sois$SOIidx <- rep(x, times = nrow(sois))
        sois$nscans <- sl
        return(sois)
    })
    if (is.null(SL)) {return()}
    SL <- do.call("rbind", SL)
    SL <- dplyr::rename(SL, rtmin = start, rtmax = end, mz = mass)
    SL$mzmin <- as.numeric(SL$mz) * (1 - ppm * 1e-06)
    SL$mzmax <- as.numeric(SL$mz) * (1 + ppm * 1e-06)
    SL$into <- sapply(SL$peaks, function(x) {sum(x[, 2])})
    SL$intb <- SL$into
    SL$maxo <- sapply(SL$peaks, function(x) {max(x[, 2])})
    SL$sample <- i
    return(SL)
}

#'@export
SOIfiltbyILv1 <- function(IL, SOIList, par){
    # SOIList as SOI@SOIList
    # IL as full IL() object
    ILanot <- IL@annotation
    filt <- unlist(sapply(seq_along(ILanot),function(i){
        a <- ILanot[[i]]
        r <- unlist(sapply(seq(nrow(a)),function(j){
            y <- a[j,]
            which(SOIList$start==y$start &
                      SOIList$end==y$end &
                      SOIList$mass==y$mass )
        }))
        # This can be avoided once OriginalSOI index is corrected
        if(length(r)>0){ 
            return(r)  
        }else{return()}
    }))
    soifilt <- SOIList[unique(filt),]
    return(soifilt)
}

SOIfiltbyILv2 <- function(IL,SOIList,par){
    # This function is for when OriginalSOI index is corrected
    # should add a flag if IL SOIid does not match 
    ILanot <- IL@annotation
    filt <- unlist(sapply(seq(length(ILanot)),function(i){
        a <- ILanot[[i]]
        r <- a$metadata[[1]]$originalSOI
        r <- unique(as.numeric(r))
            return(r)  
    }))
    soifilt <- SOIList[unique(filt),]
    return(soifilt)
    # nsoi <- length(RHermesExp@data@SOI)
    # RHermesExp@data@SOI[[nsoi+1]] <- RHermesExp@data@SOI[[nsoi]]
    # RHermesExp@data@SOI[[nsoi+1]]@SOIList <- RES
    # return(RHermesExp)
}

#'@export
consensusSOI <- function(RHermesExp,SOIids,minSOI=NULL,rtwin=5){
    if(length(SOIids)==1){stop("SOIids value too low")}
    if(is.null(minSOI)){stop("Please select a value for minSOI")}
    if(length(SOIids)<minSOI){stop("minSOI value too large, please select a valid value for minSOI")}
    slist <- lapply(SOIids,SOI,struct=RHermesExp)
    sdf <- lapply(seq(length(slist)),function(x){
        r <- slist[[x]]@SOIList
        r$soiN <- x
        return(r)
    })
    sapply(sdf,nrow)
    sdf <- do.call("rbind",sdf)
    RES <- lapply(seq(nrow(sdf)),function(x) return())
    exlist <- c()
    for(i in seq(nrow(sdf))){
        if(i %in% exlist){next}
        x <- sdf[i,]
        idx <- which(abs(sdf$start-x$start)<rtwin &
                         abs(sdf$end-x$end)<rtwin &
                         sdf$formula==x$formula)
        if(any(idx%in%exlist)) {idx <- idx[-which(idx%in%exlist)]}
        usoi <- sdf[idx,]
        if(length(unique(usoi$soiN))>=minSOI){
            newsoi <- usoi[which.max(usoi$MaxInt),]
            newsoi$start <- min(usoi$start) #mean?
            newsoi$end <- max(usoi$end) #mean?
            newsoi$length <- newsoi$end-newsoi$start
            # sum $peaks of all to unify?
            newsoi <- newsoi[,1:17]
            exlist <- c(exlist,idx)
            # too slow? potser eliminar de sdf directament?
            RES[[i]] <- newsoi
        }
        
    }
    RES <- do.call("rbind",RES)
    return(RES)
    # nsoi <- length(RHermesExp@data@SOI)
    # RHermesExp@data@SOI[[nsoi+1]] <- RHermesExp@data@SOI[[nsoi]]
    # RHermesExp@data@SOI[[nsoi+1]]@SOIList <- RES
    # return(RHermesExp)
}