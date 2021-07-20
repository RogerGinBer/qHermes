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
<<<<<<< Updated upstream
    pks <- extractToHermes(cur_MSnExp)
    mzs <- unlist(pks[[3]][, 1])
    ints <- unlist(pks[[3]][, 2])
    h <- pks[[2]]
=======
    raw_data <- qHermes:::extractToHermes(cur_MSnExp)
    mzs <- unlist(raw_data[[3]][, 1])
    ints <- unlist(raw_data[[3]][, 2])
    h <- raw_data[[2]]
>>>>>>> Stashed changes
    scantime <- h$retentionTime
    valsPerSpect <- h$originalPeaksCount
    scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
    SL <- lapply(seq(nrow(target_list)), function(x) {
<<<<<<< Updated upstream
      tgt <- target_list[x, ]
      sr <- which(scantime >= (tgt$start - rtwin) & scantime <= (tgt$end + rtwin))
      sr <- c(min(sr),max(sr))
      eic <- .Call("getEIC", mzs, ints, scanindex, 
                   as.double(as.numeric(c(tgt$mmin, tgt$mmax))), 
                   as.integer(sr),
                   as.integer(length(scanindex)), PACKAGE = "xcms")
      sois <- soi_from_eic(eic=eic$intensity,tol=tol)
      sl <- apply(sois, 1, function(x) x[2] - x[1])
      sois <- sois[which(sl > 0),,drop=F]
      sl <- sl[which(sl > 0)]
      # sois <- sois[which(sl > (tgt$nscans * 0.5)),,drop=F]
      if (nrow(sois) == 0) {return()}
      sois <- as.data.frame(sois)
      sois$start <- scantime[sr[1]:sr[2]][sois$scstart]
      sois$end <- scantime[sr[1]:sr[2]][sois$scend]
      sois$rt <- sapply(seq(nrow(sois)),function(j) median(scantime[sr[1]:sr[2]][sois$scstart[j]:sois$scend[j]]))
      sois$length <- sois$end - sois$start
      sois$formula <- tgt$formula
      sois$peaks <- apply(sois, 1, function(row) {
        data.frame(rt = scantime[row[1]:row[2]],
                   rtiv = eic$intensity[row[1]:row[2]])
      })
      sois$mass <- tgt$mass
      sois <- dplyr::select(sois, start, rt, end, length, formula, 
                            peaks, mass)
      sois$SOIidx <- rep(x,times=nrow(sois))
      sois$nscans <- sl
      return(sois)
=======
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
>>>>>>> Stashed changes
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


SOIfiltbyIL <- function(IL, SOIList, par){
    # SOIList as SOI@SOIList
    # IL as full IL() object
    ILList <- IL@IL
    ILanot <- IL@annotation
    # could this be simplified by do.call(rbind,Ilanot)?
    filt <- unlist(sapply(seq(nrow(ILList)),function(i){
        x <- ILList[i,]
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