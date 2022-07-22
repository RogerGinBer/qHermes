#### fastSOI-related ####

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
#' @importFrom dplyr select rename filter
#' @importFrom BiocParallel SerialParam bpparam bplapply
#' @importFrom MSnbase fileNames filterFile
#' @importFrom xcms chromPeaks chromPeakData
#' @importFrom S4Vectors DataFrame
#' @export
fastSOI <- function(MSnExp, ChemFormulaParam, minint = 1000,
                    peakwidth = c(5, 30)){
    ppm <- ChemFormulaParam@ppm
    ionf <- ChemFormulaParam@ionFormulas
    message("Calculating fast SOIs")
    ionf$s <- ionf$m * (1 - ppm * 1e-6)
    ionf$e <- ionf$m * (1 + ppm * 1e-6)
    
    ##Processing each file
    files <- fileNames(MSnExp)
    SOIs <- bplapply(seq_along(files), function(i) {
        message("Current file: ", files[i], " ", Sys.time())
        cur_MSnExp <- filterFile(MSnExp, i)
        pks <- extractToHermes(cur_MSnExp)
        mzs <- pks[[3]]$mz
        ints <- pks[[3]]$rtiv
        
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
            sois <- soi_from_eic(eic$intensity)
            if(nrow(sois) == 0){return()}
            sois <- as.data.frame(sois)
            
            sois$start <- scantime[sois$scstart]
            sois$end <- scantime[sois$scend]
            sois$length <- sois$end - sois$start
            sois$formula <- x[1]
            sois$peaks <- apply(sois, 1, function(row){
                x <- row[1]:row[2]
                l <- list(rt = scantime[x], rtiv = eic$intensity[x])
                class(l) <- "data.frame"
                l
            })
            sois$mass <- as.numeric(x[2]) 
            sois <- sois[, c("start", "end", "length", "formula", "peaks", "mass")]
            return(sois)
        })
        if(is.null(SL)){return()}
        SL <- do.call("rbind", SL)
        SL <- filter(SL, length >= peakwidth[1])
        suppressMessages({SL <- groupShort(SL, maxlen = peakwidth[2])})
        SL <- rename(SL, rtmin = start, rtmax = end, mz = mass)
        SL <- filter(SL, length >= peakwidth[1])
        SL$mzmin <- as.numeric(SL$mz) * (1 - ppm * 1e-6)
        SL$mzmax <- as.numeric(SL$mz) * (1 + ppm * 1e-6)
        SL$into <- sapply(SL$peaks, function(x){sum(x[,2])})
        SL$intb <- SL$into
        SL$maxo <- sapply(SL$peaks, function(x){max(x[,2])})
        SL$rt <- sapply(SL$peaks, function(x){x[which.max(x[,2]),1]})
        SL$sample <- i
        SL
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

## soi_from_eic: similar to ROI detection in XCMS, but with a scan gap
## tolerance implemented.

soi_from_eic <- function(eic, thr = 1000, tol = 5){
    above <- which(eic > thr)
    if (length(above) == 0) 
        return(matrix(nrow = 0, ncol = 2))
    m <- diff(above) < tol
    in_soi <- FALSE
    st <- c()
    end <- c()
    for (i in seq_along(m)) {
        if (m[i] & !in_soi) {
            st <- c(st, above[i])
            in_soi <- TRUE
        } else if (!m[i] & in_soi) {
            end <- c(end, above[i])
            in_soi <- FALSE
        }
    }
    if (in_soi)
        end <- c(end, above[i])
    if (length(st) == 0 | length(end) == 0)
        return(matrix(nrow = 0, ncol = 2))
    return(cbind(scstart = st, scend = end))
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


# Internally edefine groupShort so that it does not run centWave peak-picking
# (time consuming)

groupShort <- function(Groups, maxlen){
  message("Shortening and selecting long groups:")
  SG <- filter(Groups, length <= maxlen)
  LG <- filter(Groups, length > maxlen)
  LG <- lapply(seq_len(nrow(LG)), parallelGroupShort, LG, maxlen)
  LG <- do.call(rbind, LG)
  return(rbind(SG, LG))
}
parallelGroupShort <- function(i, LG, maxlen){
  ms1data <- LG$peaks[[i]]
  curGR <- LG[i,]
  
  #Divide long traces into equal-sized smaller traces
  times <- seq(from = curGR[1, 1][[1]], to = curGR[1, 2][[1]],
               length.out = ceiling(curGR[1, 3][[1]] / maxlen) + 1)
  deltat <- times[2] - times[1]
  NewGR <- data.table(start = times[-length(times)], end = times[-1],
                      length = deltat, formula = curGR$formula,
                      peaks = curGR$peaks, mass = curGR$mass)
  
  #Data point redistribution within each new SOI
  NewGR[, "peaks"] <- apply(NewGR, 1, function(x) {
    pks <- x[5][[1]]
    return(pks[between(pks$rt, x[1][[1]], x[2][[1]]), ])
  })
  return(NewGR)
}


#### fastSOIfromList-related ####

#' @title fastSOIfromList
#' @description Quantify SOIs in multiple files using an RHermes SOI list as a
#'   template.
#' @details The function uses the fastSOI approach to quickly detect SOIs in
#'   multiple files using a target list of all SOIs in a given RHermesExp
#'   SOIList.
#' @inheritParams fastSOI
#' @inheritParams annotateFeaturesFromList
#' @param tol Numeric. Gap (datapoints with <thr intensity) tolerance. 
#' @param thr Numeric. Intensity threshold for signal detection.
#'@export
fastSOIfromList <- function (MSnExp, RHermesExp, SOI_id = 1, RTtol = 3, tol = 3, 
                                thr = 1000) {
    ppm <- RHermesExp@metadata@ExpParam@ppm
    target_list <- RHermes::SOI(RHermesExp, SOI_id)@SOIList
    target_list$mmin <- target_list$mass * (1 - ppm * 1e-06)
    target_list$mmax <- target_list$mass * (1 + ppm * 1e-06)
    files <- fileNames(MSnExp)
    SOIs <- bplapply(seq_along(files), single_fastSOI_from_list,
                     MSnExp = MSnExp, target_list = target_list, ppm = ppm,
                     rtwin = RTtol, tol = tol, thr = thr, files = files,
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
            is_filled = rep(FALSE, times = nrow(SOIs)),
            row.names = seq(nrow(SOIs))
        )
    }
    return(MSnExp)
}

fastSOIfromDF <- function (MSnExp, target_list, ppm = 5, SOI_id = 1, rtwin = 3,
                           tol = 3, thr = 1000) {
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
               "into", "intb", "maxo", "sn", "sample")
    SOIs$sn <- 1000
    cp_mat <- SOIs[, names] %>% as.matrix() %>% apply(., 2, as.numeric)
    row.names(cp_mat) <- seq(nrow(SOIs))
    chromPeaks(MSnExp) <- cp_mat
    chromPeakData(MSnExp) <- S4Vectors::DataFrame(
        ms_level = rep(as.integer(1), times = nrow(SOIs)),
        is_filled = rep(TRUE, times = nrow(SOIs)),
        row.names = seq(nrow(SOIs))
    )
  }
  return(MSnExp)
}

single_fastSOI_from_list <- function(i, MSnExp, target_list, ppm, rtwin, tol, thr, files){
    message("Current file: ", files[i], " ", Sys.time())
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
        sois <- sois[, c("start", "rt", "end", "length", "formula", "peaks",
                         "mass")]
        sois$SOIidx <- rep(x, times = nrow(sois))
        sois$nscans <- sl
        return(sois)
    })
    if (is.null(SL)) {return()}
    SL <- do.call("rbind", SL)
    SL <- rename(SL, rtmin = start, rtmax = end, mz = mass)
    SL$mzmin <- as.numeric(SL$mz) * (1 - ppm * 1e-06)
    SL$mzmax <- as.numeric(SL$mz) * (1 + ppm * 1e-06)
    SL$into <- sapply(SL$peaks, function(x) {sum(x[, 2])})
    SL$intb <- SL$into
    SL$maxo <- sapply(SL$peaks, function(x) {max(x[, 2])})
    SL$sample <- i
    return(SL)
}

#' @title filterSOIByIL
#' @description Filter a SOI list and keep only those SOIs that were included in
#'   the inclusion list.
#' @details The rationale of this function is that one usually wants to quantify
#'   only those signals for which he has acquired MS2 data.
#' @inheritParams annotateFeaturesFromList
#' @param IL_id Numeric. Index of the inclusion list that will be used to filter the SOIs.

#' @importFrom RHermes SOI IL 
#'@export
filterSOIByIL <- function(RHermesExp, SOI_id, IL_id){
    ILanot <- IL(RHermesExp, IL_id)@annotation
    SOIList <- SOI(RHermesExp, SOI_id)@SOIList
    filt <- unlist(sapply(seq_along(ILanot), function(i) {
        inclusionGroup <- ILanot[[i]]
        SOI_number <- inclusionGroup$metadata[[1]]$originalSOI
        return(unique(as.numeric(SOI_number)))
    }))
    RHermesExp@data@SOI[[SOI_id]]@SOIList <- SOIList[unique(filt), ]
    return(RHermesExp)
}



# usage
# MsnExp  <- readMSData (XXXXX) - Raw data onDisk of qHermes files
# cwp <- CentWaveParam() - peakwidth always c(5,30)?
# rtp <-  ObiwarpParam() - Indicate centersample?
# #'@export
# getCorrectMatrix <- function(MsnExp, cwp, rtp){
#     MsnExp <- findChromPeaks(MsnExp, param = cwp)
#     MsnExp <- adjustRtime(MsnExp, param = rtp)
#     
#     rtRaw <- rtime(MsnExp, bySample = TRUE, adjusted = FALSE)
#     rtAdj <- rtime(MsnExp, bySample = TRUE, adjusted = TRUE)
#     rtr <- rtRaw[[1]] #raw data instrumeant scan rate
#     stepRt <- round(min(diff(rtr)), 1) #+0.1
#     rtmax <- ceiling(max(unlist(rtRaw)))
#     tpoints <- seq(0, rtmax, stepRt)
#     xpfiles <- fileNames(MsnExp)
#     
#     rtCorr <- lapply(seq(xpfiles), function(i) {
#         rta <- rtAdj[[i]]
#         rtr <- rtRaw[[i]]
#         r <- sapply(tpoints, function(j) {
#             minj <- j - stepRt
#             idx <- which(rtr >= minj & rtr < j)
#             if (length(idx) == 0) {
#                 return(c(j, 0))
#             } else {
#                 rtdif <- median(rta[idx] - rtr[idx])
#                 return(c(j, rtdif))
#                 # rtdif is the value to be subtracted
#                 # to raw RT to match corrected RT
#                 # rtAdj = rtRaw + rtdif
#             }
#         })
#         r <- t(r)
#         r
#     })
#     ## there is still some gaps
#     ## need to do some imputation before application to qHermes
#     idx <- seq(rtCorr)
#     idx <- idx[-rtp@centerSample]
#     for (i in idx) {
#         rtc <- rtCorr[[i]]
#         rt0 <- which(rtc[,2] == 0)
#         rt0 <- rt0[which(rt0 > 1 & rt0 < nrow(rtc))]
#         for (j in rt0) {
#             a <- j - 1
#             if (a %in% rt0) {
#                 while(a %in% rt0){
#                     a <- a - 1
#                     if (a == 0) {next}
#                 }
#             }
#             b <- j + 1
#             if (b %in% rt0) {
#                 while (b %in% rt0) {
#                     b <- b + 1
#                     if (b == 0) {
#                         next
#                     }
#                 }
#             }
#             a <- rtc[a, ]
#             b <- rtc[b, ]
#             #impute
#             D1 <- sqrt(((b[2] - a[2]) ^ 2) + ((b[1] - a[1]) ^ 2))
#             rtc[j, 2] <- (a[2] + (((b[1] - a[1]) / D1) * (b[2] - a[2])))
#         }
#         rtCorr[[i]] <- rtc
#     }
#     return(rtCorr)
# }

# #' @export
# SOIFeatureDefinitions <- function(object){
#     def <- featureDefinitions(object) %>% as.data.frame()
#     pk_data <- chromPeakData(object) %>% as.data.frame()
#     def$formula <-
#         lapply(def$peakidx, function(id) {
#             pk_data$formula[id[1]]
#         })
#     def$anot <- lapply(def$peakidx, function(id) {
#         pk_data$anot[id[1]]
#     })
#     def <-
#         cbind(def, featureValues(object, value = "maxo") %>% as.data.frame())
#     return(def)
# }
