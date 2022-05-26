
findROIpeaks <- function(MSnExp, ppm = 3){
    files <- fileNames(MSnExp)
    ROIs <- lapply(seq_along(files), function(i){
        cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
        data <- qHermes:::extractToHermes(cur_MSnExp)
        pks <- data[[1]]
        h <- data[[2]]
        
        ## Get m/z and intensity values
        mzs <- pks[[1]]
        ints <- pks[[2]]
        
        ## Define the values per spectrum:
        valsPerSpect <- h$originalPeaksCount
        scant <- h$retentionTime
        
        roi <- centWave_orig(mzs, ints, scant, valsPerSpect, ppm = ppm, 
                                prefilter = c(3, 100),
                                peakwidth = c(10, 60),
                                snthresh = 6,
                                noise = 0)
        roi <- do.call(rbind, roi)

        roi <- apply(roi, 2, as.numeric)
        roi <- as.data.frame(roi)
        
        #Johannes, is this the correct way to recover the RT values from scmin and scmax?
        roi <-  mutate(roi, rtmin = scant[scmin]) %>%
                mutate(., rtmax = scant[scmax]) %>%
                mutate(., rt = (rtmax + rtmin)/2) %>% 
                rename(maxo = intensity) %>% mutate(into = maxo, intb = maxo)
        roi$sample <- i
        roi$sn <- 1e3
        return(roi)
    })
    ROIs <- do.call("rbind", ROIs)
    row.names(ROIs) <- NULL
    
    MSnExp <- as(MSnExp, "XCMSnExp")
    if (nrow(ROIs) > 0) {
        names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
                   "intb", "maxo", "sn", "sample")
        ROIs$sn <- 1e3
        cp_mat <- ROIs[, names] %>% as.matrix() %>% apply(., 2, as.numeric)
        row.names(cp_mat) <- seq(nrow(ROIs))
        chromPeaks(MSnExp) <- cp_mat
        chromPeakData(MSnExp) <- S4Vectors::DataFrame(
            ms_level = rep(as.integer(1), nrow(ROIs)),
            is_filled = rep(FALSE, nrow(ROIs)),
            row.names = seq(nrow(ROIs))
        )
    }
    return(MSnExp)
}




#Modified version of the internal Centwave detection that only generates ROIs
centWave_orig <- function(mz, int, scantime, valsPerSpect,
                           ppm = 25, peakwidth = c(20,50), snthresh = 10,
                           prefilter = c(3,100),
                           integrate = 1, fitgauss = FALSE,
                           noise = 0, ## noise.local=TRUE,
                           sleep = 0, verboseColumns = FALSE, roiList = list(),
                           firstBaselineCheck = TRUE, roiScales = NULL,
                           extendLengthMSW = FALSE) {
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to match. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect) ## Get index vector for C calls
    if (!is.double(mz))
        mz <- as.double(mz)
    if (!is.double(int))
        int <- as.double(int)
    if (!is.logical(firstBaselineCheck))
        stop("Parameter 'firstBaselineCheck' should be logical!")
    if (length(firstBaselineCheck) != 1)
        stop("Parameter 'firstBaselineCheck' should be a single logical !")
    if (length(roiScales) > 0)
        if (length(roiScales) != length(roiList) | !is.numeric(roiScales))
            stop("If provided, parameter 'roiScales' has to be a numeric with",
                 " length equal to the length of 'roiList'!")

    basenames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                   "into", "intb", "maxo", "sn")
    verbosenames <- c("egauss", "mu", "sigma", "h", "f", "dppm", "scale",
                      "scpos", "scmin", "scmax", "lmin", "lmax")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(scantime))) / 2)

    if (length(z <- which(scalerange == 0)))
        scalerange <- scalerange[-z]
    if (length(scalerange) < 1) {
        warning("No scales? Please check peak width!")
        if (verboseColumns) {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                            length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames))
            colnames(nopeaks) <- c(basenames)
        }
        return(invisible(nopeaks))
    }

    if (length(scalerange) > 1)
        scales <- seq(from = scalerange[1], to = scalerange[2], by = 2)
    else
        scales <- scalerange

    minPeakWidth <-  scales[1]
    noiserange <- c(minPeakWidth * 3, max(scales) * 3)
    maxGaussOverlap <- 0.5
    minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
    minCentroids <- minPtsAboveBaseLine
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth / 2)
    scanrange <- c(1, length(scantime))

    ## If no ROIs are supplied then search for them.
    if (length(roiList) == 0) {
        message("Detecting mass traces at ", ppm, " ppm ... ", appendLF = FALSE)
        ## flush.console();
        ## We're including the findmzROI code in this function to reduce
        ## the need to copy objects etc.
        ## We could also sort the data by m/z anyway; wouldn't need that
        ## much time. Once we're using classes from MSnbase we can be
        ## sure that values are correctly sorted.
        withRestarts(
            tryCatch({
                tmp <- capture.output(
                    roiList <- .Call("findmzROI",
                                     mz, int, scanindex,
                                     as.double(c(0.0, 0.0)),
                                     as.integer(scanrange),
                                     as.integer(length(scantime)),
                                     as.double(ppm * 1e-6),
                                     as.integer(minCentroids),
                                     as.integer(prefilter),
                                     as.integer(noise),
                                     PACKAGE ='xcms' )
                )
            },
            error = function(e){
                if (grepl("m/z sort assumption violated !", e$message)) {
                    invokeRestart("fixSort")
                } else {
                    simpleError(e)
                }
            }),
            fixSort = function() {
                ## Force ordering of values within spectrum by mz:
                ##  o split values into a list -> mz per spectrum, intensity per
                ##    spectrum.
                ##  o define the ordering.
                ##  o re-order the mz and intensity and unlist again.
                ## Note: the Rle split is faster than the "conventional" factor split.
                splitF <- Rle(1:length(valsPerSpect), valsPerSpect)
                mzl <- as.list(S4Vectors::split(mz, f = splitF))
                oidx <- lapply(mzl, order)
                mz <<- unlist(mapply(mzl, oidx, FUN = function(y, z) {
                    return(y[z])
                }, SIMPLIFY = FALSE, USE.NAMES = FALSE), use.names = FALSE)
                int <<- unlist(mapply(as.list(split(int, f = splitF)), oidx,
                                      FUN=function(y, z) {
                                          return(y[z])
                                      }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
                               use.names = FALSE)
                rm(mzl)
                rm(splitF)
                tmp <- capture.output(
                    roiList <<- .Call("findmzROI",
                                      mz, int, scanindex,
                                      as.double(c(0.0, 0.0)),
                                      as.integer(scanrange),
                                      as.integer(length(scantime)),
                                      as.double(ppm * 1e-6),
                                      as.integer(minCentroids),
                                      as.integer(prefilter),
                                      as.integer(noise),
                                      PACKAGE ='xcms' )
                )
            }
        )
        message("OK")
        ## ROI.list <- findmzROI(object,scanrange=scanrange,dev=ppm * 1e-6,minCentroids=minCentroids, prefilter=prefilter, noise=noise)
        if (length(roiList) == 0) {
            warning("No ROIs found! \n")
            if (verboseColumns) {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                                length(verbosenames))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames))
                colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
        }
    }
    return(roiList)
}


