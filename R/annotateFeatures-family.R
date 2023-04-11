#### annotateFeaturesFromList-related ####

#' @title annotateFeaturesFromList
#' @description Finds a correspondance between XCMS features and RHermes SOIs
#'   in order to filter the former and obtain a cleaner peak matrix
#' @param XCMSnExp XCMS object that contains a list of features (ie. after
#'   running `groupChromPeaks()`).
#' @param RHermesExp RHermes object that contains a SOI (Scans of Interest)
#'   list.
#' @param SOI_id Index of the SOI that will be used as reference for the matching.
#' @param MS2Exp Numeric. If provided, index of the MS2Exp object inside the
#'   RHermesExp object from which the MS2 scans will be imported and added to
#'   the feature table.
#' @param RTtol Numeric. Retention time shift tolerance in seconds.
#' @param filter Logical. Whether to remove those features that can't be matched
#'   to a SOI.
#' @param quantifyMissingSOI Logical. Whether to quantify via EIC (Extracted Ion
#'   Chromatogram) those SOIs that were not matched with a feature.
#' @importFrom RHermes SOI
#' @importFrom xcms featureDefinitions
#' @importFrom S4Vectors DataFrame
#' @export
annotateFeaturesFromList <- function(XCMSnExp, RHermesExp, SOI_id, MS2Exp = NA,
                                     RTtol = 10, filter = TRUE,
                                     quantifyMissingSOI = FALSE){
    
    #Match XCMS features with SOIs
    output <- matchPeaksToSOI(XCMSnExp, RHermesExp, SOI_id, RTtol, quantifySOI = quantifyMissingSOI)
    XCMSnExp <- output[[1]]
    pks <- output[[2]]
    message(length(which(!is.na(pks$soi))), " XCMS peaks with SOI")
    
    #Filter all features without matching SOI
    if(filter) pks <- subset(pks, !is.na(pks$soi))
    
    #Extract MS2 information and match with features
    if(!is.na(MS2Exp)){
        pks <- retrieveMS2Info(pks, RHermesExp, MS2Exp, RTtol)
        # message(length(which(pks$foundMS2)), " XCMS peaks with MS2")
        # message("   - Of which ",
        #         length(which(!pks$putativeID[!is.na(pks$soi) & pks$foundMS2] %in%
        #                          c("No significant hits", "Missing reference spectra"))),
        #         " have a putative ID")
    }
    
    #Update XCMSnExp object
    class(pks$peakidx) <- "list"
    # class(pks$soi) <- "numeric"
    if(!is.na(MS2Exp)) class(pks$putativeID) <- "character"
    featureDefinitions(XCMSnExp) <- DataFrame(pks)
    return(XCMSnExp)
}

#' @importFrom xcms featureDefinitions processHistory manualChromPeaks
#'   groupChromPeaks chromPeakData
#' @importFrom RHermes SOI
#' @importFrom dplyr between mutate rename select
#' @importFrom stats na.omit
#' 
matchPeaksToSOI <- function(XCMSnExp, RHermesExp, SOI_id, RTtol = 3, quantifySOI = FALSE) {
    pks <- as.data.frame(featureDefinitions(XCMSnExp))
    sois <- SOI(RHermesExp, SOI_id)@SOIList
    ppm <- RHermesExp@metadata@ExpParam@ppm
    
    ##Perform correspondance XCMS ChromPeak --> SOI 
    matching <- lapply(seq_len(nrow(pks)), function(i){
        mz_range <- pks$mzmed[i] * c(1 - ppm * 1e-6,
                                     1 + ppm * 1e-6)  
        id <- which(
            between(sois$mass, mz_range[1], mz_range[2]) &
                    sois$start < pks$rtmed[i] &
                    sois$end > pks$rtmed[i]
        )
        return(id)
    })
    pks$soi <- lapply(seq_along(matching), function(i){
        if (length(matching[[i]]) == 0) return(NA)
        matching[[i]]
    })
    
    #Matching summary
    representedSOI_ID <- na.omit(unique(unlist(pks$soi)))
    message(length(representedSOI_ID), " out of ", nrow(sois),
            " (", round(length(representedSOI_ID)/nrow(sois) * 100, 2),"%) ",
            "SOIs were represented in the feature list") 
    message(length(which(!is.na(pks$soi))), " XCMS features with SOI out of ",
                nrow(pks))
    
    #Fill unmatched SOI -- Optional
    if (quantifySOI) {
        message("Quantifying SOIs that could not be matched to a peak")
        target_soi <- sois[-representedSOI_ID, c("start", "end", "mass")] %>%
            mutate(mzmin = mass * (1 - 1e-6 * ppm),
                   mzmax = mass * (1 + 1e-6 * ppm)) %>%
            rename(rtmin = start, rtmax = end) %>%
            select("mzmin", "mzmax", "rtmin", "rtmax")
        
        pdp_idx <- max(which(sapply(processHistory(XCMSnExp),
                                    function(x){x@type == "Peak grouping"})))
        pdp <- processHistory(XCMSnExp)[[pdp_idx]]@param
        
        npeak <- nrow(chromPeaks(XCMSnExp))
        XCMSnExp <- manualChromPeaks(XCMSnExp, target_soi)
        npeak_new <- nrow(chromPeaks(XCMSnExp))
        chromPeakData(XCMSnExp)[seq(npeak + 1, npeak_new), "is_filled"] <- TRUE
        XCMSnExp <- groupChromPeaks(XCMSnExp, pdp)
        
        pks <- as.data.frame(featureDefinitions(XCMSnExp))
        matching <- lapply(1:nrow(pks), function(i){
            mz_range <- pks$mzmed[i] * c(1 - ppm * 1e-6,
                                         1 + ppm * 1e-6)  
            id <- which(
                between(sois$mass, mz_range[1], mz_range[2]) &
                    sois$start < pks$rtmed[i] &
                    sois$end > pks$rtmed[i]
            )
            return(id)
        })
        pks$soi <- lapply(seq_along(matching), function(i){
            if(length(matching[[i]])==0){return(NA)}
            matching[[i]]
        })
        
        pks$filled <- sapply(pks$peakidx, function(id, cpd){
            all(cpd[id, "is_filled"])
        }, cpd = chromPeakData(XCMSnExp))
        
        representedSOI_ID <- na.omit(unique(unlist(pks$soi)))
        message(length(representedSOI_ID), " out of ", nrow(sois),
                " (", round(length(representedSOI_ID)/nrow(sois)*100, 2),"%) ",
                "SOIs were represented in the feature list after filling")
    } else {
        pks$filled <- FALSE
    }
    pks$formula <- lapply(seq_along(matching), function(i){
        if (length(matching[[i]]) == 0) {return(NA)}
        sois$formula[matching[[i]]]
    })
    pks$isotope <- lapply(seq_along(matching), function(i){
        if (length(matching[[i]]) == 0) {return(NA)}
        return("M0")
    })
    
    #Add confirmed isotope annotations to feature list
    if ("isodf" %in% colnames(sois)) {
        for (i in which(sois$isofound != 0)){
            iso_table <- sois$isodf[[i]]
            for (j in seq_len(nrow(iso_table))){
                mz_range <- iso_table$mass[j] * c(1 - ppm * 1e-6,
                                                  1 + ppm * 1e-6) 
                hits <- which(between(pks$mzmed, mz_range[1], mz_range[2]) &
                              between(pks$rtmed, sois$start[[i]], sois$end[[i]]))
                pks$isotope[hits] <- iso_table$iso[[j]]
            }
        }
    }
    pks$isotope <- as.character(pks$isotope)
    pks$SOIData <- lapply(matching, function(i){
        if(length(i) == 0){return(NA)}
        sois[i,]
    })
    return(list(XCMSnExp, pks))
}


#' @importFrom dplyr between 
retrieveMS2Info <- function(pks, RHermesExp, MS2ExpID, RTtol = 10){
    MS2Features <- RHermesExp@data@MS2Exp[[MS2ExpID]]@Ident$MS2Features
    # ppm <- RHermesExp@metadata@ExpParam@ppm
    matching <- lapply(seq_len(nrow(pks)), function(i){
        # mz_range <- pks$mz[i] * c(1 - ppm * 1e-6,
        #                           1 + ppm * 1e-6)  
        mz_range <- c(pks$mzmed[i] - 0.2,
                      pks$mzmed[i] + 0.2)  
        id <- which(
            between(MS2Features$precmass, mz_range[1], mz_range[2]) &
                abs(MS2Features$apex - pks$rtmed[i]) < RTtol
        )
        return(id)
    })
    pks$foundMS2 <- FALSE
    for (i in seq_along(matching)){
        if (length(matching[[i]]) == 0) next
        pks$foundMS2[i] <- TRUE
    }
    pks$MS2Data <- lapply(matching, function(i){
        if (length(i) == 0) {return(NA)}
        MS2Features[i,]
    })
    pks$putativeID <- lapply(pks$MS2Data, function(x){
        if (!is.na(x)[[1]]) {
            if (is.data.frame(x$results[[1]])) {
                return(strsplit(x$results[[1]]$formula[[1]], "#")[[1]][3])
            } else {
                return(x$results[[1]])
            }
        } else {
            return("No available ID")
    }})
    message(length(which(!is.na(pks$soi) & pks$foundMS2)),
            " XCMS features with SOI and MS2")
    message("   - Of which ",
            length(which(!pks$putativeID[!is.na(pks$soi) & pks$foundMS2] %in%
                             c("No significant hits", "Missing reference spectra"))),
            " have a putative ID")
    return(pks)
}

#### annotateFeatures-related ####

#' @title annotateFeatures
#' @description Annotate features using a list of molecular formulas and adducts.
#' @inheritParams annotateFeaturesFromList
#' @param ChemFormulaParam A ChemFormulaParam object that contains the lists of
#'   molecular formulas and adducts to annotate the features.
#' @param filter Logical. Whether to remove those features that can't be matched
#'   to a SOI.
#' @importFrom xcms featureDefinitions
#' @importFrom S4Vectors DataFrame
#' @export
annotateFeatures <- function(XCMSnExp, ChemFormulaParam, filter = TRUE){
    #Match XCMS features with annotation
    feature <- featureDefinitions(XCMSnExp)
    pks <- .matchPeaksCFP(feature, ChemFormulaParam)

    #Filter all features without matching annotation
    if(filter) pks <- subset(pks, !vapply(pks$anot, is.null,
                             FUN.VALUE = logical(1)))
    
    #Update XCMSnExp object
    featureDefinitions(XCMSnExp) <- DataFrame(pks)
    return(XCMSnExp)
}

.matchPeaksCFP <- function(feature, ChemFormulaParam){
    ppm <- ChemFormulaParam@ppm
    ionf <- ChemFormulaParam@ionFormulas
    ionf <- ionf[order(ionf$m),]
    match <- vector("list", nrow(feature))
    for (i in seq_len(nrow(feature))) {
        mz <- feature$mzmed[[i]]
        anot <- ionf$f[RHermes:::binarySearch(matrix(ncol = 1, ionf$m), mz, ppm)]
        if (length(anot)) match[[i]] <- anot
    }
    feature$anot <- match
    return(feature)
}
