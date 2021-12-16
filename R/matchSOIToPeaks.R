#' @title matchPeaksToSOI
#' @description Finds a correspondance between XCMS Chrompeaks and RHermes SOIs
#'   in order to filter the former and obtain a cleaner peak matrix
#' @param XCMSnExp XCMS object
#' @param RHermesExp RHermes object
#' @param SOI_id Index of the SOI that will be used as reference for the
#'   matching
#' @export
#' 
matchPeaksToSOI <- function(XCMSnExp, RHermesExp, SOI_id){
    
    pks <- as.data.frame(chromPeaks(XCMSnExp))
    sois <- SOI(RHermesExp, SOI_id)@SOIList
    ppm <- RHermesExp@metadata@ExpParam@ppm
    tol <- 3 #Move to params
    
    ##Perform correspondance SOI --> XCMS ChromPeak
    # matching <- lapply(1:nrow(sois), function(i){
    #     mz_range <- sois$mass[i] * c(1-ppm*1e-6, 1+ppm*1e-6)  
    #     rt_range <- c(sois$start[i]-tol, sois$end[i]+tol)
    #     id <- which(
    #         between(pks$mz, mz_range[1], mz_range[2]) &
    #         pks$rtmin > rt_range[1] &
    #         pks$rtmax < rt_range[2]
    #     )
    #     return(id)
    # })
    
    ##Perform correspondance XCMS ChromPeak --> SOI 
    matching <- lapply(1:nrow(pks), function(i){
            mz_range <- pks$mz[i] * c(1 - ppm * 1e-6,
                                      1 + ppm * 1e-6)  
            id <- which(
                between(sois$mass, mz_range[1], mz_range[2]) &
                sois$start < pks$rt[i] &
                sois$end > pks$rt[i]
            )
            return(id)
        })
    
    ##Add information to ChromPeak matrix
    pks$soi <- NA
    pks$formula <- NA
    for(i in seq_along(matching)){
        if(length(matching[[i]])==0) next
        pks$soi[i] <- list(matching[[i]])
        pks$formula[i] <- list(sois$formula[matching[[i]]])
    }
    pks$SOIData <- lapply(matching, function(i){
        if(length(i) == 0){return(NA)}
        sois[i,]
    })
    
    newpks <- subset(pks, !is.na(pks$soi))
    
    representedSOI_ID <- na.omit(unique(unlist(pks$soi)))
    message(length(representedSOI_ID), " out of ", nrow(sois),
            " (",round(length(representedSOI_ID)/nrow(sois)*100,2),"%) ",
            "SOIs were represented in the peak list")
    
    representedSOI <- sois[representedSOI_ID,]
    unrepresentedSOI <- sois[-representedSOI_ID,]
    return(list(pks = pks,
                newpks = newpks,
                repSOI = representedSOI,
                norepSOI = unrepresentedSOI))
}

retrieveMS2Info <- function(XCMSnExp, RHermesExp, MS2ExpID){
    MS2Features <- RHermesExp@data@MS2Exp[[MS2ExpID]]@Ident$MS2Features
    pks <- as.data.frame(chromPeaks(XCMSnExp))
    ppm <- RHermesExp@metadata@ExpParam@ppm
    matching <- lapply(1:nrow(pks), function(i){
        # mz_range <- pks$mz[i] * c(1 - ppm * 1e-6,
        #                           1 + ppm * 1e-6)  
        mz_range <- c(pks$mz[i] - 0.2,
                      pks$mz[i] + 0.2)  
        id <- which(
            between(MS2Features$precmass, mz_range[1], mz_range[2]) &
                abs(MS2Features$apex - pks$rt[i]) < 5
        )
        return(id)
    })
    pks$foundMS2 <- FALSE
    for(i in seq_along(matching)){
        if(length(matching[[i]])==0) next
        pks$foundMS2[i] <- TRUE
    }
    pks$MS2Data <- lapply(matching, function(i){
        if(length(i) == 0){return(NA)}
        MS2Features[i,]
    })
    
}


plotpks <- function(peakMat){
    ggplot(peakMat) + 
        geom_segment(aes(x=rtmin,xend=rtmax,y=mz,yend=mz, color = log10(maxo)))+
        geom_point(aes(x=rt,y=mz), size = 0.1, color = "tomato", alpha = 0.2)+
        facet_grid(sample~.) +
        theme_minimal()
}

plot_particular_peak <- function(id, XCMSnExp){
    
}



