matchSOIToPeaks <- function(XCMSnExp, RHermesExp, SOI_id){
    
    pks <- as.data.frame(chromPeaks(XCMSnExp))
    sois <- SOI(RHermesExp, SOI_id)@SOIList
    ppm <- RHermesExp@metadata@ExpParam@ppm
    tol <- 10
    
    ##Perform correspondance SOI --> XCMS ChromPeak
    matching <- lapply(1:nrow(sois), function(i){
        mz_range <- sois$mass[i] * c(1-ppm*1e-6, 1+ppm*1e-6)  
        rt_range <- c(sois$start[i]-tol, sois$end[i]+tol)
        id <- which(
            between(pks$mz, mz_range[1], mz_range[2]) &
            pks$rtmin > rt_range[1] &
            pks$rtmax < rt_range[2]
        )
        return(id)
    })
    
    ##Add information to ChromPeak matrix
    pks$soi <- NA
    pks$formula <- NA
    for(i in seq_along(matching)){
        if(length(matching[[i]])==0) next
        pks$soi[matching[[i]]] <- i
        pks$formula[matching[[i]]] <- sois$formula[i]
    }
    
    newpks <- subset(pks, !is.na(pks$soi))
    
    representedSOI_ID <- na.omit(unique(unlist(pks$soi)))
    message(length(representedSOI_ID), " out of ", nrow(sois),
            " (",round(length(representedSOI_ID)/nrow(sois)*100,2),"%) ",
            "SOIs were represented in the peak list")
    
    representedSOI <- sois[representedSOI_ID,]
    unrepresentedSOI <- sois[-representedSOI_ID,]
    return(list(pks = newpks,
                repSOI = representedSOI,
                norepSOI = unrepresentedSOI))
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

