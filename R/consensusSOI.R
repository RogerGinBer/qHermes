#'@title consensusSOI
#'@description Group a set of SOI lists into a consensus list that contains all
#'  SOIs common in a subset of them. The concept is similar to the grouping step
#'  from peaks to features present in other tools.
#'@param RHermesExp An RHermesExp object containing the SOI lists
#'@param SOIids Which SOI lists are to be grouped.
#'@param minSOI Minimum number of SOI lists where a particular SOI has to be
#'  present in order to be preserved in the final list.
#'@param rtwin Retention time tolerance to match the SOIs from the different
#'  lists.
#'@details consensusSOI is commonly used when multiple quality control samples
#'  are available in order to reduce the number of SOIs and ensure only
#'  consistent signals are included into the inclusion list.
#'@examples 
#'if(F){
#'  #Find SOIs common in at least two of the three selected SOI lists
#'  consensusSOI(RHermesExp, 1:3, 2) 
#' }
#'@export
#'@importFrom data.table rbindlist
#'@importFrom dplyr mutate
consensusSOI <- function(RHermesExp, SOIids, minSOI = NULL, rtwin = 5){
    if(length(SOIids) == 1){stop("Can't run consensus with a single SOI list")}
    if(is.null(minSOI)){
        minSOI <- ceiling(length(SOIids) * 0.8)
        message("No value selected for minSOI. Using 80% rule (minSOI = ",
                minSOI, ")")
    }
    if(length(SOIids) < minSOI){
        stop("minSOI value too large, please select a valid value for minSOI")
    }
    slist <- lapply(SOIids, SOI, struct = RHermesExp)
    sdf <- lapply(seq(length(slist)), function(x, id){
        r <- slist[[x]]@SOIList
        r$soiN <- id[x]
        return(r)
    }, id = SOIids)
    sdf <- data.table::rbindlist(sdf, fill = TRUE)
    sdf <- dplyr::mutate(sdf,
                   apex = as.numeric(
                       lapply(peaks, function(x){x$rt[which.max(x$rtiv)]}))
    )
    
    message("Calculating consensus SOIs")
    RES <- lapply(seq(nrow(sdf)),function(x) return())
    exlist <- c()
    for(i in seq(nrow(sdf))){
        if(i %in% exlist){next}
        x <- sdf[i,]
        idx <- which(abs(sdf$apex - x$apex) < rtwin & sdf$formula == x$formula)
        if(any(idx %in% exlist)) {idx <- idx[-which(idx %in% exlist)]}
        usoi <- sdf[idx,]
        if(length(unique(usoi$soiN)) >= minSOI){
            newsoi <- usoi[which.max(usoi$MaxInt),]
            newsoi$start <- min(usoi$start) #mean?
            newsoi$end <- max(usoi$end) #mean?
            newsoi$length <- newsoi$end - newsoi$start
            # sum $peaks of all to unify?
            # newsoi <- newsoi[,1:18]
            newsoi$soiN <- paste(usoi$soiN, collapse="")
            newsoi$nsample <- length(unique(usoi$soiN))
            exlist <- c(exlist, idx)
            # too slow? potser eliminar de sdf directament?
            RES[[i]] <- newsoi
        }
    }
    RES <- data.table::rbindlist(RES, fill = TRUE)
    
    message("Recalculating SOI plotting dataframe")
    RES <- data.table::data.table(RES) 
    data.table::setkey(RES,formula)
    plist <- lapply(unique(RES$formula), RHermes:::preparePlottingDF, RES)
    plist <- do.call(rbind, plist)
    plist$isov <- rep("M0", nrow(plist))    
    
    ## Create a new SOI list object
    nsoi <- length(RHermesExp@data@SOI)
    originalNames <- vapply(slist, function(x) x@filename, character(1))
    RHermesExp@data@SOI[[nsoi + 1]] <- RHermesSOI(
        SOIList = RES,
        PlotDF = as.data.table(plist),
        SOIParam = SOIParam(),
        filename = unique(originalNames)
    )
    
    return(RHermesExp)
}
