#'@export
findSOIGroups <- function(MSnExp, sampleGroups = NULL, minfrac = 0.8){
    if(is.null(sampleGroups)){
        stop("Please introduce your sample class as a vector in sampleGroups")
    }
    if(!between(minfrac, 0, 1)) stop("minfrac should be between 0 and 1")
    
    # maybe we will need to add sampleclass info as in standardxcms for MsnExp
    pks <- data.table(xcms::chromPeaks(MSnExp))
    # target_list <- RHermes::SOI(struct,1)@SOIList
    # can add the SOIList from a RHermesExp or not ? 

    filenames <- tools::file_path_sans_ext(MSnbase::sampleNames(MSnExp))
    # consider parallelize this
    df <- bplapply(unique(pks$SOIidx),
                   function(n, pks, filenames, sampleGroups){
        peakidx <- which(pks$SOIidx==n) 
        if(length(peakidx)==0){return()}
        npeaks <- pks[peakidx,]
        npeaks$peakidx <- peakidx

        grped <- lapply(seq(length(filenames)), function(i){
            nipks <- npeaks[which(npeaks$sample==i), ]
            if(nrow(nipks) > 1){
                nres <- nipks[1,]
                nres$mz <- mean(nipks$mz)
                nres$mzmin <- min(nipks$mzmin)
                nres$mzmax <- max(nipks$mzmax)
                nres$rt <- mean(nipks$rt)
                nres$rtmin <- min(nipks$rtmin)
                nres$rtmax <- max(nipks$rtmax)
                nres$into <- sum(nres$into)
                nres$intb <- sum(nres$intb)
                nres$maxo <- max(nres$into)
                nres$npeaks <- nrow(nipks)
            }
            if(nrow(nipks) == 1){
                nres <- nipks[1, ]
                nres$npeaks <- 1
            }
            if(nrow(nipks) == 0){return()}
            nres$peakidx <- list(nipks$peakidx)
            nres <- dplyr::rename(nres, mzmed = mz, rtmed = rt)
        })
        grped <- do.call("rbind",grped)
        grped$sampclass <- sampleGroups[grped$sample]
        sampleGroupNames <- unique(sampleGroups)
        sampleGroupTable <- table(sampleGroups)
        col_nms <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                     "npeaks", sampleGroupNames)
        res_mat <- matrix(nrow = 0, ncol = length(col_nms),
                          dimnames = list(character(), col_nms))
        gcount <- rep(0, length(sampleGroupNames))
        names(gcount) <- sampleGroupNames
        tt <- table(grped$sampclass)
        gcount[names(tt)] <- as.numeric(tt)
        res_mat <- rbind(res_mat,
                         c(median(grped$mzmed),
                           min(grped$mzmin),
                           max(grped$mzmax),
                           median(grped$rtmed),
                           min(grped$rtmin),
                           max(grped$rtmax),
                           sum(grped$npeaks),
                           gcount)
        )
        res <- as.data.frame(res_mat)
        res$peakidx <- list(sort(unlist(grped$peakidx)))
        return(res)
    }, BPPARAM = bpparam(),
    pks = pks, filenames = filenames, sampleGroups = sampleGroups)
    df <- do.call("rbind", df)
    
    #Minimal fraction filtering
    filtgroup <- "Sample"
    sg <- t(as.matrix(table(sampleGroups)))
    uniqueClass <- colnames(sg)
    minfrac <- apply(df[, uniqueClass, drop = FALSE], 1, function(x)
        any(sapply(uniqueClass,
                   function(y) as.numeric(x[y])/sg[,y]) >= minfrac)
    )
    df <- df[minfrac,]

    featureDefinitions(MSnExp) <- S4Vectors::DataFrame(
        df,
        row.names = seq(nrow(df)) # rename this as SOIxx?
    )
    return(MSnExp)
}