findgroupsfromSOI <- function(MsnExp, SOIList, sampleGroups=NULL, rtwin=10){
  if(is.null(sampleGroups)){stop("Please introduce your sample class as a vector in sampleGroups")}
  # maybe we will need to add sampleclass info as in standardxcms for MsnExp
  pks <- data.table(xcms::chromPeaks(MsnExp))
  # target_list <- RHermes::SOI(struct,1)@SOIList
  target_list <- SOIList
  sampleGroups <- as.character(sampleGroups)
  sampleGroupNames <- unique(sampleGroups)
  sampleGroupTable <- table(sampleGroups)
  nSampleGroups <- length(sampleGroupTable)
  filenames <- tools::file_path_sans_ext(xcms::sampleNames(MsnExp))
  # consider parallelize this
  x <- lapply(seq(length(unique(pks$SOIidx))),function(n){
    n <- unique(pks$SOIidx)[n]
    peakidx <- which(pks$SOIidx==n)
    npeaks <- pks[peakidx,]
    npeaks$peakidx <- peakidx
    nsoi <- target_list[n,]
    
    grped <- lapply(seq(length(filenames)),function(i){
      nipks <- npeaks[which(npeaks$sample==i),]
      # if(nrow(nipks)>1){ 
      #   nipks <- nipks[which(nipks$rtmin>=(nsoi$start-rtwin) &
      #                          nipks$rtmax<=(nsoi$end+rtwin)),]
      #   
      # }
      if(nrow(nipks)>1){
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
      if(nrow(nipks)==1){
        nres <- nipks[1,]
        nres$npeaks <- 1
      }
      if(nrow(nipks)==0){return()}
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
  })
  df <- do.call("rbind",x)
  featureDefinitions(MSnExp) <- S4Vectors::DataFrame(
    df,
    row.names = seq(nrow(df)) # rename this as SOIxx
  )
  return(MSnExp)
}