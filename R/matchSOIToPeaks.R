#' @export
mergeRHermesXCMS <- function(XCMSnExp, RHermesExp, SOI, MS2Exp = NA, RTtol = 10){
    pks <- as.data.frame(chromPeaks(XCMSnExp))
    
    #Match XCMS peaks with SOIs
    pks <- matchPeaksToSOI(pks, RHermesExp, SOI, RTtol)
    message(length(which(!is.na(pks$soi))), " XCMS peaks with SOI")
    
    #Extract MS2 information and match with peaks
    if(!is.na(MS2Exp)){
        pks <- subset(pks, !is.na(pks$soi))
        pks <- retrieveMS2Info(pks, RHermesExp, MS2Exp, RTtol)
        # message(length(which(pks$foundMS2)), " XCMS peaks with MS2")
        message(length(which(!is.na(pks$soi) & pks$foundMS2)),
                " XCMS peaks with SOI and MS2")
        message("   - Of which ",
                length(which(!pks$putativeID[!is.na(pks$soi) & pks$foundMS2] %in%
                                 c("No significant hits", "Missing reference spectra"))),
                " have a putative ID")
    }
    
    #Filter all peaks without matching SOI
    pks <- subset(pks, !is.na(pks$soi))
    
    #Update XCMSnExp object
    names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
               "intb", "maxo", "sn", "sample")
    cp_mat <- pks[, names] %>% as.matrix() %>% apply(., 2, as.numeric)
    row.names(cp_mat) <- row.names(pks)
    chromPeaks(XCMSnExp) <- cp_mat
    otherColumns <- colnames(pks)[!colnames(pks) %in% names]
    chromPeakData(XCMSnExp) <- cbind(chromPeakData(XCMSnExp),
                                   S4Vectors::DataFrame(pks[, otherColumns],
                                                  row.names = row.names(pks)))
    return(XCMSnExp)
}

#' @title matchPeaksToSOI
#' @description Finds a correspondance between XCMS Chrompeaks and RHermes SOIs
#'   in order to filter the former and obtain a cleaner peak matrix
#' @param pks Peak matrix (ie. chromPeaks()) from XCMS object
#' @param RHermesExp RHermes object
#' @param SOI_id Index of the SOI that will be used as reference for the
#'   matching
#' 
matchPeaksToSOI <- function(pks, RHermesExp, SOI_id, RTtol = 3){
    sois <- SOI(RHermesExp, SOI_id)@SOIList
    ppm <- RHermesExp@metadata@ExpParam@ppm
    
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
    # pks$soi <- NA
    # pks$formula <- NA
    # for(i in which(sapply(matching, length) != 0)){
    #     # if(length(matching[[i]])==0) next
    #     pks$soi[i] <- list(matching[[i]])
    #     pks$formula[i] <- list(sois$formula[matching[[i]]])
    # }
    pks$soi <- lapply(seq_along(matching), function(i){
        if(length(matching[[i]])==0){return(NA)}
        matching[[i]]
    })
    pks$formula <- lapply(seq_along(matching), function(i){
        if(length(matching[[i]])==0){return(NA)}
        sois$formula[matching[[i]]]
    })
    pks$SOIData <- lapply(matching, function(i){
        if(length(i) == 0){return(NA)}
        sois[i,]
    })
    
    representedSOI_ID <- na.omit(unique(unlist(pks$soi)))
    message(length(representedSOI_ID), " out of ", nrow(sois),
            " (", round(length(representedSOI_ID)/nrow(sois)*100, 2),"%) ",
            "SOIs were represented in the peak list")
    return(pks)
}

retrieveMS2Info <- function(pks, RHermesExp, MS2ExpID, RTtol = 10){
    MS2Features <- RHermesExp@data@MS2Exp[[MS2ExpID]]@Ident$MS2Features
    ppm <- RHermesExp@metadata@ExpParam@ppm
    matching <- lapply(1:nrow(pks), function(i){
        # mz_range <- pks$mz[i] * c(1 - ppm * 1e-6,
        #                           1 + ppm * 1e-6)  
        mz_range <- c(pks$mz[i] - 0.2,
                      pks$mz[i] + 0.2)  
        id <- which(
            between(MS2Features$precmass, mz_range[1], mz_range[2]) &
                abs(MS2Features$apex - pks$rt[i]) < RTtol
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
    pks$putativeID <- lapply(pks$MS2Data, function(x){
        if(!is.na(x)[[1]]){
            if(is.data.frame(x$results[[1]])){
                return(strsplit(x$results[[1]]$formula[[1]], "#")[[1]][3])
            } else {
                return(x$results[[1]])
            }
        } else {
            return("No available ID")
    }})
    return(pks)
}

getQuantitativeMatrix <- function(XCMSnExp){
    ft <- as.data.frame(featureDefinitions(XCMSnExp))
    ft <- cbind(ft, featureValues(XCMSnExp))
    pks <- as.data.frame(chromPeakData(XCMSnExp))
    pk_data <- lapply(seq(nrow(ft)), function(row){
        peakid <- ft$peakidx[[row]]
        cur_pks <- pks[peakid,]
        return(list(
            soi = cur_pks$soi[[1]],
            formula = cur_pks$formula[[1]],
            putativeID = cur_pks$putativeID[[1]]
        ))
    })
    pk_data <- do.call(rbind, pk_data)
    ft <- cbind(pk_data, ft)
    # ft$soi <- unlist(ft$soi)
    # ft$formula <- unlist(ft$formula)
    # ft$putativeID <- unlist(ft$putativeID)
    return(ft)
}

extract_MS2_from_feature <- function(XCMSnExp_THS, feature_df){
    pkdata <- as.data.frame(chromPeakData(XCMSnExp_THS))
    feature_df$MS2data <- lapply(feature_df$peakidx, function(peaks){
        ms2_df <- do.call(rbind, pkdata$MS2Data[peaks])
        if(!is.data.frame(ms2_df)){return(NA)}
        ms2_df <- ms2_df[!duplicated(ms2_df$apex),]
        return(ms2_df)
    })
    return(feature_df)
}

calculate_MS2_correlation_from_feature <- function(feature_df){
    all_ms2 <- lapply(seq(nrow(feature_df)), function(x){
        ms2_df <- feature_df$MS2data[[x]]
        if(is.na(ms2_df)){return()}
        ms2_df$originalFeature <- x
        return(ms2_df)
    })
    all_ms2 <- do.call(rbind, all_ms2)
    all_ms2 <- all_ms2[sapply(all_ms2$ssdata, is.data.frame), ]
    
    dist <- pblapply(seq(nrow(all_ms2)), function(x){
        query <- all_ms2$ssdata[[x]]
        sapply(seq(nrow(all_ms2)), function(y){
            if(y >= x) return(0)
            pattern <- all_ms2$ssdata[[y]]
            RHermes:::calculate_MS2_distance(query, pattern, method = "cosine")
        })
    })
    dist_matrix <- do.call(rbind, dist)
    net <- graph_from_adjacency_matrix(dist_matrix, mode = "undirected")
    net <- simplify(net)
    netd3 <- igraph_to_networkD3(net, group = membership(cluster_fast_greedy(net)))
    
    newscript <- 'alert("mz: " + (d.var1.mz) + "\\n" + "Int:" + (d.var1.int));'
    
    fn <- forceNetwork(Links = netd3$links, Nodes = netd3$nodes, 
                 Source = 'source', Target = 'target', 
                 NodeID = 'name', Group = 'group',
                 zoom = TRUE, opacity = 1, clickAction = newscript)
    fn$x$nodes$var1 <- lapply(all_ms2$ssdata, function(ss){round(ss)})
    fn 
    
    fts_cluster <- all_ms2 %>% 
        group_by(gr = membership(cluster_fast_greedy(net))) %>%
        summarize(size = n(), ft = list(unique(originalFeature)))
    feature_groups <- lapply(fts_cluster$ft, function(x){feature_df[x,]})
    feature_groups <- feature_groups[sapply(feature_groups, nrow) > 1]
    lapply(feature_groups, function(group){
        plots <- lapply(1:nrow(group), function(id){plot_MS2_from_feature(group, id)})
    })
}

#' @export
mergeRHermesXCMSFeatures <- function(XCMSnExp, RHermesExp, SOI, MS2Exp = NA, RTtol = 10){
    pks <- as.data.frame(featureDefinitions(XCMSnExp))
    
    pks <- rename(pks, mz = mzmed, rt = rtmed)
    
    #Match XCMS features with SOIs
    pks <- matchPeaksToSOI(pks, RHermesExp, SOI, RTtol)
    message(length(which(!is.na(pks$soi))), " XCMS features with SOI")
    
    #Extract MS2 information and match with features
    if(!is.na(MS2Exp)){
        pks <- subset(pks, !is.na(pks$soi))
        pks <- retrieveMS2Info(pks, RHermesExp, MS2Exp, RTtol)
        # message(length(which(pks$foundMS2)), " XCMS peaks with MS2")
        message(length(which(!is.na(pks$soi) & pks$foundMS2)),
                " XCMS features with SOI and MS2")
        message("   - Of which ",
                length(which(!pks$putativeID[!is.na(pks$soi) & pks$foundMS2] %in%
                                 c("No significant hits", "Missing reference spectra"))),
                " have a putative ID")
    }
    
    #Filter all features without matching SOI
    pks <- subset(pks, !is.na(pks$soi))
    
    
    pks <- rename(pks, mzmed = mz, rtmed = rt)
    
    #Update XCMSnExp object
    class(pks$peakidx) <- "list"
    # class(pks$soi) <- "numeric"
    if(!is.na(MS2Exp)) class(pks$putativeID) <- "character"
    featureDefinitions(XCMSnExp) <- S4Vectors::DataFrame(pks)
    return(XCMSnExp)
}

