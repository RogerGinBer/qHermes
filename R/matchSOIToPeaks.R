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
    chromPeakData(XCMSnExp) <- cbind(chromPeakData(XCMSnExp),
                                   S4Vectors::DataFrame(pks[, 12:17],
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



plotpks <- function(peakMat){
    ggplot(peakMat) + 
        geom_segment(aes(x=rtmin,xend=rtmax,y=mz,yend=mz, color = log10(maxo)))+
        geom_point(aes(x=rt,y=mz), size = 0.1, color = "tomato", alpha = 0.2)+
        facet_grid(sample~.) +
        theme_minimal()
}

plot_particular_peak <- function(id, XCMSnExp){
    
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

plot_MS2_from_feature <- function(feature_df, id){
    entryMS2 <- feature_df$MS2data[[id]]
    molecmass <- entryMS2$precmass[[1]]
    lapply(entryMS2$ssdata, function(query){
        plot_MS2_from_ssdata(query, molecmass)
    })
}

plot_MS2_from_ssdata <- function(query, molecmass){
    maxint <- max(query$int)
    query$int <- query$int/maxint * 100
    bestdf <- query[query$int > 10, ]
    bestdf$mz <- round(bestdf$mz, 4)
    moldf <- data.frame(mz = molecmass)
    
    subtitle <- ""
    title <- ""
    
    pl <- ggplot() +
        geom_segment(data = query, aes(x = .data$mz, xend = .data$mz,
                                       y = 0, yend = .data$int),
                     color = "black") +
        geom_point(data = moldf, aes(x = .data$mz, y = 0), shape = 17,
                   size = 2) +
        theme_minimal() + ylab("Relative Intensity") +
        theme(plot.margin = unit(c(1, 0.7, 1, 0.8), "cm"),
              text = element_text(size = 11, family = "Segoe UI Light"),
              plot.title = element_text(hjust = 0.5)) +
        geom_text(data = bestdf, aes(x = .data$mz, y = .data$int + 5,
                                     label = .data$mz),
                  family = "Segoe UI Light", check_overlap = TRUE) +
        scale_x_continuous(limits = c(min(query$mz, molecmass) - 20,
                                      max(query$mz, molecmass) + 20))+
        ggtitle(title, subtitle)
    
    # ggplotly(pl, height = 400)
    return(pl)
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


# library(dplyr)
# newpks <- newpks[order(newpks$maxo, decreasing = T),]
# i <- 5
# whichSOI <- newpks$soi[[i]]
# ms2plot <- ggplot(newpks$MS2Data[[i]]$ssdata[[1]]) +
#     geom_segment(aes(x=mz,xend=mz,y=0,yend=int)) + theme_minimal()
# 
# ggplot(sois$peaks[[whichSOI]]) +
#     geom_point(aes(x=rt, y = log10(rtiv))) +
#     geom_point(aes(x=rt, y = log10(rtiv)), color = "red",
#                data = filter(sois$peaks[[whichSOI]],
#                              between(rt, newpks$rtmin[[i]], newpks$rtmax[[i]])))+
#     geom_vline(aes(xintercept=rt, color = rtiv),
#                data = newpks$MS2Data[[i]]$profile[[1]], lty = 3) +
#     annotation_custom(grob=ggplotGrob(ms2plot), 
#                       ymin = log10(newpks$maxo[[i]])*2/3,
#                       ymax=log10(newpks$maxo[[i]]),
#                       xmin=(sois$end[[whichSOI]]-sois$start[[whichSOI]])/2 +
#                           sois$start[[whichSOI]],
#                       xmax=sois$end[[whichSOI]]) +
#     ggtitle(paste(newpks$formula[[i]], round(newpks$mz,4)),
#             subtitle = paste("Putative identification:",
#                              strsplit(newpks$MS2Data[[i]]$results[[1]]$formula[[1]], "#")[[1]][3]))

