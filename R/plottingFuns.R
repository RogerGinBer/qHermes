#'@export
plotFeature <- function(XCMSnExp, feature, soi = NULL){
    cp <- chromPeaks(XCMSnExp)
    cpdata <- chromPeakData(XCMSnExp)%>% as.data.frame
    if(is.null(soi)){
        ft <- featureDefinitions(XCMSnExp)[feature,] %>% as.data.frame
        pks <- cpdata[ft$peakidx[[1]], ] 
        pknames <- fileNames(XCMSnExp)[cp[ft$peakidx[[1]], "sample"]] %>%
            basename
        mint <- ft$minrt[1]
    } else {
        pks <- dplyr::filter(cpdata, SOIidx == soi)
        if(nrow(pks) == 0){return()}
        pknames <- fileNames(XCMSnExp)[pks$sample] %>% basename
        mint <- min(pks$rtmin)
    }
    
    plot_data <- lapply(1:nrow(pks), function(x){
        pk <- pks$peaks[[x]]
        pk <- pk[pk$rtiv != 0,]
        pk$samp <- pknames[x]
        pk <- rbind(data.frame(rt = min(pk$rt), rtiv = 0, samp = pknames[x]),
                    pk,
                    data.frame(rt = max(pk$rt), rtiv = 0, samp = pknames[x]))
    })
    plot_data <- do.call(rbind, plot_data)
    
    p1 <- ggplotly(ggplot(plot_data) +
                       geom_polygon(aes(x=rt+mint, y=rtiv, fill = samp),
                                    alpha = 0.5) +
                       xlab("Retention time (s)") + ylab("Intensity") +
                       ggtitle(pks$formula[1]))
    
    integral <- plot_data %>% group_by(samp) %>%
        summarise(intensity = sum(rtiv))
    
    p2 <- ggplotly(ggplot(integral) +
                       geom_point(aes(x=intensity, y=0, color = samp),
                                  show.legend = FALSE)+
                       scale_x_log10() +
                       theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank())) %>% hide_legend()
    for (i in seq_len(nrow(integral))){
        p2$x$data[[i]]$text <- c(p2$x$data[[i]]$text, "") 
        p2$x$data[[i]]$showlegend <- FALSE
    }
    subplot(p1, p2, nrows = 2, heights = c(0.9, 0.1), which_layout = 1)
} 

#' @export
plotFeaturesMerged <- function(XCMSnExp, id){
    ft <- as.data.frame(featureDefinitions(XCMSnExp))
    pks <- chromPeaks(XCMSnExp)
    pkid <- ft$peakidx[[id]]
    cur_pks <- pks[pkid, , drop = FALSE]
    chrom <- xcms::chromatogram(xcmsExp,
                       mz = c(min(cur_pks[,2]), max(cur_pks[,3])),
                       rt = c(min(cur_pks[,5]), max(cur_pks[,6])))
    
    chromData <- lapply(chrom@.Data, function(XChrom){
        df <- as.data.frame(XChrom) 
        df$intensity[is.na(df$intensity)] <- 0
        df$sample <- fromFile(XChrom)
        return(df)
    })
    chromData <- do.call(rbind, chromData)
    fnames <- sub("\\.mz[xX]?ML$", "", basename(fileNames(XCMSnExp)))
    chromData$sample <- fnames[chromData$sample]
    print(ggplot(chromData) +
        geom_point(aes(x=rtime, y=intensity, color = as.factor(sample))) +
        # geom_area(aes(x=rtime, y=intensity, color = as.factor(sample), fill = as.factor(sample)), alpha = 0.2) +
        labs(color = "Sample")+
        scale_y_log10()+
        ggtitle(label = paste("mz:", round(min(cur_pks[,2]),4), "-", round(max(cur_pks[,3]),4)),
                subtitle = paste("PutativeID:", ft$putativeID[[id]])) +
        theme_minimal() +
        theme(legend.position="bottom") +
        guides(color = guide_legend(nrow = 3, byrow = TRUE)))
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

plotpks <- function(peakMat){
    ggplot(peakMat) + 
        geom_segment(aes(x=rtmin,xend=rtmax,y=mz,yend=mz, color = log10(maxo)))+
        geom_point(aes(x=rt,y=mz), size = 0.1, color = "tomato", alpha = 0.2)+
        facet_grid(sample~.) +
        theme_minimal()
}












