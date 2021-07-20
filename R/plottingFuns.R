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