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

# #'@importFrom igraph  graph_from_adjacency_matrix simplify 
# calculate_MS2_correlation_from_feature <- function(feature_df){
#     all_ms2 <- lapply(seq(nrow(feature_df)), function(x){
#         ms2_df <- feature_df$MS2data[[x]]
#         if(is.na(ms2_df)){return()}
#         ms2_df$originalFeature <- x
#         return(ms2_df)
#     })
#     all_ms2 <- do.call(rbind, all_ms2)
#     all_ms2 <- all_ms2[sapply(all_ms2$ssdata, is.data.frame), ]
#     
#     dist <- lapply(seq(nrow(all_ms2)), function(x){
#         query <- all_ms2$ssdata[[x]]
#         sapply(seq(nrow(all_ms2)), function(y){
#             if(y >= x) return(0)
#             pattern <- all_ms2$ssdata[[y]]
#             RHermes:::calculate_MS2_distance(query, pattern, method = "cosine")
#         })
#     })
#     dist_matrix <- do.call(rbind, dist)
#     net <- graph_from_adjacency_matrix(dist_matrix, mode = "undirected")
#     net <- simplify(net)
#     netd3 <- igraph_to_networkD3(net, group = membership(cluster_fast_greedy(net)))
#     
#     newscript <- 'alert("mz: " + (d.var1.mz) + "\\n" + "Int:" + (d.var1.int));'
#     
#     fn <- forceNetwork(Links = netd3$links, Nodes = netd3$nodes, 
#                        Source = 'source', Target = 'target', 
#                        NodeID = 'name', Group = 'group',
#                        zoom = TRUE, opacity = 1, clickAction = newscript)
#     fn$x$nodes$var1 <- lapply(all_ms2$ssdata, function(ss){round(ss)})
#     fn 
#     
#     fts_cluster <- all_ms2 %>% 
#         group_by(gr = membership(cluster_fast_greedy(net))) %>%
#         summarize(size = n(), ft = list(unique(originalFeature)))
#     feature_groups <- lapply(fts_cluster$ft, function(x){feature_df[x,]})
#     feature_groups <- feature_groups[sapply(feature_groups, nrow) > 1]
#     lapply(feature_groups, function(group){
#         plots <- lapply(1:nrow(group), function(id){plot_MS2_from_feature(group, id)})
#     })
# }