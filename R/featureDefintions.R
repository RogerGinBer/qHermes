#' @export
SOIFeatureDefinitions <- function(object){
    def <- featureDefinitions(object) %>% as.data.frame()
    pk_data <- chromPeakData(object) %>% as.data.frame()
    def$formula <- lapply(def$peakidx, function(id){pk_data$formula[id[1]]})
    def$anot <- lapply(def$peakidx, function(id){pk_data$anot[id[1]]})
    def <- cbind(def, featureValues(object, value = "maxo") %>% as.data.frame())
}
