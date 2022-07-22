#### isotopic fidelity-related ####

#' @importFrom enviPat check_chemform isowrap
#' @importFrom xcms featureDefinitions featureSpectra
#' @importFrom MetaboCoreUtils standardizeFormula
#' @export
featureFidelity <- function(XCMSnExp){
    ft <- as.data.frame(featureDefinitions(XCMSnExp))
    ft$id <- seq(nrow(ft))
    data(isotopes, package = "enviPat", envir = environment())
    fidelity <- apply(ft, 1, function(x){
        fs <- x["formula"]
        lapply(fs, function(formula, id = x[["id"]]){
            formula <- standardizeFormula(formula)
            checked <- check_chemform(isotopes, chemforms = formula)
            pattern <- isowrap(isotopes, checked = checked, resmass = F,
                               resolution = 120000, threshold = 1)[[1]]
            spec <- featureSpectra(XCMSnExp, msLevel = 1,
                                   return.type = "Spectra",
                                   features = id,
                                   method = "closest_rt")[1]
            hits <- closest(pattern[,1], mz(spec)[[1]], tolerance = 0, ppm = 3)
            query <- intensity(spec)[[1]][hits]
            target <- pattern[,2]
            if(all(is.na(query))) return(-Inf)
            query[is.na(query)] <- 0
            
            #Filter out low intensity isotopes (theoretical below thrInt)
            thrInt <- 1000
            tooLow <- (query[1] * (target / 100)) <  thrInt
            query <- query[!tooLow]
            target <- target[!tooLow]
            if(length(query) == 0) return(-Inf)
            .isofidelity(query, target)
        })
    })
    fidelity
}

.isofidelity <- function(query, pattern){
    query <- query/sum(query)
    pattern <- pattern/sum(pattern)
    probs <- mapply(function(f, p){
        .erfc(log(f / p) / (sqrt(2) * log((f + 0.03) / f)))
    }, query, pattern)
    sum(log(probs[!is.nan(probs)]))
}

.erfc <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)

#### Isotopic labelling quantification ####

#'@importFrom MSnbase filterFile
isoLabelling <- function(XCMSnExp,
                         labelled_files,
                         unlabelled_file) {
    targetFeatures <- featureDefinitions(XCMSnExp)
    targetFeatures$M0 <-
        featureValues(XCMSnExp, "maxint", value = "maxo")[, unlabelled_file[1]]
    labXCMS <- filterFile(XCMSnExp, labelled_files)
    featData <- featureData(labXCMS)@data
    featData$spectrum <- seq(nrow(featData))
    res <- apply(targetFeatures, 1, function(feat) {
        idx <- sapply(order(labelled_files), function(f) {
            cur <- featData[featData$fileIdx == f, ]
            cur$spectrum[which.min(abs(cur$retentionTime - feat[["rtmed"]]))]
        })
        spec <- spectra(labXCMS[idx])
        spec <- spec[order(idx)]
        spec <- lapply(spec, function(sp) {
            data.frame(mz = sp@mz,
                       intensity = sp@intensity)
        })
        frc <- sapply(spec, .calculateEnrichment, feature = feat)
        frc
    })
    matrix(res, ncol = length(labelled_files), byrow = TRUE)
}


.calculateEnrichment <- function(spectra, feature){
    if(is.na(feature[["M0"]])) return(0)
    maxC <- ceiling(feature[["mzmed"]]/12) #Conservative estimate of maximum carbon atoms based on mz
    match <- MetaboAnnotation::matchMz(spectra,
                              feature[["mzmed"]] + seq(0, maxC) * 1.003355,
                              MzParam())
    if(nrow(matches(match))==0) return(0)
    lab_product <- sum(query(match)[matches(match)$query_idx, "intensity"] * 
                       (matches(match)$target_idx - 1))
    if(lab_product != 0){
        # lab_product <- lab_product / feature[["M0"]]
        lab_product <- lab_product / sum(query(match)[matches(match)$query_idx, "intensity"])
    }
    lab_product
}

