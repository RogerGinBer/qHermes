soi_from_eic <- function(eic, thr = 1000, tol = 5){
    above <- which(eic > thr)
    if(length(above) == 0){return(matrix(nrow=0,ncol=2))}
    m <- diff(above) < tol
    in_soi <- F
    st <- c()
    end <- c()
    for(i in seq_along(m)){
        if(m[i] & !in_soi){
            st <- c(st, above[i])
            in_soi <- T
        } else if(!m[i] & in_soi){
            end <- c(end, above[i])
            in_soi <- F
        }
    }
    if(in_soi){end <- c(end, above[i])}
    if(length(st) == 0 | length(end) == 0){return(matrix(nrow=0,ncol=2))}
    return(cbind(scstart = st, scend = end))
}

#'@export
fastSOI <- function(MSnExp, ppm = 10, minint = 1000, DBfile = NA, 
                        adfile = NA){
    
    struct <- RHermesExp()
    struct <- setExpParam(struct, ExpParam(ppm = ppm))

    #Selecting formulas and adducts
    if(!is.na(DBfile)){
        if(!is.na(adfile)){
            struct <- setDB(struct, db = "custom", filename = DBfile,
                            adductfile = adfile)
        } else {
            struct <- setDB(struct, db = "custom", filename = DBfile)
        }
    } else {
        struct <- setDB(struct)
    }
    
    ##Calculating all possible ionic formulas
    F_DB <- struct@metadata@ExpParam@DB[,c("MolecularFormula", "EnviPatMass")]
    #Could break if colname isn't exactly "MolecularFormula"
    F_DB <- dplyr::distinct_at(F_DB, "MolecularFormula",
                                .keep_all = T)
    colnames(F_DB) <- c("fms", "m")

    message("Calculating ionic formulas")
    ionf <- RHermes:::IonicForm(F_DB, struct@metadata@ExpParam@adlist)[[1]]
    
    message("Calculating fast SOIs")
    ionf$s <- ionf$m*(1-ppm*1e-6)
    ionf$e <- ionf$m*(1+ppm*1e-6)
    
    ##Processing each file
    files <- fileNames(MSnExp)
    SOIs <- bplapply(seq_along(files),
                      function(i) {
        cur_MSnExp <- MSnbase::filterFile(MSnExp, i)
        pks <- XCHermes:::extractToHermes(cur_MSnExp)
        mzs <- unlist(pks[[3]][,1])
        ints <- unlist(pks[[3]][,2])

        ## Define the values per spectrum:
        h <- pks[[2]]
        valsPerSpect <- h$originalPeaksCount
        scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
        scantime <- h$retentionTime
        
        ## Calculate SOIs from ionic formula list
        SL <- apply(ionf, 1, function(x){
            eic <- .Call("getEIC", mzs, ints, scanindex,
                            as.double(as.numeric(c(x[5], x[6]))),
                            as.integer(c(1,length(scanindex))),
                            as.integer(length(scanindex)),
                            PACKAGE = "xcms")
            if(all(eic$intensity < minint)){return()}
            sois <- XCHermes:::soi_from_eic(eic$intensity)
            if(nrow(sois) == 0){return()}
            sois <- as.data.frame(sois)
            sois$formula <- x[1]
            sois$mz <- x[2]; sois$mzmin <- x[5]; sois$mzmax <- x[6]
            sois$peaks <- apply(sois, 1, function(row){
                data.frame(rt = scantime[row[1]:row[2]],
                           rtiv = eic$intensity[row[1]:row[2]])
            })
            sois$rtmin <- scantime[sois$scstart]
            sois$rtmax <- scantime[sois$scend]
            sois$into <- sapply(sois$peaks, function(x){sum(x[,2])})
            sois$intb <- sois$into
            sois$maxo <- sapply(sois$peaks, function(x){max(x[,2])})
            sois$rt <- sapply(sois$peaks, function(x){x[which.max(x[,2]),1]})

            return(sois)
        })
        if(is.null(SL)){return()}
        SL <- do.call("rbind", SL)
        SL$length <- SL$rtmax-SL$rtmin
        SL <- SL[SL$length > 5,] #Filter by min time
        SL$sample <- i
        
        ## Do some filtering here
        
        return(SL)
    },  BPPARAM = bpparam())
    SOIs <- do.call("rbind", SOIs)
    row.names(SOIs) <- NULL
    
    MSnExp <- as(MSnExp, "XCMSnExp")
    if (nrow(SOIs) > 0) {
        names <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
                    "intb", "maxo", "sample")
        chromPeaks(MSnExp) <- SOIs[, names] %>% as.matrix() %>%
            apply(., 2, as.numeric)
        chromPeakData(MSnExp) <- S4Vectors::DataFrame(SOIs[, 3:11],
                                                      ms_level = as.integer(1),
                                                      is_filled = FALSE)
    }
    return(MSnExp)
}

