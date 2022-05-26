#' @title ChemFormulaParam
#' @description Main ionic formula structure
#' @slot DB Character, path to the molecular formula database
#' @slot adlist Character, list of adducts that will be used to generate ionic
#'   formulas
#' @slot ppm Numeric, instrumental mass error (in parts per million, ppm).
#' @slot ionFormulas Dataframe. NOT to be manually provided by the user. It is
#'   automatically filled when the object is created, using the provided DB and
#'   adlist.
#' @export ChemFormulaParam
ChemFormulaParam <- setClass("ChemFormulaParam",
        slots = list(DB = "character",
                        adlist = "character",
                        ppm = "numeric",
                        ionFormulas = "data.frame",
                        mode = "character"),
        prototype = list(DB = character(),
                         adlist = character(),
                         ppm = 5,
                         ionFormulas = data.frame(),
                         mode = "RHermes")
        )


setMethod("initialize", signature = "ChemFormulaParam",
    function(.Object, ...){
        .Object <- callNextMethod()
        validObject(.Object)
        if(length(.Object@adlist) == 0){
            message(paste("No adducts provided, proceeding with M+H, M+Na,",
                          "M+K and M+NH4"))
            .Object@adlist <- c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+")
        }
        struct <- RHermesExp()
        ads <- adducts_MetaboCore_to_Hermes()
        struct <- setDB(struct, db = "custom", DBfile = .Object@DB)
        
        #Retrieve molecular formula database
        DB <- database_importer("custom",
                                    filename = .Object@DB,
                                    minmass = -Inf,
                                    maxmass = Inf,
                                    filter = FALSE)

        #Join all positive and negative ads
        ad <- adductTables(3,3)
        ion <- struct@metadata@ExpParam@ion
        struct@metadata@ExpParam@adlist <- rbind(ad[[1]], ad[[2]])
        row.names(struct@metadata@ExpParam@adlist) <- NULL
        
        #Filter adducts with given list
        struct@metadata@ExpParam@adlist <- ads[ads$adduct %in% .Object@adlist, ]
        
        F_DB <- DB[,c("MolecularFormula", "EnviPatMass")]
        #Could break if colname isn't exactly "MolecularFormula"
        F_DB <- dplyr::distinct_at(F_DB, "MolecularFormula",
                                   .keep_all = T)
        colnames(F_DB) <- c("fms", "m")
        
        IF_DB <- IonicForm(F_DB, struct@metadata@ExpParam@adlist)
    
        IF_DB <- IF_DB[order(IF_DB$m),]
        .Object@ionFormulas <- IF_DB
        .Object
})


#' @importFrom MetaboCoreUtils adducts
adducts_MetaboCore_to_Hermes <- function(){
    ads <- rbind(MetaboCoreUtils::adducts("positive"),
                MetaboCoreUtils::adducts("negative"))
    
    #Reorder columns
    ads <- ads[,c("name", "charge", "mass_multi", "mass_add", "positive",
                  "formula_add", "formula_sub")]
    
    #Adapt names and format
    names(ads) <- c("adduct", "Charge", "Mult", "massdiff", "Ion_mode",
                    "formula_add", "formula_ded")
    ads$Ion_mode <- ifelse(ads$Ion_mode, "positive", "negative")
    ads$Mult <- ads$Mult * abs(ads$Charge)
    ads$massdiff <- ads$massdiff * abs(ads$Charge)

    #Adapt atoms to add and subtract
    incomplete <- grepl("[A-Za-z]$", ads$formula_add)
    ads$formula_add[incomplete] <- paste0(ads$formula_add[incomplete], "1")
    
    incomplete <- grepl("[A-Za-z]$", ads$formula_ded)
    ads$formula_ded[incomplete] <- paste0(ads$formula_ded[incomplete], "1")
    ads$formula_ded[ads$formula_ded == ""] <- "C0"
    
    ads
}


adductTables <- function(ch_max = 1, mult_max = 1) {
    data(adducts, package = "enviPat", envir = environment())
    adducts$Mass[49] <- adducts$Mass[49] * (-1)  #Fixed wrong one
    
    negative.envi <- adducts[which(adducts$Ion_mode == "negative"), ]
    positive.envi <- adducts[which(adducts$Ion_mode == "positive"), ]
    
    ad_positive <- which(
        positive.envi$Charge %in% as.character(seq(ch_max)) &
            positive.envi$Mult %in% seq(mult_max)
    )
    ad_negative <- which(
        negative.envi$Charge %in% as.character(seq(from = -1, to = -ch_max)) &
            negative.envi$Mult %in% seq(mult_max)
    )
    
    positive.envi <- positive.envi[ad_positive, ]
    negative.envi <- negative.envi[ad_negative, ]
    
    negative.ad <- negative.envi[, -c(2, 9)]
    colnames(negative.ad)[c(1, 4)] <- c("adduct", "massdiff")
    negative.ad$massdiff <- negative.ad$massdiff * abs(negative.ad$Charge)
    
    positive.ad <- positive.envi[, -c(2, 9)]
    colnames(positive.ad)[c(1, 4)] <- c("adduct", "massdiff")
    positive.ad$massdiff <- positive.ad$massdiff * abs(positive.ad$Charge)
    
    return(list(negative.ad, positive.ad))
}

#' @importFrom utils data read.csv read.csv2 write.csv
database_importer <- function(template = "hmdb",
                              filename = "./app/www/norman.xls",
                              minmass = 70, maxmass = 750, keggpath = "",
                              filter = TRUE) {
    db <- import_DB(template, filename, keggpath)
    if (filter) db <- filter_DB(db)
    db <- calculate_envipat_mass(db)
    db <- filter(db, dplyr::between(as.numeric(db$EnviPatMass),
                                    minmass, maxmass))
    return(db)
}

import_DB <- function(template, filename, keggpath){
    if (template == "hmdb") {
        db <- read.csv(system.file("extdata", "hmdb.csv", package = "RHermes"))
        db <- db[!grepl("\\.", db$MolecularFormula), ]
    } else if (template == "norman") {
        db <- readxl::read_excel(system.file("extdata", "norman.xlsx",
                                             package = "RHermes"))
        names(db)[names(db) == "MOLECULAR_FORMULA"] <- "MolecularFormula"
        names(db)[names(db) == "NAME"] <- "Name"
    } else if (template == "custom") {
        if (grepl(pattern = "csv", x = filename)) {
            db <- read.csv(filename)
            if(ncol(db) == 1){
                db <- read.csv2(filename)
            }
        } else {
            db <- readxl::read_excel(filename)
        }
        if (!all(c("MolecularFormula", "Name") %in% names(db))) {
            stop("The colnames aren't adequate. The file must contain the
            'MolecularFormula' and 'Name' columns at least")
        }
    } else if (template == "kegg_p") {
        suppressWarnings(
            splitpath <- split(seq_along(keggpath),
                               seq_len(ceiling(length(keggpath)/10)))
        )
        comp <- lapply(splitpath, function(x) {
            pathdata <- KEGGREST::keggGet(keggpath[x])
            lapply(pathdata, function(y) {
                names(y$COMPOUND)
            })
        })
        comp <- unique(unlist(comp))
        suppressWarnings(splitcomp <- split(seq_along(comp),
                                            seq_len(ceiling(length(comp)/10))))
        db <- lapply(splitcomp, function(x) {
            data <- KEGGREST::keggGet(comp[x])
            data <- lapply(data, function(y) {
                c(y[1], gsub(y[2][[1]][[1]], pattern = ";", replacement = ""),
                  y[3])
            })
            as.data.frame(do.call(rbind, data))
        })
        db <- as.data.frame(do.call(rbind, db))
        names(db)[c(2, 3)] <- c("Name", "MolecularFormula")
    } else {
        stop("You haven't entered a valid template. Options are: 'hmdb',
        'norman' and 'custom'")
    }
    return(db)
}

filter_DB <- function(db){
    #Clean molecular formulas that contain unknown elements
    db <- db[grepl("^C.?.?", db$MolecularFormula), ]
    db <- db[!grepl("[", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("R", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("X", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("T", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(".", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(")", db$MolecularFormula, fixed = TRUE), ]
    db$MolecularFormula <- as.character(db$MolecularFormula)
    return(db)
}

calculate_envipat_mass <- function(db){
    #Calculate M0 mass from formula
    isotopes <- NULL #To appease R CMD Check "no visible binding"
    data(isotopes, package = "enviPat", envir = environment())
    db$EnviPatMass <- lapply(db$MolecularFormula, function(x) {
        res <- tryCatch(enviPat::isopattern(isotopes, x, threshold = 99,
                                            verbose = FALSE),
                        error = function(cond){NA})
        # Isotopes should be loaded first
        if (is.matrix(res[[1]])) {
            res <- res[[1]][1, 1]
            res <- as.numeric(unname(res))
            return(res)
        } else {
            return(NA)
        }
    })
    if (any(is.na(db$EnviPatMass))) {
        db <- db[!is.na(db$EnviPatMass), ]
    }
    db$EnviPatMass <- as.numeric(db$EnviPatMass)
    return(db)
}


IonicForm <- function(F_DB, Ad_DB) {
    RES <- lapply(seq_len(nrow(F_DB)), calculate_ionic_forms,
                  F_DB = F_DB,
                  Ad_DB = Ad_DB)
    db <- do.call(rbind, RES)
    db <- db[!duplicated(db[, 1]), ]
    return(db)
}


calculate_ionic_forms <- function(i, F_DB, Ad_DB){
    f <- as.character(F_DB$fms[i])
    j <- apply(Ad_DB, 1, function(x) {
        current_f <- paste(f, x[1], sep = "_")
        return(current_f)
    })
    good <- which(!is.na(j))
    j <- j[!is.na(j)]
    if(length(j) == 0){return()}
    
    charge <- abs(as.numeric(Ad_DB[good, 2]))
    multiplicity <- as.numeric(Ad_DB[good,3])
    adduct_delta <- as.numeric(Ad_DB[good,4])
    
    mass <- ((F_DB$m[i] * multiplicity + adduct_delta) / charge) %>%
        round(., digits = 5)
    db <- data.table(f = j, m = mass, ch = charge, envi = j)
    return(db)
}
