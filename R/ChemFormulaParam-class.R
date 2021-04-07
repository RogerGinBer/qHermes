#' @export ChemFormulaParam
ChemFormulaParam <- setClass("ChemFormulaParam",
        slots = list(DB = "character",
                        adlist = "character",
                        ppm = "numeric",
                        ionFormulas = "data.frame"),
        prototype = list(DB = character(),
                         adlist = character(),
                         ppm = 5,
                         ionFormulas = data.frame())
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
        struct <- setDB(struct, db = "custom", filename = .Object@DB)
        struct@metadata@ExpParam@adlist <- ads[ads$adduct %in% .Object@adlist, ]
        browser()
        prepro <- RHermes:::preprocessing(struct)
        IF_DB <- prepro[[1]]
        IF_DB[[1]] <- IF_DB[[1]][order(IF_DB[[1]]$m),]
        .Object@ionFormulas <- IF_DB[[1]]
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
