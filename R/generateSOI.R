#'@export
findSOIpeaks <- function(files, DBfile, adfile){
    struct <- RHermesExp()
    struct <- setCluster(struct, BiocParallel::SerialParam())
    
    #This would be what the end user should do:
    # struct <- setDB(struct, db = DBfile, adductfile = adfile)
    
    #For testing let's just use the "demo DB" defaults 
    #(just a few formulas and adducts)
    struct <- setDB(struct)
    
    #Idea 1: Just using RHermes functions. Straightforward but with some
    #space overhead (the object can get really big, really fast if there are
    #many samples at once)
    struct <- processMS1(struct, files)
    struct <- findSOI(struct, getSOIpar(), fileID = seq_along(files))
    
    SOIs <- do.call("rbind", lapply(seq_along(files), function(soi){
        soi_list <- struct@data@SOI[[soi]]@SOIList
        soi_list$file <- files[[soi]]
        return(soi_list)
    }))
    return(SOIs)
}

