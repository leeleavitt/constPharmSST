# First Load up the softwares
main_dir <- getwd()
rPkgs <- list.files("./R", full.names=T)
lapply(rPkgs, function(x, i){source(x, verbose=F)})

# Load up the pure data
if(length( ls(pattern = 'tsSuper') ) < 1 ){
    load("./tsSuper.Rdata")
}

profileLoader()
invisible(tsInteract(SETTINGS))

setwd(main_dir)

