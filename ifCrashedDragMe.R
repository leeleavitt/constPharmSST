setwd('..')
setwd('..')

# First Load up the softwares
main_dir <- getwd()
print(main_dir)
rPkgs <- list.files("./R", full.names=T)
invisible(lapply(rPkgs, function(x, i){source(x, verbose=F)}))

# Load up the pure data
main_dir <- getwd()
if(length( ls(pattern = 'tsSuper') ) < 1 ){
    cat('Loading the data, please wait..\n\n')
    flush.console()
    load("./tsSuper.Rdata")
    alarm()
}

profileLoader()
