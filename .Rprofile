#setwd("C:/Users/leele/Documents/sst/")
# Required packages
library(utils)

## designate packages to install/load
all_pkgs <-  c('RColorBrewer', 'randomForest', 'condformat', 'graphics','grDevices', 'datasets', 'stats', 'utils', 'methods', 'base')
## find packages that need to be installed
already_installed <- rownames( installed.packages() )
to_install <- setdiff(all_pkgs, already_installed)
if (length(to_install) > 0) {
    install.packages(to_install, dependencies=TRUE)
}
## now load all packages
invisible(sapply(all_pkgs, library, character.only=TRUE))

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


