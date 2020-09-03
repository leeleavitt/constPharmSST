# Here we create the rdSuper that has all rd files combine together

# The first thing we need to do is make a data base of all cells for each 
### ts_info <- read.xlsx("../rawData/FX SS library data 200218.xlsx",1) # Original
ts_info <- xlsx::read.xlsx("./Misc/rawData/FX SS library data 200723.xlsx",1) # newest
source("./R/transcriptome_pharmer.r")

# This funciton gets rid:
# rows that haven't been sequenced
# of rows that don't have an RD file
# This also gets rid of 
tsInfoClean <- tsInfoCleaner(ts_info)

# Now lets change our Directory to where the experiments are saved
mainDir <- "Y:/transcriptomePharmer/"
expDirs <- "Z:/Lee Leavitt/Cell Picking/"

setwd(expDirs)
# These are all of the unique expeirment folders to look at
uniqueExpDirs <- unique(tsInfoClean$Date.of.cDNA.prep)
# These are the unique directories in my curdir
allDirs <- list.dirs()

# Additionally we need all experiment names
rdNames <- as.character(unique(tsInfoClean$rd.name))

rdSuper <- list()
#i=12
for( i in c(1:length(uniqueExpDirs)) ){
        tryCatch({
        print(i)
        #Find out directory
        expDir <- grep(uniqueExpDirs[i], allDirs, value=T)[1]
        # go There
        setwd(expDir)
        # Load the experiment
        rdToLoad <- list.files(pattern = rdNames[i])
        exp <- get(load(rdToLoad))
 
        # Now that the experiment is loaded, we need to see which we need to keep
        cellsToAddLogic <- tsInfoClean$rd.name == rdNames[i]
        cellsToKeep <- tsInfoClean[cellsToAddLogic,][,'Cell.name',drop=T]
        cellsToKeep <- as.numeric(Reduce(c,strsplit(as.character(cellsToKeep), ',')))
        cellsToKeep <- paste0("X.", cellsToKeep)

        # Within this experiment collect the data per cell
        expTmp <- exp
        traces <- c('t.dat', 'blc')
        for( j in 1:length(traces)){
            expTmp[[ traces[j] ]] <- exp[[ traces[j] ]][,c('Time',cellsToKeep), drop=F]
        }

        data <- c('c.dat', 'bin', 'scp')
        for( j in 1:length(data)){
            expTmp[[ data[j] ]] <- exp[[ data[j] ]][cellsToKeep,,drop=F]
        }

        imgNames <- grep("img", names(expTmp), value = T)
        print(imgNames)
        expTmp  <- expTmp[ c(traces, 'w.dat', data, imgNames) ]
        nameForSuper <- sub('[.]Rdata', '', rdToLoad, ignore.case = T)
        
        rdSuper[[ nameForSuper ]] <- expTmp
    }, error = function(e){
        cat("\nCould not add experiment: ",expDir,"\n")
        }
    )
    setwd(expDirs)
}

rdSuper <- imageShrinker(rdSuper)

setwd(mainDir)

save(rdSuper, file='.Misc/rawData/rdSuper.Rdata')