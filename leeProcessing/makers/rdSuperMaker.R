# The first thing we need to do is make a data base of all cells for each 
require(xlsx)
ts_info <- read.xlsx("./leeProcessing/rawData/FX SS library data 200218.xlsx",1)
source("./leeProcessing/transcriptome_pharmer.r")
tsInfoClean <- ts_info_cleaner(ts_info)

# Now lets change our Directory to where the experiments are saved
mainDir <- "C:/Users/leele/Documents/sst/"
expDirs <- "Z:/Lee Leavitt/Cell Picking/"

setwd(expDirs)
# These are all of the unique expeirment folders to look at
uniqueExpDirs <- unique(tsInfoClean$Date.of.cDNA.prep)
# These are the unique directories in my curdir
allDirs <- list.dirs()

# Additionally we need all experiment names
rdNames <- as.character(unique(tsInfoClean$rd.name))

rdSuper <- list()
i=12
for( i in c(1:12, 14:length(uniqueExpDirs))){
    print(i)
    #Find out directory
    expDir <- grep(uniqueExpDirs[i], allDirs, value=T)[1]
    # go There
    setwd(expDir)
    # Load the experiment
    rdToLoad <- list.files(pattern=rdNames[i])
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
    nameForSuper <- sub('[.]Rdata', '', rdToLoad)
    
    rdSuper[[ nameForSuper ]] <- expTmp
    setwd(expDirs)
}

setwd(mainDir)
save(rdSuper, file='./leeProcessing/rawData/rdSuper.Rdata')