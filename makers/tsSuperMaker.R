# Making the DataSet for Transcriptome pharming
require(xlsx)
require(procPharm)
source('./leeProcessing/transcriptome_pharmer.r')

# Load up the suepr data
load("./leeProcessing/rawData/rdSuper.Rdata")
names(rdSuper)[9] <- 'RD.180613.43.m.p1'
#rdSuper <- rdSuper[-c(20,23)]

# Load and Clean the transcriptome info
# The first thing we need to do is make a data base of all cells for each 
ts_info <- read.xlsx("./leeProcessing/rawData/FX SS library data 200218.xlsx",1)
ts_info <- tsInfoCleaner(ts_info)

# load in the tsData
ts_data <- read.csv('./leeProcessing/rawData/ts_data.csv', check.names=F)
row.names(ts_data) <- make.names(ts_data$Gene.name, unique=T)
ts_data <- t(ts_data[-c(1:2)])

# What I need to do is to add the transcritome data to each experiment.
# I need to add the Gnomex Label to each cell in rdSuper

# LOOP THROUGH ALL add the gnomex id to each cells c.dat
caExperiments <- names(rdSuper)
for( i in 1:length(caExperiments)){
    expLogic <- ts_info[,'rd.name'] == caExperiments[i]
    tsInfoSel <- ts_info[expLogic, ]

    for( j in 1:dim(tsInfoSel)[1]){
        multCellsLogic <- nchar(as.character(tsInfoSel[j,]$Cell.name))
        if(multCellsLogic > 4 ){
            # Now that the experiment is loaded, we need to see which we need to keep
            cellsToKeep <- tsInfoSel[j,'Cell.name',drop=T]
            cellsToKeep <- as.numeric(Reduce(c,strsplit(as.character(cellsToKeep), ',')))
            cellsToKeep <- paste0("X.", cellsToKeep)
        }else{
            cellsToKeep <- paste0('X.', as.character(tsInfoSel[j,'Cell.name',drop=T]))
        }
        
        cellsToKeep
        rdName <- as.character(tsInfoSel[j,'rd.name'])
        rdSuper[[ rdName ]][[ 'c.dat' ]][cellsToKeep,'Gnomex.Label'] <- tsInfoSel[j,'Gnomex.Label']
    }
}

# Also Add the mouse.gene.go to the ts super
# this is packaged inside of kevins blob
load('./Misc/cell.type.blob.032020.Rdata')

# Lets get the gene info that we need for synonyms and what not
# actually very useful table
gene_info <- read.delim("leeProcessing/rawData/gene_result.txt", sep = "\t")

tsSuper <- list()
# The Calcium Imaging Experiment
tsSuper[['RD']] <- rdSuper
# The information describing each library
tsSuper[['ts_info']] <- ts_info
# The transcriptome Data
tsSuper[['ts_data']] <- ts_data
# Gene information
tsSuper[['gene_info']] <- gene_info
# Gene desriptions
tsSuper[['gene.desc']] <- mouse.gene.desc
# Gene go terms,
tsSuper[['gene.go']] <- mouse.gene.go
save(tsSuper, file='./leeProcessing/tsSuper.Rdata')
