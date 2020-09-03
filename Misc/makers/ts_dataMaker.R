require(procPharm)
source('./R/transcriptome_pharmer.r')

# Now that I have created the rdSuper I need to begin the visualization
#load('./leeProcessing/rdSuper.Rdata')

# Load and Clean the transcriptome info
# The first thing we need to do is make a data base of all cells for each 
ts_info <- xlsx::read.xlsx("./Misc/rawData/FX SS library data 200723.xlsx",1)
ts_info <- tsInfoCleaner(ts_info)

# Load up the transcriptome data
ts_data <- read.csv('./Misc/rawData/dat.out.090120.csv')
row.names(ts_data) <- make.names(ts_data$Gene.name, unique=T)

gnomexLabs <- as.character(ts_info$Gnomex.Label)

for(i in 1:length(gnomexLabs)){
    labToRename <- grep(paste0(gnomexLabs[i], "$"), names(ts_data))
    colnames(ts_data)[labToRename] <- gnomexLabs[i]
}
names(ts_data)

ts_data <- ts_data[c('Gene.name',gnomexLabs)]
dim(ts_data)

write.csv(ts_data, './Misc/rawData/ts_data.csv')









