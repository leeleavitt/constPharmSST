require(xlsx)
require(procPharm)
source('./leeProcessing/transcriptome_pharmer.r')

# Now that I have created the rdSuper I need to begin the visualization
load('./leeProcessing/rdSuper.Rdata')

# Load and Clean the transcriptome info
# The first thing we need to do is make a data base of all cells for each 
ts_info <- read.xlsx("./leeProcessing/FX SS library data 200218.xlsx",1)
ts_info <- tsInfoCleaner(ts_info)

# Load up the transcriptome data
ts_data <- read.csv('./leeProcessing/dat.out.032020.csv')
row.names(ts_data) <- make.names(ts_data$Gene.name, unique=T)

# I fucking hat the collumnNames, I need these collumn names to match the 
# ts_info$Gnomex.Label
gnomexLabs <- as.character(ts_info$Gnomex.Label)

for(i in 1:length(gnomexLabs)){
    labToRename <- grep(paste0(gnomexLabs[i], "$"), names(ts_data))
    colnames(ts_data)[labToRename] <- gnomexLabs[i]
}
names(ts_data)

ts_data <- ts_data[c('Gene.name',gnomexLabs)]
dim(ts_data)

write.csv(ts_data, 'ts_data.csv')









