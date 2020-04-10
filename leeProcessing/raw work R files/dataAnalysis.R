#install.packages('condformat')\
#load('./Misc/cell.type.blob.032020.Rdata')

#load('./LeeProcessing/tsSuper.Rdata')
source("./leeProcessing/transcriptome_pharmer.r")

# Function to Select Cells
tsInfoRed <- dataSelector(expSel = T)

libraryNames <- row.names(tsInfoRed)
newLibraryNames <- paste0(tsInfoRed$Gnomex.Label,'__',tsInfoRed$label_cellType,'__', tsInfoRed$label_experiment)

# I need to make a function that looks at the input 
# If the input doesn't return genes then us Gene.Go.Finder
#hkGenes <- as.character(tsSuper$gene.desc[tsSuper$gene.desc$Gene.stable.ID %in% hkg,'Gene.name'])

#toSearch <- c('calca', 'mrgprd','trpa1', 'trpm8', 'trpv1','kcna','p2rx', 'p2ry', 'cacna1h')
#toSearch <- c('calca', 'mrgprd','trpa1', 'trpm8', 'trpv1','pvalb', 'ntrk2', 'galnt','^th$','kcna', 'cacna1h', 'cacna1i', 'scn')

scan("./leeProcessing/searchTerms.txt", "character", sep='\n', )
#toSearch <- c('calca', 'mrgprd','trpa1', 'trpm8', 'trpv1','pvalb', 'ntrk2', '^th$','chrna', 'chrnb')
#toSearch <- c(hkGenes, toSearch)
toSearch <- searchSelector()

# Function to find genes
genes <- geneFinder(toSearch)

geneDF <- tsSuper$ts_data[libraryNames, genes]
row.names(geneDF) <- newLibraryNames

heatMapper(geneDF)


