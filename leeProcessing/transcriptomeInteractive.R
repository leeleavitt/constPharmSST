# Load up other software
rPkgs <- list.files("./leeProcessing/R", full.names=T)
lapply(rPkgs, function(x) source(x, verbose=F))

#load("./leeProcessing/tsSuper.Rdata")

readkeygraph <- function(prompt){
    getGraphicsEvent(prompt = prompt, 
                 onMouseDown = NULL, onMouseMove = NULL,
                 onMouseUp = NULL, onKeybd = onKeybd,
                 consolePrompt = "uh")
    Sys.sleep(0.01)
    return(keyPressed)
}

onKeybd <- function(key){
    keyPressed <<- key
}

#' Function to interactively change the data selected
#' 
#' 
#' 
#' 
tsInteract <- function(){}

# Open the interactive window
graphics.off()
dev.new(width = 5 , height = 5)
tsWindow <- dev.cur()
par( mar=c(0,0,0,0),xaxt='n', yaxt='n' )
plot(0,0, xlim = c(-10,10), ylim=c(-10,10), pch='')
text(0,0, "Welcome to the\nTranscriptome Surfer\nMake this window\nfocused when using\nkeyboard", cex=2, font=2)

# Open the heatmap window
dev.new(width = 10, height = 8)
hmWindow <- dev.cur()

# Open the boxplot window
dev.new(width = 9, height = 5.5)
bpWindow <- dev.cur()

# Open the peakfunc
dev.new(width = 12, height = 4)
pfWindow <- dev.cur()

# Open the biplot window
dev.new(width=7, height=7)
biPlotWindow <- dev.cur()

dev.set(tsWindow)

############################################
# Initialize default values.
tsInfoReduce <- tsSuper$ts_info

keyPressed <- 'z'
# Whether the heatmap should be plotted
heatMapFlag <- FALSE
# Renaming of the Cell_types
renameFlag <- FALSE
geneDfFlag <- FALSE
tsSvdBiPlotFlag <- FALSE
boxPlotFlag <- FALSE
labelFlag <- FALSE
cellsSelection <- NULL

# Make empty vectors for the genes and geneSubset
genes <- c()
geneSubset <- c()

# Now initialize our transcriptome information
tsInfoReduce <- tsSuper$ts_info[order(tsSuper$ts_info$label_cellType_numeric),,drop=F ]
tsInfoReduce$label_cellType <- factor(tsInfoReduce$label_cellType)

# This is a setting which is default selected cells.
cellsSelection <- as.character(tsSuper$ts_info$label_cellType)
# This is the name of the librarys, this is consistent across ts_info and ts_data
libraryNames <- row.names(tsInfoReduce)

# Options for how to represent the cells
cellRepOptions <- c("Gnomex.Label", "rd.name", 'Cell.name')
labelReps <- grep("^label_", colnames(tsInfoReduce), value=T)
cellRepOptions <- c(cellRepOptions, labelReps)
# Default represenation of the cells.
cellRep <- c('Gnomex.Label', 'label_cellType', 'label_experiment')
# Make the new cell representaion
newLibraryNames <- c()

newLibraryNames <- apply(tsInfoReduce[, cellRep ], 1, paste0, collapse="__")


# # 
# newOrder <- as.character(tsInfoReduce$Gnomex.Label)
# tsSuper$ts_data <- tsSuper$ts_data[newOrder,]

length(genes)
while(keyPressed != 'q'){
    # This makes the new labels and updates the functions with the new labels
    if(labelFlag){
        # combine these labels
        labelConcat <- apply(tsInfoReduce[libraryNames, labels, drop=F], 1,paste0,collapse ='__')
        # convert to factor, with level specified.
        labelConcat <- factor(labelConcat, levels= unique(labelConcat))
        labelConcat <- sort(labelConcat)
        libraryNames <- names(labelConcat)
        formals(tsHeatMap)$labels <- labelConcat
        formals(tsBoxPlot)$labels <- labelConcat
        formals(tsSVDBiPlot)$labels <- labelConcat

        boxplotFlag <- T
        heatMapFlag <- T
        tsSvdBiPlotFlag <- T
        renameFlag <- T
        labelFlag <- F
    }
    
    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(renameFlag){
        if(length(cellRep) > 1){
            newLibraryNames <- apply(tsInfoReduce[libraryNames, cellRep],1, paste0, collapse='__')    
        }else{
            newLibraryNames <- as.character(tsInfoReduce[libraryNames, cellRep])
        }
        geneDfFlag <- TRUE
    }
    
    # This creates the new gene data frame
    if(geneDfFlag){
        if(length(genes) > 0 & length(libraryNames) > 0){
            geneDF <- tsSuper$ts_data[libraryNames, genes[ genes %in% geneSubset ],drop=F]
            geneDF <- matrixWrangler(geneDF)
            row.names(geneDF) <- newLibraryNames
            #geneDF <- apply(geneDF, 2, rev)
            formals(tsBoxPlot)$geneDF <- geneDF
        }
        geneDfFlag <- FALSE
    }

    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(heatMapFlag){
        if(length(genes) > 0 & length(libraryNames) > 0){            
            dev.set(hmWindow)
            tsHeatMap(geneDF)    
            dev.set(tsWindow)
            heatMapFlag <- FALSE
        }
    }

    # This updates the biplot
    if(tsSvdBiPlotFlag){
        dev.set(biPlotWindow)
        tsSVDBiPlot(geneDF)

        dev.set(tsWindow)
        tsSvdBiPlotFlag <- FALSE
    }

    # This plots the boxplot
    if(boxPlotFlag){
        dev.set(bpWindow)
        tsBoxPlot(geneDF)
        dev.set(tsWindow)
        boxPlotFlag <- FALSE
    }
    
    # After the flags have been parsed, lets move onto which key was pressed
    keyPressed <- readkeygraph("Press q to EXIT")

    if(keyPressed == 'F1'){
        dev.set(hmWindow)
        clickLoc <- locator(n=1)
        
        # Find Genes
        xLoc <- clickLoc$x
        geneLocs <- seq(0,1,length.out = dim(geneDF)[2])
        geneSelected <- which.min(abs(geneLocs - xLoc))
        geneForBox <- colnames(geneDF)[geneSelected]
        formals(tsBoxPlot)$gene <- geneForBox
        dev.set(bpWindow)
        tryCatch(
            tsBoxPlot(geneDF),
                error=function(e)NULL
            )

        # Find Cell selected
        yLoc <- clickLoc$y
        cellLocs <- seq(1,0,length.out = dim(geneDF)[1])
        cellSelected <- which.min(abs(cellLocs - yLoc))
        cellForPeakunc <- libraryNames[cellSelected]
        
        rdName <- as.character(tsSuper$ts_info[cellForPeakunc,'rd.name'])
        cellName <- as.character(tsSuper$ts_info[cellForPeakunc,'Cell.name'])

        if(!nchar(cellName) > 5){
            dev.set(pfWindow)
            tryCatch(
                PeakFunc7(tsSuper$RD[[rdName]], paste0("X.", cellName), dat.n = names(tsSuper$RD[rdName])),
                error=function(e)NULL
            )
        }
        
        # Replo the heatmap with where you clicked
        dev.set(hmWindow)
        formals(tsHeatMap)$geneSelected <- geneSelected
        formals(tsHeatMap)$cellSelected <- cellSelected
        
        tsHeatMap(geneDF)
        
        dev.set(tsWindow)
    }

    #' @param d Select cells 
    if(keyPressed == 'd'){
        # Function to return ts_info_Reduced
        # Also returns selected cell_types
        dataSelectorReturn <- dataSelector(cellsSelection = cellsSelection)
        # Unpack the return
        tsInfoReduce <- dataSelectorReturn[[1]]
        cellsSelection <- dataSelectorReturn[[2]]
        libraryNames <- row.names(tsInfoReduce)
        
        #formals(tsHeatMap)$labels <- NA

        if(length(genes) > 0){
            geneDfFlag <- TRUE
            heatMapFlag <- TRUE
            renameFlag <- TRUE
            labelFlag <- TRUE
            tsSvdBiPlotFlag <- TRUE
        }else{
            cat("\nNo genes have been selected please press g first.\n")
        }
    }

    #' @param g Select genes
    if(keyPressed == 'g'){
        toSearch <- searchSelector()
        genes <- geneFinder(toSearch)
        geneSubset <- genes
        cat('\nThese are the genes that have come up after your search\nRemember you can press "G" to cleanup those genes that you have selected\n\n')
        if(length(genes) < 200){
            cat(genes, sep=" ")
        }else{
            cat(sample(genes)[1:200], sep=" ")
            cat('\n')
        }
        renameFlag <- TRUE
        heatMapFlag <- TRUE
        tsSvdBiPlotFlag <- TRUE
    }

    #' @param G Cleanup genes
    if(keyPressed == 'G'){
        cat('\nSelect a subset of genes from this list of genes\nRemeber you can press cancel to return the same genes\n')
        tryCatch(
            geneSubset <- select.list(genes, geneSubset, multiple=T),
            error=function(e){
                cat('\nThere are no genes selected, make sure to press "g" first\n')
            }
        )
        # If cancel was selected return all genes
        if(length(geneSubset) == 0){
            geneSubset <- genes
        }
        renameFlag <- TRUE
        heatMapFlag <- TRUE
        tsSvdBiPlotFlag <- TRUE
    }

    #' @param h make html heatmap
    if(keyPressed == 'h'){
        cat('\nI\'m busy plotting, my friend\n')
        geneDF <- apply(geneDF, 2, rev)
        heatMapper(geneDF)
        print(get('cf', .GlobalEnv))
    }

    #' @param r Represent Cells with Different Names
    if(keyPressed == 'r'){
        if(length(tsInfoReduce) > 0){
            cat("\nThis function displays the represenation of the cells with the collumn\nvalues from the ts_info\n")
            cellRep <- select.list(cellRepOptions, cellRep, multiple = T, "Select Cell Represenation")

            if(length(cellRep) == 0){
                cellRep <- cellRepOptions
            }
        }else{
            cat("\nPress d to select data to work with\n")
        }
        renameFlag <- TRUE
        heatMapFlag <- TRUE
        tsSvdBiPlotFlag <- TRUE
    }

    #' @param s Save the Current geneDF as a csv to continue work outside
    if(keyPressed == 's'){
        if(!is.null(geneDF)){
            cat('\nPlease enter the name of the file\n')
            fileName <- scan(n=1,what='character', quiet=TRUE)
            fileName <- paste0("./savedFiles/",fileName,'.csv')
            geneDF <- apply(geneDF, 2, rev)
            write.csv(geneDF,file = fileName)
        }
    }

    #' @param ctrl-S Save the gene list you have updated 
    if(keyPressed == 'ctrl-S'){
        geneFileName <- scan(n=1, what='character')
        geneFileName <- paste0("./searches/", geneFileName, '.txt')
        write.table(genes[ genes %in% geneSubset ], file = geneFileName, sep='\n', quote=F, col.names=FALSE, row.names = FALSE)
    }
    
    #' @param normalization Save the gene list you have updated 
    if(keyPressed == 'n'){
        matrixWranglerOptions <-  c("row", "column", "none", "log")
        cat('\nHow would you like to scale the heatmap?\n')
        newOptions <- select.list(matrixWranglerOptions, multiple=T, "Normalization")
        formals(matrixWrangler)$scale <- newOptions

        if('log' %in% newOptions){
            formals(tsBoxPlot)$log <- TRUE
        }else{
            formals(tsBoxPlot)$log <- FALSE
        }
        heatMapFlag <- TRUE 
        tsSvdBiPlotFlag <- T
        geneDfFlag <- TRUE
    }

    #' @param  labels to observe groupings for the cells
    if(keyPressed == 'l'){
        # decide what labels to work with
        labelTypes <- grep("^label", names(tsInfoReduce), value=T)
        # select the label/labels
        labels <- select.list(labelTypes, multiple = T, 'Select label/s')
        cellRep <- labels
        renameFlag <- TRUE
        if( length(labels) > 0){
            labelFlag <- T
        }else{
            formals(tsHeatMap)$labels <- NA
            formals(tsBoxPlot)$labels <- NA
            formals(tsSVDBiPlot)$labels <- NA
            heatMapFlag <- T
            tsSvdBiPlotFlag <- T
        }
    }

    #' @param Singular-vector chooser
    if(keyPressed == 'v'){
        bringToTop(-1)
        cat("Select the singular vectors that you\nwould like to observe\non the biplot\n")
        Sys.sleep(0.5)
        singularVectors <- singularVectorPicker(geneDF, 10)
        formals(tsSVDBiPlot)$SV <- singularVectors
        tsSvdBiPlotFlag <- TRUE
    }

    #' @param f random forest 
    if(keyPressed == 'f'){
        # Do a quick random forest
        if( exists('labelConcat') ){
            if(nlevels(labelConcat) > 1 ){
                # Now once we enter the randomForest, should we do all vs all, or 1 vs all?
                cat("\n1 vs All? [yes or no]\n")
                compQuestion <- select.list(c('yes','no'), title='1 vs all?')                
                # If you answer yes, now select the label for this
                if(compQuestion == 'yes' | length(compQuestion) > 0){
                    labelForComparison <- select.list(levels(labelConcat), multiple=T, title = "Select you label")
                    # This creates a factor for 1 vs all
                    rfLabel <- labelConcat == labelForComparison
                    formals(tsBoxPlot)$labelForComparison <- labelForComparison
                }else{
                    rfLabel <- labelConcat
                }
                # Walk through the forest
                rft <- randomForest::randomForest(geneDF, rfLabel)
                # Rank them my importance
                imp <- rft$importance[order(rft$importance[,1], decreasing = TRUE),]
                # Make genes this newly ranked order
                genes <- names(imp)
                
                # This is a place where i have 
                if(length(imp) > 200){
                    cat("How Many genes should I return?\n")
                    toReturn <- scan(what = 'integer', n=1)
                    importantGenes <- imp[1:toReturn]
                }else{
                    importantGenes <- imp
                }
                geneSubset <- names(importantGenes)
                heatMapFlag <- TRUE
                geneDfFlag <- TRUE
                tsSvdBiPlotFlag <- TRUE
            }else{
                cat("\nThere are not enough labels to define you cells\n")
            }
        }else{
            cat("\nYou have defined labels yet\n")
        }
    }

}



