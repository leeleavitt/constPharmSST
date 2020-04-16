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
tryCatch(
    windows(width = 5, height =2, , xpos = 0, ypos = 0),
    error=function(e)
    dev.new(width = 5 , height = 2)
)
tsWindow <- dev.cur()
par( mar=c(0,0,0,0),xaxt='n', yaxt='n' )
plot(0,0, xlim = c(-10,10), ylim=c(-10,10), pch='')
text(0,0, "Welcome to the\nTranscriptome Surfer\nMake this window\nfocused when using\nkeyboard", cex=2, font=2)

# Open the heatmap window
tryCatch(
    windows(width = 12, height = 7.8, xpos = 0, ypos = 250),
    error=function(e)
    dev.new(width = 12, height = 7.8)
)
hmWindow <- dev.cur()

# Open the boxplot window
tryCatch(
    windows(width = 9, height = 4.5, xpos = 500, ypos = 0),
    error=function(e)
    dev.new(width = 9 , height = 4.5)
)
bpWindow <- dev.cur()

# Open the peakfunc
tryCatch(
    windows(width = 12, height = 4, xpos = 760, ypos = 0),
    error=function(e)
    dev.new(width = 12 , height = 4)
)
pfWindow <- dev.cur()

# Open the biplot window
tryCatch(
    windows(width = 7.8, height = 7.8, xpos = 1150, ypos = 250),
    error=function(e)
    dev.new(width = 7 , height = 7)
)
biPlotWindow <- dev.cur()

dev.set(tsWindow)

############################################
# Initialize default values.

keyPressed <- 'z'

# Make empty vectors for the genes and geneSubset
if(is.null(SETTINGS)){
    SETTINGS <- list()
    SETTINGS[[ 'genes' ]] <- c()
    SETTINGS[[ 'geneSubset' ]] <- c()

    # Now initialize our transcriptome information
    SETTINGS[[ 'tsInfoReduce' ]] <- tsSuper$ts_info[order(tsSuper$ts_info$label_cellType_numeric),,drop=F ]
    SETTINGS[[ 'tsInfoReduce' ]]$label_cellType <- factor(SETTINGS[[ 'tsInfoReduce' ]]$label_cellType)

    # This is a setting which is default selected cells.
    SETTINGS[[ 'cellsSelection' ]] <- as.character(tsSuper$ts_info$label_cellType)
    # This is the name of the librarys, this is consistent across ts_info and ts_data
    SETTINGS[[ 'libraryNames' ]] <- row.names(SETTINGS[[ 'tsInfoReduce' ]])

    # Options for how to represent the cells
    SETTINGS[[ 'cellRepOptions' ]] <- c("Gnomex.Label", "rd.name", 'Cell.name')
    SETTINGS[[ 'labelReps' ]] <- grep("^label_", colnames(SETTINGS[[ 'tsInfoReduce' ]]), value=T)
    SETTINGS[[ 'cellRepOptions' ]] <- c(SETTINGS[[ 'cellRepOptions' ]], SETTINGS[[ 'labelReps' ]])
    # Default represenation of the cells.
    SETTINGS[[ 'cellRep' ]] <- c('Gnomex.Label', 'label_cellType', 'label_experiment')
    # Make the new cell representaion
    SETTINGS[[ 'newLibraryNames' ]] <- c()
    SETTINGS[[ 'newLibraryNames' ]] <- apply(SETTINGS[[ 'tsInfoReduce' ]][, SETTINGS[[ 'cellRep' ]] ], 1, paste0, collapse="__")
    
    #Define the label's and the label for comparison
    SETTINGS[[ 'labelConcat' ]] <- NA
    SETTINGs[[ 'labelForComparison' ]] <- NA
    
    # Define the default flag values
    SETTINGS[[ 'heatMapFlag' ]] <- FALSE
    SETTINGS[[ 'renameFlag' ]] <- FALSE
    SETTINGS[[ 'geneDfFlag' ]] <- FALSE
    SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- FALSE
    SETTINGS[[ 'boxPlotFlag' ]] <- FALSE
    SETTINGS[[ 'labelFlag' ]] <- FALSE
}

length(SETTINGS[[ 'genes' ]])
while(keyPressed != 'q'){
    # This makes the new labels and updates the functions with the new labels
    if(SETTINGS[[ 'labelFlag' ]]){
        # combine these labels
        SETTINGS[[ 'labelConcat' ]] <- apply(SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], labels, drop=F], 1,paste0,collapse ='__')
        # convert to factor, with level specified.
        SETTINGS[[ 'labelConcat' ]] <- factor(SETTINGS[[ 'labelConcat' ]], levels= unique(SETTINGS[[ 'labelConcat' ]]))
        SETTINGS[[ 'labelConcat' ]] <- sort(SETTINGS[[ 'labelConcat' ]])
        SETTINGS[[ 'libraryNames' ]] <- names(SETTINGS[[ 'labelConcat' ]])
        formals(tsHeatMap)$labels <- SETTINGS[[ 'labelConcat' ]]
        formals(tsBoxPlot)$labels <- SETTINGS[[ 'labelConcat' ]]
        formals(tsSVDBiPlot)$labels <- SETTINGS[[ 'labelConcat' ]]

        boxplotFlag <- TRUE
        SETTINGS[[ 'heatMapFlag' ]] <- TRUE
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
        SETTINGS[[ 'renameFlag' ]] <- TRUE
        SETTINGS[[ 'labelFlag' ]] <- FALSE
    }
    
    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(SETTINGS[[ 'renameFlag' ]]){
        if(length(SETTINGS[[ 'cellRep' ]]) > 1){
            SETTINGS[[ 'newLibraryNames' ]] <- apply(SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], SETTINGS[[ 'cellRep' ]]],1, paste0, collapse='__')    
        }else{
            SETTINGS[[ 'newLibraryNames' ]] <- as.character(SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], cellRep])
        }
        SETTINGS[[ 'geneDfFlag' ]] <- TRUE
    }
    
    # This creates the new gene data frame
    if(SETTINGS[[ 'geneDfFlag' ]]){
        if(length(SETTINGS[[ 'genes' ]]) > 0 & length(SETTINGS[[ 'libraryNames' ]]) > 0){
            geneDF <- tsSuper$ts_data[SETTINGS[[ 'libraryNames' ]], SETTINGS[[ 'genes' ]][ SETTINGS[[ 'genes' ]] %in% SETTINGS[[ 'geneSubset' ]] ],drop=F]
            geneDF <- matrixWrangler(geneDF)
            row.names(geneDF) <- SETTINGS[[ 'newLibraryNames' ]]
            #geneDF <- apply(geneDF, 2, rev)
            formals(tsBoxPlot)$geneDF <- geneDF
        }
        SETTINGS[[ 'geneDfFlag' ]] <- FALSE
    }

    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(SETTINGS[[ 'heatMapFlag' ]]){
        if(length(SETTINGS[[ 'genes' ]]) > 0 & length(SETTINGS[[ 'libraryNames' ]]) > 0){            
            dev.set(hmWindow)
            tsHeatMap(geneDF)    
            dev.set(tsWindow)
            SETTINGS[[ 'heatMapFlag' ]] <- FALSE
        }
    }

    # This updates the biplot
    if(SETTINGS[[ 'tsSvdBiPlotFlag' ]]){
        dev.set(biPlotWindow)
        tsSVDBiPlot(geneDF)

        dev.set(tsWindow)
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- FALSE
    }

    # This plots the boxplot
    if(SETTINGS[[ 'boxPlotFlag' ]]){
        dev.set(bpWindow)
        tsBoxPlot(geneDF)
        dev.set(tsWindow)
        SETTINGS[[ 'boxPlotFlag' ]] <- FALSE
    }
    
    ## After the flags have been parsed, lets move onto which key was pressed
    keyPressed <- readkeygraph("Press q to EXIT")

    #' @param F1 makes the heatmap interactive
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
        cellForPeakunc <- SETTINGS[[ 'libraryNames' ]][cellSelected]
        
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
        dataSelectorReturn <- dataSelector(cellsSelection = SETTINGS[[ 'cellsSelection' ]])
        # Unpack the return
        SETTINGS[[ 'tsInfoReduce' ]] <- dataSelectorReturn[[1]]
        SETTINGS[[ 'cellsSelection' ]] <- dataSelectorReturn[[2]]
        SETTINGS[[ 'libraryNames' ]] <- row.names(SETTINGS[[ 'tsInfoReduce' ]])
        
        #formals(tsHeatMap)$labels <- NA

        if(length(SETTINGS[[ 'genes' ]]) > 0){
            SETTINGS[[ 'geneDfFlag' ]] <- TRUE
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE
            SETTINGS[[ 'renameFlag' ]] <- TRUE
            SETTINGS[[ 'labelFlag' ]] <- TRUE
            SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
        }else{
            cat("\nNo genes have been selected please press g first.\n")
        }
    }
    
    #' @param f random forest 
    if(keyPressed == 'f'){
        # Do a quick and easy random forest
        if(nlevels(SETTINGS[[ 'labelConcat' ]]) > 1 ){
            # Now once we enter the randomForest, should we do all vs all, or 1 vs all?
            cat("\n1 vs All? [yes or no]\n")
            compQuestion <- select.list(c('yes','no'), title='1 vs all?')                
            # If you answer yes, now select the label for this
            if(compQuestion == 'yes' | length(compQuestion) > 0){
                SETTINGS[[ 'labelForComparison' ]] <- select.list(levels(SETTINGS[[ 'labelConcat' ]]), preSelect = SETTINGS[[ 'labelForComparison' ]], title = "Select you label")
                # This creates a factor for 1 vs all
                rfLabel <- SETTINGS[[ 'labelConcat' ]] == SETTINGS[[ 'labelForComparison' ]]
                formals(tsBoxPlot)$labelForComparison <- SETTINGS[[ 'labelForComparison' ]]
            }else{
                rfLabel <- SETTINGS[[ 'labelConcat' ]]
            }
            # Walk through the forest
            rft <- randomForest::randomForest(geneDF, rfLabel)
            # Rank them my importance
            imp <- rft$importance[order(rft$importance[,1], decreasing = TRUE),]
            # Make genes this newly ranked order
            SETTINGS[[ 'genes' ]] <- names(imp)
            
            # This is a place where i have 
            if(length(imp) > 200){
                cat("How Many genes should I return?\n")
                toReturn <- scan(what = 'integer', n=1)
                importantGenes <- imp[1:toReturn]
            }else{
                importantGenes <- imp
            }
            SETTINGS[[ 'geneSubset' ]] <- names(importantGenes)
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE
            SETTINGS[[ 'geneDfFlag' ]] <- TRUE
            SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
        }else{
            cat("\nThere are not enough labels to define you cells\nOR\nYou haven't defined labels yet\n")
        }
    }

    #' @param g Select genes
    if(keyPressed == 'g'){
        toSearch <- searchSelector()
        SETTINGS[[ 'genes' ]] <- geneFinder(toSearch)
        SETTINGS[[ 'geneSubset' ]] <- SETTINGS[[ 'genes' ]]
        cat('\nThese are the genes that have come up after your search\nRemember you can press "G" to cleanup those genes that you have selected\n\n')
        if(length(SETTINGS[[ 'genes' ]]) < 200){
            cat(SETTINGS[[ 'genes' ]], sep=" ")
        }else{
            cat(sample(SETTINGS[[ 'genes' ]])[1:200], sep=" ")
            cat('\n')
        }
        SETTINGS[[ 'renameFlag' ]] <- TRUE
        SETTINGS[[ 'heatMapFlag' ]] <- TRUE
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
    }

    #' @param G Cleanup genes
    if(keyPressed == 'G'){
        cat('\nSelect a subset of genes from this list of genes\nRemeber you can press cancel to return the same genes\n')
        tryCatch(
            SETTINGS[[ 'geneSubset' ]] <- select.list(SETTINGS[[ 'genes' ]], SETTINGS[[ 'geneSubset' ]], multiple=T),
            error=function(e){
                cat('\nThere are no genes selected, make sure to press "g" first\n')
            }
        )
        # If cancel was selected return all genes
        if(length(SETTINGS[[ 'geneSubset' ]]) == 0){
            SETTINGS[[ 'geneSubset' ]] <- SETTINGS[[ 'genes' ]]
        }
        SETTINGS[[ 'renameFlag' ]] <- TRUE
        SETTINGS[[ 'heatMapFlag' ]] <- TRUE
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
    }

    #' @param h make html heatmap
    if(keyPressed == 'h'){
        cat('\nI\'m busy plotting, my friend\n')
        geneDF <- apply(geneDF, 2, rev)
        heatMapper(geneDF)
        print(get('cf', .GlobalEnv))
    }
    
    #' @param  labels to observe groupings for the cells
    if(keyPressed == 'l'){
        # decide what labels to work with
        labelTypes <- grep("^label", names(SETTINGS[[ 'tsInfoReduce' ]]), value=T)
        # select the label/labels
        labels <- select.list(labelTypes, multiple = T, 'Select label/s')
        SETTINGS[[ 'cellRep' ]] <- labels
        SETTINGS[[ 'renameFlag' ]] <- TRUE
        if( length(labels) > 0){
            SETTINGS[[ 'labelFlag' ]] <- TRUE
        }else{
            formals(tsHeatMap)$labels <- NA
            formals(tsBoxPlot)$labels <- NA
            formals(tsSVDBiPlot)$labels <- NA
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE
            SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
        }
    }
    
    #' @param n normalization Save the gene list you have updated 
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
        SETTINGS[[ 'heatMapFlag' ]] <- TRUE 
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
        SETTINGS[[ 'geneDfFlag' ]] <- TRUE
    }

    #' @param r Represent Cells with Different Names
    if(keyPressed == 'r'){
        if(length(SETTINGS[[ 'tsInfoReduce' ]]) > 0){
            cat("\nThis function displays the represenation of the cells with the collumn\nvalues from the ts_info\n")
            SETTINGS[[ 'cellRep' ]] <- select.list(SETTINGS[[ 'cellRepOptions' ]], SETTINGS[[ 'cellRep' ]], multiple = T, "Select Cell Represenation")

            if(length(SETTINGS[[ 'cellRep' ]]) == 0){
                SETTINGS[[ 'cellRep' ]] <- SETTINGS[[ 'cellRepOptions' ]]
            }
        }else{
            cat("\nPress d to select data to work with\n")
        }
        SETTINGS[[ 'renameFlag' ]] <- TRUE
        SETTINGS[[ 'heatMapFlag' ]] <- TRUE
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
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
        write.table(SETTINGS[[ 'genes' ]][ SETTINGS[[ 'genes' ]] %in% SETTINGS[[ 'geneSubset' ]] ], file = geneFileName, sep='\n', quote=F, col.names=FALSE, row.names = FALSE)
    }
    
    #' @param Singular-vector chooser
    if(keyPressed == 'v'){
        bringToTop(-1)
        cat("Select the singular vectors that you\nwould like to observe\non the biplot\n")
        Sys.sleep(0.5)
        singularVectors <- singularVectorPicker(geneDF, 10)
        formals(tsSVDBiPlot)$SV <- singularVectors
        SETTINGS[[ 'tsSvdBiPlotFlag' ]] <- TRUE
    }

}



