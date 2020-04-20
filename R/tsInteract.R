#' Function to interactively change the data selected
#' @param F1 makes the heatmap interactive
#' @param d Select cells 
#' @param f random forest 
#' @param g Select genes
#' @param G Cleanup selected genes
#' @param h make html heatmap
#' @param l labels to observe groupings for the cells
#' @param n normalization Save the gene list you have updated 
#' @param r Represent Cells with Different Names
#' @param s Save the Current geneDF as a csv to continue work outside
#' @param ctrl-S Save the gene list you have updated 
#' @param v Singular-vector chooser
tsInteract <- function(SETTINGS){
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
        formals(tsBoxPlot)$SETTINGS <- SETTINGS

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
        SETTINGS[[ 'labelForComparison' ]] <- ''
        
        # Define the default flag values
        SETTINGS[[ 'heatMapFlag' ]] <- FALSE
        SETTINGS[[ 'renameFlag' ]] <- FALSE
        SETTINGS[[ 'geneDfFlag' ]] <- FALSE
        SETTINGS[[ 'biPlotFlag' ]] <- FALSE
        SETTINGS[[ 'boxPlotFlag' ]] <- FALSE
        SETTINGS[[ 'labelFlag' ]] <- FALSE

        # Returned biplot dimensions
        SETTINGS[[ 'biPlotDims' ]] <- NA
        
        # Scaling for geneDF
        SETTINGS[[ 'newOptions' ]] <- 'log'
    }else{
        formals(matrixWrangler)$scale <- SETTINGS[[ 'newOptions' ]]
        SETTINGS[[ 'heatMapFlag' ]] <- TRUE
        SETTINGS[[ 'renameFlag' ]] <- TRUE
        SETTINGS[[ 'geneDfFlag' ]] <- TRUE
        SETTINGS[[ 'biPlotFlag' ]] <- TRUE
        SETTINGS[[ 'labelFlag' ]] <- TRUE
    }

    while(keyPressed != 'q'){
        # This makes the new labels and updates the functions with the new labels
        if(SETTINGS[[ 'labelFlag' ]]){
            # combine these labels
            SETTINGS[[ 'labelConcat' ]] <- apply(
                SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], SETTINGS[[ 'cellRep' ]], drop=F], 
                1,
                paste0,
                collapse ='__'
            )
            # convert to factor, with level specified.
            SETTINGS[[ 'labelConcat' ]] <- factor(SETTINGS[[ 'labelConcat' ]], levels= unique(SETTINGS[[ 'labelConcat' ]]))
            SETTINGS[[ 'labelConcat' ]] <- sort(SETTINGS[[ 'labelConcat' ]])
            SETTINGS[[ 'libraryNames' ]] <- names(SETTINGS[[ 'labelConcat' ]])
            formals(tsHeatMap)$labels <- SETTINGS[[ 'labelConcat' ]]
            formals(tsBoxPlot)$labels <- SETTINGS[[ 'labelConcat' ]]
            formals(tsSVDBiPlot)$labels <- SETTINGS[[ 'labelConcat' ]]

            SETTINGS[[ 'boxplotFlag' ]] <- TRUE
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
            SETTINGS[[ 'renameFlag' ]] <- TRUE
            SETTINGS[[ 'labelFlag' ]] <- FALSE
        }
        
        # Update the heatmap, this also makes the matrix of 
        # genes vs cell_types
        if(SETTINGS[[ 'renameFlag' ]]){
            if(length(SETTINGS[[ 'cellRep' ]]) > 1){
                SETTINGS[[ 'newLibraryNames' ]] <- apply(SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], SETTINGS[[ 'cellRep' ]]],1, paste0, collapse='__')    
            }else{
                SETTINGS[[ 'newLibraryNames' ]] <- as.character(SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], SETTINGS[[ 'cellRep' ]]])
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
        if(SETTINGS[[ 'biPlotFlag' ]]){
            dev.set(biPlotWindow)
            tryCatch(
                SETTINGS[[ 'biPlotDims' ]] <- tsSVDBiPlot(geneDF),
                error=function(e)NULL
            )

            dev.set(tsWindow)
            SETTINGS[[ 'biPlotFlag' ]] <- FALSE
        }

        # This plots the boxplot
        if(SETTINGS[[ 'boxPlotFlag' ]]){
            dev.set(bpWindow)
            formals(tsBoxPlot)$SETTINGS <- SETTINGS
            tryCatch(
                tsBoxPlot(geneDF),
                error = function(e) NULL
            )
            dev.set(tsWindow)
            SETTINGS[[ 'boxPlotFlag' ]] <- FALSE
        }
        
        ## After the flags have been parsed, lets move onto detecting the keypress
        keyPressed <- readkeygraph("Press q to EXIT")

        #' @param F1 makes the heatmap interactive
        if(keyPressed == 'F1'){
            dev.set(hmWindow)
            clickLoc <- locator(n=1)
            
            # Find Genes
            xLoc <- clickLoc$x
            geneLocs <- seq(0,1,length.out = dim(geneDF)[2])
            geneSelected <- which.min(abs(geneLocs - xLoc))
            
            # Update the boxplot
            geneForBox <- colnames(geneDF)[geneSelected]
            formals(tsBoxPlot)$gene <- geneForBox
            SETTINGS[[ 'boxPlotFlag' ]] <- TRUE

            # Update the biPlot
            formals(tsSVDBiPlot)$geneSelected <- geneForBox
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE

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
            }else{
                dev.set(pfWindow)
                tryCatch(
                    {
                    cellName <- as.numeric(Reduce(c,strsplit(as.character(cellName), ',')))
                    cellName <- paste0("X.", cellName)
                    print(cellName)
                    LinesEvery.6(tsSuper$RD[[rdName]], cellName, dat.n = names(tsSuper$RD[rdName]))
                    }
                    ,error=function(e)NULL
                )

            }
            
            # Replo the heatmap with where you clicked
            dev.set(hmWindow)
            formals(tsHeatMap)$geneSelected <- geneForBox
            formals(tsHeatMap)$cellSelected <- cellSelected
            
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE
            
            dev.set(tsWindow)
        }

        #' @param F2 makes the biPlot interactive
        if(keyPressed == 'F2'){
            dev.set(biPlotWindow)
            clickLoc <- identify(SETTINGS[[ 'biPlotDims' ]],n=1, labels='' )
            geneSelected <- row.names(SETTINGS[[ 'biPlotDims' ]])[clickLoc]
            formals(tsSVDBiPlot)$geneSelected <- geneSelected
            SETTINGS[[ 'biPlotDims' ]] <- tsSVDBiPlot(geneDF)

            # Update the boxplot with this selection
            formals(tsBoxPlot)$gene <- geneSelected
            SETTINGS[[ 'boxPlotFlag' ]] <- TRUE

            # update the heatMap
            formals(tsHeatMap)$geneSelected <- geneSelected
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE

            dev.set(tsWindow)
        }
        

        if(keyPressed == 'b'){
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
        }

        #' @param d Select cells 
        if(keyPressed == 'c'){
            # Function to return ts_info_Reduced
            # Also returns selected cell_types
            dataSelectorReturn <- dataSelector(cellsSelection = SETTINGS[[ 'cellsSelection' ]])
            # Unpack the return
            SETTINGS[[ 'tsInfoReduce' ]] <- dataSelectorReturn[[1]]
            SETTINGS[[ 'cellsSelection' ]] <- dataSelectorReturn[[2]]
            SETTINGS[[ 'libraryNames' ]] <- row.names(SETTINGS[[ 'tsInfoReduce' ]])
            
            formals(tsBoxPlot)$SETTINGS <- SETTINGS

            formals(tsBoxPlot)$labelForComparison <- NA
            formals(tsHeatMap)$labelForComparison <- NA


            if(length(SETTINGS[[ 'genes' ]]) > 0){
                SETTINGS[[ 'geneDfFlag' ]] <- TRUE
                SETTINGS[[ 'heatMapFlag' ]] <- TRUE
                SETTINGS[[ 'renameFlag' ]] <- TRUE
                SETTINGS[[ 'labelFlag' ]] <- TRUE
                SETTINGS[[ 'biPlotFlag' ]] <- TRUE
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
                if(compQuestion == 'yes' | compQuestion != ''){
                    SETTINGS[[ 'labelForComparison' ]] <- select.list(
                            levels(SETTINGS[[ 'labelConcat' ]]), 
                            multiple = TRUE, 
                            preselect = SETTINGS[[ 'labelForComparison' ]], 
                            title = "Select your label"
                        )
                    rfLabel <- SETTINGS[[ 'labelConcat' ]] %in% SETTINGS[[ 'labelForComparison' ]]
                    formals(tsBoxPlot)$labelForComparison <- SETTINGS[[ 'labelForComparison' ]]
                    formals(tsHeatMap)$labelForComparison <- SETTINGS[[ 'labelForComparison' ]] 
                }else{
                    rfLabel <- SETTINGS[[ 'labelConcat' ]]
                }

                # This is a place where i have 
                if(dim(geneDF)[2] > 200){
                    bringToTop(-1)
                    alarm()
                    cat("How Many genes should I return?\n")
                    toReturn <- scan(what = 'integer', n=1)
                    if(length(toReturn) == 0){
                        toReturn <- 199
                    }
                }else{
                    toReturn <- dim(geneDF)[2]
                }

                # Walk through the forest
                rft <- randomForest::randomForest(geneDF, rfLabel)
                # Rank them my importance
                imp <- rft$importance[order(rft$importance[,1], decreasing = TRUE),]
                importantGenes <- imp[1:toReturn]

                # Make genes this newly ranked order
                SETTINGS[[ 'genes' ]] <- names(imp)                
                SETTINGS[[ 'geneSubset' ]] <- names(importantGenes)

                # Turn on all flags
                SETTINGS[[ 'heatMapFlag' ]] <- TRUE
                SETTINGS[[ 'geneDfFlag' ]] <- TRUE
                SETTINGS[[ 'biPlotFlag' ]] <- TRUE
                SETTINGS[[ 'boxPlotFlag' ]] <- TRUE
            }else{
                cat("\nThere are not enough labels to define you cells\nOR\nYou haven't defined labels yet\n")
            }
        }

        #' @param g Select genes
        if(keyPressed == 'g'){
            toSearch <- searchSelector()
            SETTINGS[[ 'genes' ]] <- unique(geneFinder(toSearch))
            SETTINGS[[ 'geneSubset' ]] <- SETTINGS[[ 'genes' ]]
            cat("\nYour search returned\n\n", length(SETTINGS[[ 'genes' ]]),' Genes\n')
            cat('\nThese are the genes that have come up after your search\nRemember you can press "G" to cleanup those genes that you have selected\n\n')
            if(length(SETTINGS[[ 'genes' ]]) < 50){
                cat(SETTINGS[[ 'genes' ]], sep=" ")
            }else{
                cat(sample(SETTINGS[[ 'genes' ]])[1:50], sep=" ")
                cat('\n')
            }
            SETTINGS[[ 'renameFlag' ]] <- TRUE
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
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
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
        }

        #' @param h make html heatmap
        if(keyPressed == 'h'){
            cat('\nI\'m busy plotting, my friend\n')
            geneDF <- apply(geneDF, 2, rev)
            heatMapper(geneDF)
            print(get('cf', .GlobalEnv))
        }
        
        #' @param l labels to observe groupings for the cells
        if(keyPressed == 'l'){
            # decide what labels to work with
            labelTypes <- grep("^label", names(SETTINGS[[ 'tsInfoReduce' ]]), value=T)
            # select the label/labels
            SETTINGS[[ 'cellRep' ]] <- select.list(labelTypes, multiple = T, 'Select label/s')
            SETTINGS[[ 'renameFlag' ]] <- TRUE
            if( length(SETTINGS[[ 'cellRep' ]]) > 0){
                SETTINGS[[ 'labelFlag' ]] <- TRUE
            }else{
                formals(tsHeatMap)$labels <- NA
                formals(tsBoxPlot)$labels <- NA
                formals(tsSVDBiPlot)$labels <- NA
                SETTINGS[[ 'heatMapFlag' ]] <- TRUE
                SETTINGS[[ 'biPlotFlag' ]] <- TRUE
            }
        }
        
        #' @param n normalization Save the gene list you have updated 
        if(keyPressed == 'n'){
            matrixWranglerOptions <-  c("row", "column", "none", "log")
            cat('\nHow would you like to scale the heatmap?\n')
            SETTINGS[[ 'newOptions' ]] <- select.list(matrixWranglerOptions, multiple=T, "Normalization")
            formals(matrixWrangler)$scale <- SETTINGS[[ 'newOptions' ]]

            if('log' %in% SETTINGS[[ 'newOptions' ]]){
                formals(tsBoxPlot)$log <- TRUE
            }else{
                formals(tsBoxPlot)$log <- FALSE
            }
            SETTINGS[[ 'heatMapFlag' ]] <- TRUE 
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
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
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
        }

        #' @param s Save the Current geneDF as a csv to continue work outside
        if(keyPressed == 's'){
            if(!is.null(geneDF)){
                cat('\nPlease enter the name of the file\n')
                bringToTop(-1)
                fileName <- scan(n=1,what='character', quiet=TRUE)
                fileName <- paste0("./savedCsv/",fileName,'.csv')
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
        
        #' @param v Singular-vector chooser
        if(keyPressed == 'v'){
            bringToTop(-1)
            cat("Select the singular vectors that you\nwould like to observe\non the biplot\n")
            Sys.sleep(0.5)
            singularVectors <- singularVectorPicker(geneDF, 10)
            formals(tsSVDBiPlot)$SV <- singularVectors
            SETTINGS[[ 'biPlotFlag' ]] <- TRUE
        }

        if(keyPressed == 'q'){
            save(SETTINGS, file = "SETTINGS.Rdata")
            graphics.off()
            cat("\nWould you like to rename your profile?\n")
            
            renameLogic <- select.list(c('yes', 'no'), title = "Rename Profile?")
            if(renameLogic == 'yes'){
                oldFolder <- paste0("./profiles/", rev(strsplit(getwd(), '/')[[1]])[1])
                # Go back two directories
                setwd("..")
                setwd("..")
                
                # Get the new profile name
                cat("\n!!!!!CLOSE ALL FILES AND FOLDERS!!!\n")
                alarm()
                Sys.sleep(.5)
                alarm()
                cat('\nEnter the name of your profile, buddy\n')
                newFolder <- scan(what = 'character', n=1, quiet = T, sep=">")
                
                # New folder
                newFolder <- paste0("./profiles/", newFolder)
                dir.create(newFolder)
                # Old folder
                oldFiles <- list.dirs(oldFolder, recursive=TRUE)[-1]
                
                # Make the savedCsv folder
                file.copy(
                    oldFiles,
                    newFolder,
                    recursive=TRUE
                )
                
                # Now save the Settings to the correct folder
                setwd(newFolder)
                save(SETTINGS, file = "SETTINGS.Rdata")
                setwd("..")
                setwd("..")

                # When done, delete the older folder!
                unlink(oldFolder, TRUE, TRUE)
            }else{
                setwd("..")
                setwd("..")
            }
            cat('\nWould you like to keep analyzing?\n')
            continueLogic <- select.list(c('yes', 'no'), title='Continue Analysis?')
            if(continueLogic == 'yes'){
                profileLoader()
            }
        }
    }
}

