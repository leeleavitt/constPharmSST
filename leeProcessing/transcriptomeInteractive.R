# Load up other software
source("./leeProcessing/transcriptome_pharmer.r")
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
text(0,0, "Welcome to the Transcriptome Surfer", cex=2, font=2)

# Open the heatmap window
dev.new(width = 8, height = 6)
hmWindow <- dev.cur()



dev.set(tsWindow)
# Initialize default values.
tsInfoReduce <- tsSuper$ts_info

keyPressed <- 'z'
heatMapFlag <- FALSE
renameFlag <- FALSE
cellsSelection <- NULL
genes <- c()
newGenes <- c()
libraryNames <- c()
cellRepOptions <- c("Gnomex.Label", "Experiment", "Cell.type", "rd.name", 'Cell.name')
cellRep <- c('Gnomex.Label', 'Cell.type', 'Experiment')
newLibraryNames <- c()

for(i in 1:length(cellRep)){
    newLibraryNames <- paste0(newLibraryNames,tsInfoReduce[, cellRep[i] ], '__')
}


length(genes)
while(keyPressed != 'q'){

    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(heatMapFlag){
        if(length(genes) > 0 & length(libraryNames) > 0){
            geneDF <- tsSuper$ts_data[libraryNames, genes]
            row.names(geneDF) <- newLibraryNames
            cat('\nI\'m busy plotting, my friend\n')
            dev.set(hmWindow)
            heatmap(geneDF, NA, NA)           
            dev.set(tsWindow)
            heatMapFlag <- FALSE
        }
    }
    
    # After the flags have been parsed, lets move onto which key was pressed
    keyPressed <- readkeygraph("Press q to EXIT")

    #' @param d Select cells 
    if(keyPressed == 'd'){
        # Function to return ts_info_Reduced
        # Also returns selected cell_types
        dataSelectorReturn <- dataSelector(cellsSelection = cellsSelection)
        # Unpack the return
        tsInfoReduce <- dataSelectorReturn[[1]]
        cellsSelection <- dataSelectorReturn[[2]]
        libraryNames <- row.names(tsInfoReduce)

        if(length(genes) > 0){
            heatMapFlag <- TRUE
            renameFlag <- TRUE
        }else{
            cat("\nNo genes have been selected please press g first.\n")
        }
    }

    #' @param g Select genes
    if(keyPressed == 'g'){
        toSearch <- searchSelector()
        genes <- geneFinder(toSearch)
        cat('\nThese are the genes that have come up after your search\nRemember you can press "G" to cleanup those genes that you have selected\n\n')
        cat(genes, sep=" ")
        renameFlag <- TRUE
        heatMapFlag <- TRUE
    }

    #' @param G Cleanup genes
    if(keyPressed == 'G'){
        cat('\nSelect a subset of genes from this list of genes\nRemeber you can press cancel to return the same genes\n')
        tryCatch(
            newGenes <- select.list(genes, multiple=T),
            error=function(e){
                cat('\nThere are no genes selected, make sure to press "g" first\n')
            }
        )

        # If cancel was selected return all genes
        if(is.null(newGenes)){
            newGenes <- genes
        }
        renameFlag <- TRUE
        heatMapFlag <- TRUE
    }

    #' @param h make html heatmap
    if(keyPressed == 'h'){
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
        newLibraryNames <- c()

        for(i in 1:length(cellRep)){
            newLibraryNames <- paste0(newLibraryNames,tsInfoReduce[, cellRep[i] ], '__')
        }

        }else{
            cat("\nPress d to select data to work with\n")
        }

    }

    #' @param s Save the Current geneDF as a csv to continue work outside
    if(keyPressed == 's'){
        if(!in.null(geneDF)){
            cat('\nPlease enter the name of the file\n')
            fileName <- scan(n=1,what='character', quiet=TRUE)
            fileName <- paste0("./savedFiles/",fileName,'.csv')
            write.csv(geneDF,file = fileName)
        }
    }

    #' @param ctrl-S Save the gene list you have updated 
    if(keyPressed == 'ctrl-S'){
        write.csv()

    }



}

