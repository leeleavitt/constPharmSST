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
tsInfoReduce <- tsSuper$ts_info[with(tsSuper$ts_info, order(Cell.type, Experiment)),,drop=F ]
cellsSelection <- as.character(tsSuper$ts_info$Cell.type)
libraryNames <- row.names(tsInfoReduce)
cellRepOptions <- c("Gnomex.Label", "Experiment", "Cell.type", "rd.name", 'Cell.name')
cellRep <- c('Gnomex.Label', 'Cell.type', 'Experiment')
newLibraryNames <- c()
logFlag <- FALSE

for(i in 1:length(cellRep)){
    newLibraryNames <- paste0(newLibraryNames,tsInfoReduce[, cellRep[i] ], '__')
}


length(genes)
while(keyPressed != 'q'){
    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(renameFlag){
        newLibraryNames <- c()
        for(i in 1:length(cellRep)){
            newLibraryNames <- paste0(newLibraryNames,tsInfoReduce[, cellRep[i] ], '__')
        }
    }

    # Update the heatmap, this also makes the matrix of 
    # genes vs cell_types
    if(heatMapFlag){
        if(length(genes) > 0 & length(libraryNames) > 0){
            geneDF <- tsSuper$ts_data[libraryNames, genes[ genes %in% newGenes ],drop=F]
            geneDF <- apply(geneDF, 2, rev)
            row.names(geneDF) <- rev(newLibraryNames)
            dev.set(hmWindow)
            if(logFlag){
                newDF <- log(geneDF)
                newDF[newDF == '-Inf'] <- 0
                heatmap(newDF, Rowv = NA, Colv = NA)
            }else{
                heatmap(geneDF, Rowv = NA, Colv = NA)
            }     
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
        newGenes <- genes
        cat('\nThese are the genes that have come up after your search\nRemember you can press "G" to cleanup those genes that you have selected\n\n')
        cat(genes, sep=" ")
        renameFlag <- TRUE
        heatMapFlag <- TRUE
    }

    #' @param G Cleanup genes
    if(keyPressed == 'G'){
        cat('\nSelect a subset of genes from this list of genes\nRemeber you can press cancel to return the same genes\n')
        tryCatch(
            newGenes <- select.list(genes, newGenes, multiple=T),
            error=function(e){
                cat('\nThere are no genes selected, make sure to press "g" first\n')
            }
        )
        # If cancel was selected return all genes
        if(length(newGenes) == 0){
            newGenes <- genes
        }
        renameFlag <- TRUE
        heatMapFlag <- TRUE
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
        geneFileName <- scan(n=1, what='character')
        geneFileName <- paste0("./searches/", geneFileName, '.txt')
        write.table(genes[ genes %in% newGenes ], file = geneFileName, sep='\n', quote=F, col.names=FALSE, row.names = FALSE)
    }
    
    #' @param normalization Save the gene list you have updated 
    if(keyPressed == 'n'){
        heatMapNormOptions <-  c("row", "column", "none", "log")
        cat('\nHow would you like to scale the heatmap?\n')
        
        newOptions <- select.list(heatMapNormOptions, multiple=F, "Normalization")
        if(newOptions != "log"){
            formals(heatmap)$scale <- newOptions
        }else{
            logFlag <- TRUE
        }
        
        heatMapFlag <- TRUE 

    }

}


newDF <- log2(geneDF)

newDF[newDF == '-Inf'] <- 0

heatmap(newDF)
