load("./Misc/rawData/rdSuper.Rdata")

#' Function to find all images and turn all pixels in the images to 0, except the 
#' locations of the cells. 
#' @param rdSuper list of RD.experiments
#' @param window region around cell to keep
imageShrinker <- function(rdSuper, window = 40){
    # For each experiment
    for( k in 1:length(rdSuper)){    
        tmpRD <- rdSuper[[k]]
        imagesToCrop <- grep("img", names(tmpRD), value=T)
        
        # Loop through each image
        for(j in 1:length(imagesToCrop)){
            # Grab the image
            img <- imagesToCrop[[j]]
            imageToPlay <- tmpRD[[ img ]]

            # Create a logic Image, to coordinate cropping
            imageLogic <- imageToPlay[,,1]
            imageLogic[,] <- FALSE

            # Define the dimension of the images
            xDim <- 1:dim(tmpRD[[ img ]])[1]
            yDim <- 1:dim(tmpRD[[ img ]])[2]

            # list to collect the list dimensions
            xLocs <- list()
            yLocs <- list()

            # These are the only cells remaining here.
            cells <- tmpRD$c.dat$id

            # Window to surround the cell
            if(dim(tmpRD[[ img ]]) > 2000){
                window <- windowSize
            }else{
                window <- windowSize / 2
            }

            # For each cell find the location surrounding the cell and update the imagelogic matrix
            for(i in 1:length(cells)){
                xLoc <- tmpRD$c.dat[cells[i], ]$center.x
                yLoc <- tmpRD$c.dat[cells[i], ]$center.y

                xLocLogic <-    ( xLoc + window ) > xDim &
                                ( xLoc - window ) < xDim 

                xLocs[[i]] <- xDim[xLocLogic]

                yLocLogic <-    ( yLoc + window ) > yDim &
                                ( yLoc - window ) < yDim 

                yLocs[[i]] <- yDim[yLocLogic]

                imageLogic[ xLocs[[i]], yLocs[[i]] ] <- TRUE
            }

            imageLogic <- t(imageLogic)
            newImage <- imageToPlay
            tmpRD[[ img ]][!imageLogic] <- 0
        }
        rdSuper[[k]] <- tmpRD
    }

    afterSize <- as.numeric(object.size(rdSuper))/1e9
    cat("\n The size after crop is ", afterSize, "GB\n")

    retunr(rdSuper)
}

# Experiment 3 is off!!!!

#current size is 664 092

# after this the file is 18,510



