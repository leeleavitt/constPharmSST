# To click in the heatmap
tsHeatMap <- function(geneDF, scale = c('none'), labels = NA, geneSelected = NA, cellSelected = NA, labelForComparison = NA){
    geneDF <- apply(geneDF, 2, rev)

    geneDF <- t(geneDF)
    par(mar=c(5,12,2,2)+0.1)
    image(geneDF, yaxt='n', xaxt='n', bty='n')

    par(xpd=T)
    # add the gene label / x labels
    geneNames <- row.names(geneDF)
    yloc <- rep(par('usr')[1] - yinch(.24), length(geneNames))
    xloc <- seq(0,1, length.out = length(geneNames))
    
    text_cex <- seq(1, .2, length.out=2000)
    #plot(xloc~exp(10*text_cex), xlab = "number of labels", ylab = "labelSize")
    text(xloc, yloc, geneNames, srt = 90, adj=1, cex = text_cex[length(geneNames)])

    # add the cell name/ cell_type label
    cellNames <- colnames(geneDF)
    xloc <- rep(par('usr')[1] - xinch(.1), length(cellNames))
    yloc <- seq(0,1, length.out = length(cellNames))
    cells_cex <- seq(1,.4, length.out=91) 
    
    if( nlevels(labels) > 2 | !is.na(labels) ){
        newPallete <- RColorBrewer::brewer.pal(n=nlevels(labels), 'Dark2')
        palette(newPallete)
        color <- rev(as.integer(labels))

        totalYDistance <- (abs(par('usr')[3]) + abs(par('usr')[4]))
        cell <- totalYDistance/length(labels)
        halfCell <- cell / 2
        
        # Add horizontal line per labelspec
        runningCellTotal <- 0
        par(xpd=F)
        for( i in length(levels(labels)):2){
            runningCellTotal <- runningCellTotal + summary(labels)[i]
            lineLocation <- (runningCellTotal * cell) - halfCell
            abline(h = lineLocation, lwd=1)
        }
        par(xpd=T)
    }else{
        color <- 'black'
        print(color)
    }

    text(
        xloc, 
        yloc, 
        cellNames, 
        adj=1, 
        cex = cells_cex[length(cellNames)], 
        col = ifelse(rev(labels) %in% labelForComparison , 'red', color) , 
        font= ifelse(rev(labels) %in% labelForComparison , 2, 1)
    )

    # add rectangular box surrouynding the selected geneName
    tryCatch(
        if( !is.na(geneSelected) ){
            geneSelected <- which(geneSelected == geneNames)
            totalXDistance <- (abs(par('usr')[1]) + abs(par('usr')[2]))
            columnWidth <- totalXDistance / length(geneNames)
            columnWidthHalf <- columnWidth / 2

            # Define top of retangle
            yBottom <- par('usr')[3]
            yTop <- par('usr')[4]
            columnMiddle <- columnWidth * (geneSelected-1)
            xLeft <- columnMiddle - columnWidthHalf
            xRight <- columnMiddle + columnWidthHalf
            rect(xLeft, yBottom, xRight, yTop, lwd=3)
        }
        , error=function(e) NULL
    )

    # add rectangular box surrouynding the selected cellName
    tryCatch(
        if( !is.na(cellSelected) ){
            totalYDistance <- (abs(par('usr')[3]) + abs(par('usr')[4]))
            rowWidth <- totalYDistance / length(cellNames)
            rowWidthHalf <- rowWidth / 2

            # Define top of retangle
            rowMiddle <- rowWidth * (length(cellNames) - cellSelected)
            yBottom <- rowMiddle - rowWidthHalf
            yTop <- rowMiddle + rowWidthHalf
            xLeft <- par('usr')[1]
            xRight <- par('usr')[2]
            rect(xLeft, yBottom, xRight, yTop, lwd=3)
        }
        , error=function(e) NULL
    )

    # Perform linear regression on genes if there are less than 100,
    tryCatch(
        if(dim(geneDF)[1] < 200 & !is.na(labelForComparison) ){
            labelsNew <- rev(labels) %in% labelForComparison
            
            geneSigLabels <- c()
            for(i in 1:dim(geneDF)[1]){
                glt <- glm(t(geneDF[i,,drop=F])~labelsNew)
                geneSigLabels[i] <- lapply(coefficients(summary(glt))[,4][-1],starfunc) #
            }

            geneSigLabels <- Reduce(c, geneSigLabels)
            
            # add the gene label / x labels
            yloc <- rep(par('usr')[2] + yinch(.24), length(geneSigLabels))
            xloc <- seq(0,1, length.out = length(geneSigLabels))
        
            text_cex <- seq(1, .2, length.out=2000)
            #plot(xloc~exp(10*text_cex), xlab = "number of labels", ylab = "labelSize")
            text(
                xloc, 
                yloc, 
                geneSigLabels, 
                srt = 90, 
                adj=0, 
                cex = text_cex[length(geneNames)],
                col = ifelse(seq(1,length(geneSigLabels)) == geneSelected, 'red', 'black')
            )
        }
        ,error= function(e) NULL
    )
}

#labels <- labelConcat
#tsHeatMap(geneDF)

# require(RColorBrewer)
# newPallete <- RColorBrewer::brewer.pal(n=nlevels(labelConcat), 'Set3')

# print(palette(newPallete))