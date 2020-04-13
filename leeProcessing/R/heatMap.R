# To click in the heatmap
tsHeatMap <- function(geneDF, scale = c('none')){
    # Log transfrom if asked.

    # Scale the dataframe the way you define    
    if('log' %in% scale){
        geneDF <- log(geneDF + 1)
    }else if('row' %in% scale){
        geneDF <- scale(geneDF)
    }else if ('column' %in% scale) {
       geneDF <- scale(t(geneDF))
       geneDF <- t(geneDF)
    }

    geneDF <- t(geneDF)
    par(mar=c(5,12,2,2)+0.1)
    image(geneDF, xaxt='n', yaxt='n', bty='n')

    par(xpd=T)
    # add the gene label / x labels
    geneNames <- row.names(geneDF)
    yloc <- rep(par('usr')[1] - yinch(.1), length(geneNames))
    xloc <- seq(0,1, length.out = length(geneNames))
    text_cex <- seq(1, .2, length.out=1000)
    text(xloc, yloc, geneNames, srt = 90, adj=1, cex = text_cex[length(geneNames)])

    # add the cell name/ cell_type label
    cellNames <- colnames(geneDF)
    xloc <- rep(par('usr')[1] - xinch(.1), length(cellNames))
    yloc <- seq(0,1, length.out = length(cellNames))
    cells_cex <- seq(1,.4, length.out=91) 
    text(xloc, yloc, cellNames, adj=1, cex = cells_cex[length(cellNames)])
}
