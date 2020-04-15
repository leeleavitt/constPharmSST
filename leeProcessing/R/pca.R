#' Fucntion to select singular vectors form the percent variability
#' @param geneDF This is the gene datafram
#' @param totalSV This is the number of singular vectors to display.
singularVectorPicker <- function(geneDF, totalSV = 10){
    dev.new(width=5, height=5)
    svPicker <- dev.cur()
    geneDFSvd <- svd(geneDF+1)
    # Variance described by each singular vectors
    percVariance <- geneDFSvd$d^2/sum(geneDFSvd$d^2)*100
    plot(
        percVariance[1:totalSV],
        ylab="Percent variability explained",
        xlab = "Singular Vectors",
        type='l',
        bty='l',
        lwd = 2
    )

    points(
        seq(1,totalSV,1),
        percVariance [1:totalSV],
        pch=15,
        col='red'
    )

    par(xpd=T)
    text(
        seq(1,totalSV,1),
        percVariance[1:totalSV] + yinch(0.2),
        paste0("SV", seq(1,totalSV,1)),
        cex=.5
    )

    svSel <- identify(
        seq(1,totalSV,1),
        percVariance [1:totalSV],
        n=2,
        pch="X",
        cex = 1,
        col='red'
    )
    Sys.sleep(2)
    dev.off(svPicker)
    return(svSel)
}

#' Function to plot the BiPlot
#' @param geneDF is the matrix of cells vs genes
#' @param SV is the singular vectors to use
#' @param labels is how you want the labels to be displayed color wise
tsSVDBiPlot <- function(geneDF, SV = NA, labels = NA){
    geneDF[is.na(geneDF)]<-0
    if( nlevels(labels) > 2 | !is.na(labels) ){
        newPallete <- RColorBrewer::brewer.pal(n=nlevels(labels), 'Dark2')
        palette(newPallete)
        color <- rev(as.integer(labels))
    }else{
        color <- 'black'
    }

    if( is.na(SV) ){
        SV <- c(1,2)
    }
    geneDFSvd <- svd(geneDF)
    pc1 <- geneDFSvd$u[,SV[1], drop=F] %*% geneDFSvd$d[ SV[1] ]
    pc2 <- geneDFSvd$u[,SV[2], drop=F] %*% geneDFSvd$d[ SV[2] ]

    pd1 <- geneDFSvd$v[,SV[1],drop=F] %*% geneDFSvd$d[ SV[1] ]
    pd2 <- geneDFSvd$v[,SV[2],drop=F] %*% geneDFSvd$d[ SV[2] ]

    # Plot the principal directions first
    colorC <- hcl(0,0,0,0.5)
    par(mar=c(5,5,5,5), bty='o')
    plot(
        pd1, 
        pd2, 
        pch= '',
        axes=F,
        ylab = '', 
        xlab=''
    )

    # Add pd labels
    pd <- cbind(pd1, pd2)
    for(i in 1:dim(pd)[1]){
        segments(0,0, pd[i,1], pd[i,2], col = 'gray90')
    }
    text(pd1, pd2, colnames(geneDF), cex = .6, col = colorC, font=2)
    axis(3)
    mtext(paste0('PD', SV[1]), 3, line=3)
    axis(4)
    mtext(paste0('PD', SV[2]), 4, line=3)

    par(new=T, xpd=T)

    # Add principal component labeling
    plot(pc1, 
        pc2, 
        pch  = '', 
        axes=F,
        xlab = '',
        ylab='')
    text(pc1, pc2, row.names(geneDF), cex=.8, col=color)

    axis(1)
    mtext(paste0('PC', SV[1]), 1, line=3)

    axis(2)
    mtext(paste0('PC', SV[2]), 2, line=3)
}

