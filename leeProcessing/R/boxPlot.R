# function for labeling significance in R
starfunc <- function(tp){
	sigc <- ""
	if(tp < .1){sigc <- "-"}
	if(tp < .05){sigc <- "*"}	
	if(tp < .01){sigc <- "**"}
	if(tp < .001){sigc <- "***"}
	return(sigc)
}

#' boxplot based on selected factors
#' @param gene the gene to view these cells
#' @param log boolean, log scale the axis
#' @param labels, this defines how the cells are grouped in each boxplot
tsBoxPlot <- function(geneDF, gene = 'Kcnc1', log = TRUE, labels = NULL){
    # Create the subset of data to work with,
    geneDF <- geneDF[, gene, drop=FALSE]
    if(log){
        #geneDF <- log(geneDF+1)
        ylab <- "log(TPM)"
    }else{
        ylab <- "TPM"
    }

    # Label Grouping, can be multiple
    if(is.null(labels)){
        labelTypes <- grep("^label", names(tsInfoReduce), value=T)
        labels <- select.list(labelTypes, multiple = T, 'Select label/s')
        labelConcat <- apply(tsInfoReduce[libraryNames, labels, drop=F], 1,paste0,collapse ='__')
        labelConcat <- factor(labelConcat, levels= unique(labelConcat))
    }else{
        labelConcat <- labels
    }

    boxLabels <- paste0(
        levels(labelConcat), 
        " : n=",
        as.character(summary(labelConcat))
    )

    # Potential xlab; paste0(gsub("^label_","",labels), collapse=" & ")
    # boxplot
    par(bty = 'l', mar = c(10, 4, 4, 1) )
    bpDims <- boxplot(
        geneDF ~ labelConcat,
        bty = 'l',
        xlab = "",
        ylab = ylab, 
        xaxs = 'n',
        names = boxLabels,
        las = 2, 
        varwidth = T,
        main = gene,
        lwd=2,
        font=2,
        boxfill = 'gray90',
        boxcol=NA,
        whisklty = 1,
        #whisklwd=2, 
        outpch=NA,
        #medlwd=2
    )

    # Added point labels
    x.jit <- jitter(as.integer(labelConcat))
    pointLabels <- as.character(tsInfoReduce[libraryNames, 'label_subClass'])
    pointLabels[is.na(pointLabels)] <- 'x'
    text(x.jit, geneDF, pointLabels, cex=.8, font=2)

    # perform kevins regression analysis
    glt <- glm(geneDF ~ labelConcat ) 
    strs <- lapply(coefficients(summary(glt))[,4],starfunc) #
    text(
        1:(length(strs)), 
        bpDims$stats[3,1:(length(strs))],
        strs,
        pos=1,
        cex=2, 
        col = 'red'
    )

    # Add gene name synonyms to the bottom left
    aliases <- as.character(tsSuper$gene_info[tsSuper$gene_info$Symbol == gene, 'Aliases'])
    xloc <- par('usr')[2]
    yloc <- par('usr')[4] + yinch(.2)
    par(xpd=T)
    text(xloc, yloc, aliases, pos=2, cex=.8)

    # Add other designations
    other_designations <- as.character(tsSuper$gene_info[tsSuper$gene_info$Symbol == gene, 'other_designations'])
    xloc <- par('usr')[1] - xinch(.7)
    yloc <- par('usr')[3] - yinch(1.5)
    par(xpd=T)
    text(xloc, yloc, other_designations, pos=4, cex=.6)

    cat("\nOther Designations\n")
    cat(strsplit(other_designations, "[|]")[[1]], sep='\n')
}