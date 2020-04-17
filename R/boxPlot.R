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
tsBoxPlot <- function(geneDF, gene = 'Kcnc1', log = TRUE, labels = NULL, labelForComparison = NA, SETTINGS){
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
        labelTypes <- grep("^label", names(SETTINGS[[ 'tsInfoReduce' ]]), value=T)
        labels <- select.list(labelTypes, multiple = T, 'Select label/s')
        labelConcat <- apply(SETTINGS[[ 'tsInfoReduce' ]][SETTINGS[[ 'libraryNames' ]], labels, drop=F], 1,paste0,collapse ='__')
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
        axes=FALSE,
        #names = '',
        las = 2, 
        varwidth = T,
        #main = gene,
        lwd=2,
        font=2,
        boxfill = 'gray90',
        boxcol=NA,
        whisklty = 1,
        #whisklwd=2, 
        outpch=NA,
        #medlwd=2
    )
    
    # Give the Plot a stylish Name top left corner
    par(xpd=T)
    text(
        par('usr')[1],
        par('usr')[4] + yinch(.5),
        gene,
        font = 2,
        cex = 2, 
        adj = 0
    )

    # Label under each boxplot
    box()
    axis(2)
    axis(1, seq(1,length(boxLabels)), labels = NA)

    text(
        seq(1,length(boxLabels)),
        rep(par('usr')[3]-yinch(.2), length(boxLabels)),
        boxLabels,  
        srt=90,
        col = ifelse(levels(labelConcat) == labelForComparison, 'red', 'black'),
        font = ifelse(levels(labelConcat) == labelForComparison, 2 , 1),
        adj = 1
    )

    # Added point labels
    x.jit <- jitter(as.integer(labelConcat))
    pointLabels <- as.character(SETTINGS[[ "tsInfoReduce" ]][SETTINGS[[ 'libraryNames' ]], 'label_subClass'])
    pointLabels[is.na(pointLabels)] <- 'x'
    text(x.jit, geneDF, pointLabels, cex=.8, font=2)

    # perform kevins regression analysis
    if(!is.na(labelForComparison)){
        toNotLabel <- which(levels(labelConcat) == labelForComparison)
        labelConcat <- relevel(labelConcat, labelForComparison)
        glt <- glm( geneDF ~ labelConcat ) 
        strs <- lapply(coefficients(summary(glt))[,4][-1],starfunc) #
        
        strSeq <- seq(1, length(levels(labelConcat))) 
        strSeq <- strSeq[strSeq!=toNotLabel]
        text(
            strSeq, 
            bpDims$stats[3,strSeq],
            strs,
            pos=1,
            cex=2, 
            col = 'red'
        )
    }

    # Add gene name synonyms to the top right
    aliases <- as.character(tsSuper$gene_info[tsSuper$gene_info$Symbol == gene, 'Aliases'])
    if(length(aliases) == 0){
        aliases <- 'NA'
    }
    aliases <- paste0('Aliases: ', aliases)
    xloc <- par('usr')[2]
    yloc <- par('usr')[4] + yinch(.3)
    par(xpd=T)
    text(xloc, yloc, aliases, pos=2, cex=.7)

    # Add other designations
    other_designations <- as.character(tsSuper$gene_info[tsSuper$gene_info$Symbol == gene, 'other_designations'])
    
    if(length(other_designations) == 0){
        other_designations <- 'NA'
    }else{
        other_designations <- paste0('Other Designations: |', other_designations)
        other_designations <- gsub('[|]','\n',other_designations)
    }
    xloc <- par('usr')[1] - xinch(.7)
    yloc <- par('usr')[3] - yinch(.8)
    par(xpd=T)
    text(xloc, yloc, other_designations, pos=4, cex=.7)

    bringToTop(-1)
    cat('\n\n\n\n')
    cat(strsplit(other_designations, "[|]")[[1]], sep='\n')
    cat('\n\n\n\n')

}