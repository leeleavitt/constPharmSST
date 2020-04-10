
# boxplot based on selected factors
gene = 'Cacna1h'
log = TRUE
tsBoxPlot <- function(){
    
}

# Create the subset of data to work with,
tsData <- tsSuper$ts_data[libraryNames, gene]
if(log){
    tsData <- log(tsData+1)
}

# Label Grouping, can be multiple
labelTypes <- grep("^label", names(tsInfoReduce), value=T)

labels <- select.list(labelTypes, multiple = T, 'Select label/s')

labelConcat <- as.character(tsInfoReduce[libraryNames, labels[1]])
if(length(labels) > 1 ){
    for( i in 2:length(labels) ){
        print(i)
        labelConcat <- paste0(labelConcat, '__', as.character(tsInfoReduce[libraryNames, labels[i] ]))
    }
}

boxLabels <- paste0(
    levels(factor(labelConcat)), 
    " : n=",
    as.character(summary(as.factor(labelConcat)))
)

# boxplot
par(bty = 'l', mar = c(10, 4, 4, 1) )
bpDims <- boxplot(
    tsData ~ factor( labelConcat ),
    bty = 'l',
    xlab = '',
    xaxs = 'n',
    names = boxLabels,
    las = 2, 
    varwidth = T
)

# Added Labels
x.jit <- jitter(as.integer(as.factor(labelConcat)))
pointLabels <- as.character(tsInfoReduce[libraryNames, 'label_subClass'])
pointLabels[is.na(pointLabels)] <- '*'
text(x.jit, tsData, pointLabels, cex=.8)
identify(x.jit, tsData)


# perform kevins regression analysis
glt <- glm(tsData ~ factor(labelConcat) ) 
strs <- lapply(coefficients(summary(glt))[,4],starfunc) #
text(
    1:(length(strs)), 
    bpDims$stats[3,1:(length(strs))],
    strs,
    pos=1,
    cex=2, 
    col = 'red'
)

tsInfoReduce['15193X5',]
