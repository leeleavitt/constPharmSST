# #toSearch <- searchSelector()
# SETTINGS[[ 'genes' ]] <- unique(geneFinder(toSearch))
# SETTINGS[[ 'geneSubset' ]] <- SETTINGS[[ 'genes' ]]
# cat("\nYour search returned\n\n", length(SETTINGS[[ 'genes' ]]),' Genes\n')
# cat('\nThese are the genes that have come up after your search\nRemember you can press "G" to cleanup those genes that you have selected\n\n')
# if(length(SETTINGS[[ 'genes' ]]) < 50){
# cat(SETTINGS[[ 'genes' ]], sep=" ")
# }else{
# cat(sample(SETTINGS[[ 'genes' ]])[1:50], sep=" ")
# cat('\n')
# }
# SETTINGS[[ 'renameFlag' ]] <- TRUE
# SETTINGS[[ 'heatMapFlag' ]] <- TRUE
# SETTINGS[[ 'biPlotFlag' ]] <- TRUE



# geneDF <- tsSuper$ts_data[SETTINGS[[ 'libraryNames' ]], SETTINGS[[ 'genes' ]][ SETTINGS[[ 'genes' ]] %in% SETTINGS[[ 'geneSubset' ]] ],drop=F]
# geneDF <- matrixWrangler(geneDF, 'log')
# row.names(geneDF) <- SETTINGS[[ 'newLibraryNames' ]]

# load("./profiles/LeeBOI/SETTINGS.Rdata")
# smallClassNames <- names(which(summary(SETTINGS$labelConcat) < 8))
# levels(SETTINGS$labelConcat) <- c(levels(SETTINGS$labelConcat), 'other')
# SETTINGS$labelConcat[SETTINGS$labelConcat %in% smallClassNames] <- 'other'
# newLevels <- unique(SETTINGS$labelConcat)
# SETTINGS$labelConcat <- factor(SETTINGS$labelConcat, levels = newLevels)


# require(MASS)
# r <- lda(geneDF, SETTINGS$labelConcat)
# pr <- predict(r)

# ld1 <- 2
# ld2 <-3

# plot(r$svd^2/sum(r$svd^2), type='b')

# plot(pr$x[,ld1], pr$x[,ld2], pch='')
# text(pr$x[,ld1], pr$x[,ld2], SETTINGS$newLibraryNames, cex=.6)

# par(new=T, xpd=T)
# plot(cr[,ld1], cr[,ld2],  pch='')
# text(cr[,ld2], cr[,ld2], row.names(cr), cex=.3)

# cr <- coef(r)
# plot(sort(abs(cr[,1]) + abs(cr[,2]), TRUE)[1:100])

# heatmap(abs(cor(geneDF)))

