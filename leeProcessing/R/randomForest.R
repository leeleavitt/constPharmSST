# library(randomForest)
# rft <- randomForest::randomForest(geneDF, labelConcat)

# imp <- rft$importance[order(rft$importance[,1], decreasing = TRUE),]

# imp[1:100]
# hist( imp[imp>0] , length(imp)/2)

# ?order