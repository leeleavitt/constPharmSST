#I've packaged several dataframes and functions in a "blob" file for convenience
load("cell.type.blob.032020.Rdata")
#cell.type.counts.032020 row.names are gene ids.  rows are genes columns are samples with descriptors. raw counts.  first column has Gene.name
#cell.details.032020 row.names are the gnomex ids.
#mouse.gene.go has all go terms for gene ids.
#hkg is a list of house keeping gene names
#geneGoFinder is function that takes a vector of strings and returns the names of all genes with those strings in their GO terms e.g.
#geneGoFinder(c("ion channel"))
#TypeSig deprecated
#tpm gene count values scaled 1 million for each sample.
#def.genes is a list of genes used to define type (e.g. Calca)
#mouse.gene.desc is a dataframe with gene descriptions.

#using tpm to run statistical tests
det <- cell.details.032020 #easier var to work with
tpm <- tpm.032020
tpm.res <- data.frame(mean=apply(tpm,1,mean))
tpm.res[,"sd"] <- apply(tpm,1,sd)
plot(tpm.res,log="xy")
#some very high means
high.genes <- row.names(tail(tpm[order(tpm.res[,1]),]))
high.genes
mouse.gene.desc[match(high.genes,mouse.gene.desc[,"Gene.name"]),]
#some mitochondrial genes
#WTF these cells are pumping out the Tac1
#note from genecards "Substance P is an antimicrobial peptide with antibacterial and antifungal properties"
#is this a measure of bacterial contamination?
plot(t(tpm[c("Tac1","Fstl1"),])+1,log="xy")
tpm.res[,"max"] <- apply(tpm,1,max)
library(randomForest)
#limit the number of genes tested based on expression levels.
sum(tpm.res[,"max"] > 20)
gi <- tpm.res[,"max"] > 20
#24072 genes 
#this still might make randomForest choke and die
rft <- randomForest(t(tpm[gi,]),as.factor(det[,"type"]))
b.names <- names(tail(rft$importance[order(rft$importance[,1]),]))
mouse.gene.desc[match(b.names,mouse.gene.desc[,"Gene.name"]),]
#not bad I'll increase the ntree setting to get a better estimate of significance.
rft <- randomForest(t(tpm[gi,]),as.factor(det[,"type"]),ntree=5000)
b.names <- names(tail(rft$importance[order(rft$importance[,1]),]))
mouse.gene.desc[match(b.names,mouse.gene.desc[,"Gene.name"]),]
tpm.res[gi,"rft.all"] <- rft$importance[,1]
#OK just limit tpm.res to the gi 
tpm.res <- tpm.res[gi,]
#so that's a measure of each genes importance in predicting type (i.e. all types)
#Let's test for genes important in individual types vs. all others.
#actually some of the samples have very low aligned counts and should be excluded from analysis.
sort(det[,"M.Aligned"])
ji <- det[,"M.Aligned"] > 5
sort(table(det[ji,"type"]))
#test for genes that inform t1 type. (note the use of gi and ji to select genes with high counts and cells with high alignment)
rft <- randomForest(t(tpm[gi,ji]),as.factor(det[ji,"type"]=="t1"),ntree=5000)
tpm.res[,"rft.t1"] <- rft$importance[,1]

test.names <- names(tail(sort(table(det[ji,"type"])),12))
#run through all the pertinent factors
for(i in test.names[3:length(test.names)])
{
	rft <- randomForest(t(tpm[gi,ji]),as.factor(det[ji,"type"]==i),ntree=5000)
	tpm.res[,paste("rft",i,sep=".")] <- rft$importance[,1]	
}

#combine the results with the raw counts
dat.out <- cbind(cell.type.counts.032020[gi,],tpm.res)
#mark ion channel genes and GPCR genes
g.ion <- geneGoFinder(c("ion channel"))
dat.out[,"ion.channel"] <- 0
dat.out[is.element(dat.out[,"Gene.name"],g.ion),"ion.channel"] <- 1

gs1 <- geneGoFinder(c("G-protein coupled receptor activity"))
dat.out[,"GPCR"] <- 0
dat.out[is.element(dat.out[,"Gene.name"],gs1),"ion.channel"] <- 1

write.csv(dat.out,file="dat.out.032020.csv")


