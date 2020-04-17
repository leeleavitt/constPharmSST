#load previous single cell data set
load("/Users/kevinchase/Box Sync/SS transcriptome/Data/Transcriptomics/cell.type.blob.061219.Rdata")

#cell.type.blob... currently consists of:
#cell.details.061219 cell.type.counts.061219 def.genes geneGoFinder hkg mouse.gene.desc mouse.gene.go tpm TypeSig
#cell.type.counts.061219 row.names are gene ids.  rows are genes columns are samples with descriptors. raw counts.  first column has Gene.name
#cell.details.061219 row.names are the gnomex ids.
#mouse.gene.go has all go terms for gene ids.
#hkg is a list of house keeping gene names
#geneGoFinder is function that takes a vector of strings and returns the names of all genes with those strings in their GO terms e.g.
#geneGoFinder(c("ion channel"))
#TypeSig deprecated
#tpm gene count values scaled 1 million for each sample.
#def.genes is a list of genes used to define type (e.g. Calca)
#mouse.gene.desc is a dataframe with gene descriptions.

#Create an updated blob with data from 17848R
new.dets <- read.csv("~/Downloads/bioinformatics-analysis-A6186/Sample.Details.17848R.csv",header=T,as.is=T)

#identify missing columns from the new data and fill in the critical ones 
#set as missing (NA) the non-critical ones
miss.names <- setdiff(names(cell.details.061219),names(new.dets))
miss.names
names(new.dets)[3] <- "Sample.Name"
new.dets[,"label"] <- paste(make.names(new.dets[,"type"]),new.dets[,"gnomex"],sep=".")
new.dets[,"label"]
miss.names <- setdiff(names(cell.details.061219),names(new.dets))

for(i in miss.names){new.dets[,i] <-NA}

row.names(new.dets) <- new.dets[,"gnomex"]

cell.details.032020 <- rbind(cell.details.061219,new.dets[,names(cell.details.061219)])

#now add the counts to the counts dataframe
new.counts <- read.table("~/Downloads/bioinformatics-analysis-A6186/all_counts.txt",header=T,row.names=1)
new.counts[1:5,1:5]
dim(cell.type.counts.061219)
dim(new.counts)
gi <- !is.element(row.names(new.counts),row.names(cell.type.counts.061219))
ng <- names(tail(sort(apply(new.counts[gi,],1,sum)),n=100))
#some expression of genes previously not detected
mouse.gene.desc[is.element(mouse.gene.desc[,"Gene.stable.ID"],ng),]
#these are not annotated and likely not really genes omit for now but worry about it a 2:00 AM when you 
#wake up and stress about stupid shit.

gi <- setdiff(row.names(cell.type.counts.061219),row.names(new.counts))
sort(apply(cell.type.counts.061219[gi,-1],1,sum))

#some more "predicted gene" squences to stress about
#trim to common gene set
gi <- intersect(row.names(cell.type.counts.061219),row.names(new.counts))
length(gi)
names(new.counts) <- cell.details.032020[sub("^X","",names(new.counts)),"label"]
cell.type.counts.032020 <- cbind(cell.type.counts.061219[gi,],new.counts[gi,])
cbind(names(cell.type.counts.032020)[-1],cell.details.032020[,"label"])
#not the same order fix that shit

setequal(names(cell.type.counts.032020)[-1],cell.details.032020[,"label"])
norder <- sub(".*\\.","",names(cell.type.counts.032020)[-1])
cell.details.032020 <- cell.details.032020[norder,]
sum(names(cell.type.counts.032020)[-1] != cell.details.032020[,"label"])

#OTAY
#now generate the normalized TPM set with gene names and scaled to 1 million total counts for each sample.
#I know this is a bad normalization but the others make me more nervous.
tpm.032020 <- cell.type.counts.032020
row.names(tpm.032020) <- make.names(tpm.032020[,1],unique=T)
tpm.032020 <- tpm.032020[,-1]
tpm.032020[1:5,1:5]
sum(sapply(tpm.032020,is.numeric))
cnt <- apply(tpm.032020,2,sum)
sort(cnt)
tpm.032020 <- sweep(tpm.032020,2,apply(tpm.032020,2,sum)/1000000,'/')
apply(tpm.032020,2,sum)

#finally save the new blob with the new name
save(list=c("cell.type.counts.032020","cell.details.032020","mouse.gene.go","hkg","geneGoFinder","TypeSig","tpm.032020","def.genes","mouse.gene.desc"),file="cell.type.blob.032020.Rdata")


#adjust the naming and the labels grrrrr
names(cell.type.counts.032020) <- sub("^t1","L1",names(cell.type.counts.032020))
names(cell.type.counts.032020) <- sub("^t2","L2",names(cell.type.counts.032020))
names(cell.type.counts.032020) <- sub("^t3","L3",names(cell.type.counts.032020))
names(cell.type.counts.032020) <- sub("^t4","L4",names(cell.type.counts.032020))
names(cell.type.counts.032020) <- sub("^t5","L5",names(cell.type.counts.032020))
names(cell.type.counts.032020) <- sub("^t6","L4",names(cell.type.counts.032020))

