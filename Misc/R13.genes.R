#load the latest transcriptome data
load("/Users/kevinchase/Box Sync/SS transcriptome/Data/Transcriptomics/cell.type.blob.032020.Rdata")
load("./Misc/cell.type.blob.032020.Rdata")
dets <- cell.details.032020
tpm <- tpm.032020
unique(dets$type)
#omit samples with low reads
gi <- dets[,"M.Assigned"] > .5 
#set type variable for approved reads
t.type <- dets[gi,"type"]
# sett all types with n < 4 to "other"
t.cnt <- tapply(rep(1,length(t.type)),t.type,sum)
t.type[is.element(t.type,names(t.cnt[t.cnt < 4]))] <- "other"
t.type <- factor(t.type,levels=names(sort(table(t.type))))
summary(t.type)

#use p.type to targe R13 vs. all others
p.type <- t.type=="R13"

#set highlight for ipsi samples.
t.ipsi <- grepl("ipsi",dets[gi,"label_experiment"])
table(t.ipsi)
table(p.type)


#matrix specific to gene type.
#calcium and sodium channel
#quick check to see the major GO terms for Calcium channels
g0 <- grep("Cacn",row.names(tpm),value=T)
tail(sort(table(as.character(mouse.gene.go[is.element(mouse.gene.go[,"Gene.name"],g0),"GO.term.name"]))))

# manually set go terms 
g0.terms <- c('integral component of membrane','ion channel activity')
g1.terms <- c('integral component of membrane','G-protein coupled receptor activity')
g1 <- geneGoFinder(g1.terms)
g0 <- geneGoFinder(g0.terms)
g4 <- g1 # just use ion channels for now, union(g1,g0) gives ion channel and GPCR
g4 <- intersect(g4,row.names(tpm)) 
gt0func <- function(x){sum(x > 0)}
g4.cnt <- apply(tpm[g4,],1,gt0func) 
g4 <- g4[g4.cnt > 4] #only genes with > 4 cells showing some expression.
g4
k.mat <- t(tpm[g4,]) #matrix of predictors.

#heatmap
hm.mat <- log(k.mat+1) 
hm.breaks <- c(0,seq(1,ceiling(max(as.vector(as.matrix(hm.mat)))),length.out=9))
library(RColorBrewer)
hm.col <- c("white",brewer.pal(4,"Blues"),brewer.pal(4,"Reds"))
#put this all together and show raw expression levels instead of log scale data.
#output to pdef

pdf("Scn.heatmap.pdf",height=12,width=16)
n <- length(hm.breaks)-1
par(fig=c(.2,1,.2,1))
heatmap(as.matrix(hm.mat),scale="none",breaks=hm.breaks,col=hm.col,margins=c(12,2),cexCol=1.5,cexRow=.5)
par(fig=c(.4,.8,0,.1), new=TRUE,mar=c(0,0,0,0))
symbols(1:n,rep(0,n),rectangles=matrix(c(rep(1,n),rep(.3,n)),ncol=2),bg=hm.col,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",inches=F)#,ylim=c(-.01,.01)
text(1:n,rep(0,n),round(exp(hm.breaks[-length(hm.breaks)]),0)-1)#,adj=5)
graphics.off()
system("open Scn.heatmap.pdf")

#randomForest prediction of type
#note use p.type for targeter R13
#or use t.type for the general types
install.packages('randomForest')
library(randomForest)
rft <- randomForest(log(k.mat[gi,]+1),as.factor(p.type),ntree=5000)
imp <- rft$importance[order(rft$importance[,1]),]
#permutation tests on this set of genes
num.perm <- 100
#matrix to hold the permuted results
perm.mat <- matrix(nrow=length(imp),ncol=num.perm)
for(i in 1:100)
{
	rft.perm <- randomForest(log(k.mat[gi,]+1),as.factor(sample(p.type)),ntree=5000) #note sample() creates a random order
	perm.mat[,i] <- rft.perm$importance[,1]
}

tcol <- rgb(.5,.5,.5,.1) # transparent grey color
#output graph to png
png("GPCR.R13.rft.png",height=7*200,width=11*200,pointsize=24)
plot(imp,xlab="Rank",ylab="Importance for predicting cell type",main="Ion channel genes predict cell type R13 (randomForest)",cex=1.5,lwd=4)
oplot <- function(x){points(jitter(rank(x)),x,col=tcol,pch=16)}
apply(perm.mat,2,oplot)
legend("topleft",c("Actual","100 Random"),pch=c(1,16),col=c("black",tcol))
imp.t <- .26 #manual threshold set based on visual inspection
imp.i <- sum(imp > imp.t)
text(seq(1,length(imp))[imp > imp.t],imp[imp > imp.t],pos=2,names(imp[imp>imp.t]),cex=1,srt=0)
dev.off()

# genes of interest
bx.names <- names(imp[imp>imp.t])

#correlation heatmap
ct <- cor(data.frame(log(k.mat[gi,bx.names]+1))) #correlation matrix 
library(gplots)
graphics.off()
heatmap.2(ct,trace="none",cellnote=as.matrix(round(ct,2)),notecol="black",keysize=.75,key.title="Correlation",cexRow=1,cexCol=1)

#boxplots of genes of interes.
starfunc <- function(tp)
{
	sigc <- ""
	if(tp < .1){sigc <- "-"}
	if(tp < .05){sigc <- "*"}	
	if(tp < .01){sigc <- "**"}
	if(tp < .001){sigc <- "***"}
	return(sigc)
}

for(i in bx.names)
{
	#main.lab <- paste(i,sub("\\[.*\\]","",mouse.gene.desc[match(i,mouse.gene.desc[,"Gene.name"]),"description"]))
	main.lab <- 'hi'
	png(paste(i,"r13.png",sep=""),height=5*200,width=11*200,pointsize=48)
	
	#regression to estimate significance
	#note relevel() puts "other" as the fist factor so it is rolled into the mean (intercept)
	#this is an unavoidable aspect of regression analysis when the intercept is in the model
	glt <- glm(log(k.mat[gi,i]+1) ~ relevel(t.type,"other")) 
	# summary(glt) gives a summary of the regression fit and statistics
	strs <- lapply(coefficients(summary(glt))[,4][-1],starfunc) #
	bx1 <- boxplot(log(k.mat[gi,i]+1) ~ as.factor(t.type),main=main.lab,xlab="",ylab="log(TPM)",range=0, varwidth=T)
	x.jit <- jitter(as.integer(as.factor(t.type)))
	points(x.jit,log(k.mat[gi,i]+1),pch=c(1,16)[as.integer(t.ipsi)+1],cex=.5)
	text(1:(length(strs)),bx1$stats[3,1:(length(strs))],strs,pos=1,cex=1.25)
	dev.off()
}

paste0(as.character(levels(t.type)), ":", as.character(summary(t.type)))
