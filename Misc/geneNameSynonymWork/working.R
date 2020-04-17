# Following instructions from here, towards the bottom
# https://www.biostars.org/p/278673/

setwd("./leeProcessing/geneNameSynonymWork")

gene_result <- read.delim("gene_result.txt", sep = "\t")

gene_result[gene_result$Symbol == 'Cacna1h',]