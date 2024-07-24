#!/usr/bin/env Rscript
library(argparse)
library(clusterProfiler)
parser <- argparse::ArgumentParser(description = "Rstrip for Go enrichment")
parser$add_argument("-g","--gene",action = "store",type = "character", 
		    help = "one colums only,each line is the gene that u want to attend to ")
parser$add_argument("-goa","--GOannotation",action = "store",type = "character",
		    help = "one colums only,each line is the gene that u want to attend to ")
parser$add_argument("-Gotb","--Gotb",action = "store",type = "character")
parser$add_argument("-out","--output_path",action = "store",type = "character")
parser$add_argument("-p","--p_value_cutoff",action = "store",type = "character")
args <- parser$parse_args()
GOannotation<- args$GOannotation
Gotb<- args$Gotb
gene<- args$gene
output_path<- args$output_path
p_value_cutoff<-args$p_value_cutoff
GOannotation <- read.delim(GOannotation,stringsAsFactors=FALSE)
Gotb <-read.delim(Gotb,stringsAsFactors=FALSE)
gene <- read.table(gene, header = FALSE)
gene <- as.factor(gene$V1)
go_enrich <- enricher(gene,pAdjustMethod="BH",pvalueCutoff = p_value_cutoff,TERM2GENE=GOannotation[c(2,1)],TERM2NAME=Gotb[1:2])
write.table(go_enrich,file=output_path,sep="\t",quote = FALSE,row.names = FALSE)
print(go_enrich)
