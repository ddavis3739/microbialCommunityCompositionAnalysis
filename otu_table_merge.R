#!/usr/bin/Rscript

if(!require("Biostrings")){
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
}
require(Biostrings)
library(dplyr)
library(tidyr)

args <- commandArgs()
print(args)

##### For bacteria -> paired based on rdp 16S ##### 
otu_table = read.table(args[6], header = F, skip = 2)

otuName = strsplit(readLines(args[6],n = 2)[2], split = '')
otuName[[1]][5] = '_'
otuName = otuName[[1]][-1]
otuName= paste(otuName, collapse = '')
otuName = strsplit(otuName, '\t')
otuName = otuName[[1]]

names(otu_table) = otuName

taxonomy <- read.delim(args[7],header = F)
names(taxonomy) <- c("OTU_ID",
                       'blank',
                       "Root",
                       "rootrank",
                       "rootconfidence",
                       paste(rep(c("taxa",
                                   "taxa_level",
                                   "confidence"),(ncol(taxonomy)-4)/3),rep(1:((ncol(taxonomy)-4)/3), each=3),sep="_"))

taxonomy = taxonomy[,-2]

taxonomy <- separate(taxonomy, OTU_ID, c("OTU_ID", "size"), sep=";")

complete <- merge(otu_table, taxonomy, by="OTU_ID")

filename = paste0('../OTUtable_merged_', args[8], '.csv')
write.csv(complete, file= filename, row.names = FALSE)
