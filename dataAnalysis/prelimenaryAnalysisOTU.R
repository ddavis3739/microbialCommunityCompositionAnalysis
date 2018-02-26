#!/usr/bin/Rscript

# load packages 

require(ggplot2); library(vegan); require(edgeR); library(car); require(Hmisc)
require(dplyr); require(tidyr); library(ggsci); 

std <- function(x) sd(x)/sqrt(length(x))

# read in OTU tables and environmental data

otu.all.tax = read.csv('OTUall.csv', header = T)

envData = read.csv('env_Variables.csv', header = T)

# change size column to numeric 

otu.all.tax$size = sapply(otu.all.tax[['size']], 
                          function(x) as.numeric(strsplit(as.character(x), split = '=')[[1]][2]))

# rename taxa levels

labs = c('domain', 'phylum', 'class', 'order', 'family', 'genus')

otu.all.tax = otu.all.tax[, -c(48, 51, 54, 57, 60, 63)]
names(otu.all.tax)[c(47, 49, 51, 53, 55, 57)] = labs

### diversity analysis

sites = names(otu.all.tax[ , 2:45])

# save tax in seperate dataframe

allTax = otu.all.tax[, (c(1, (ncol(otu.all.tax) - 12):ncol(otu.all.tax)))]

# remove taxonomic data

otu.all = otu.all.tax[, -((ncol(otu.all.tax) - 12):ncol(otu.all.tax))]

# make rownames OTUS

row.names(otu.all) = otu.all$OTU_ID
otu.all = otu.all[,-1]

# transpose 

otu.all = as.data.frame(t(otu.all))
colnames(otu.all) = otu.all.tax$OTU_ID

## Diversity Measures
# normalise to account for differences in read numbers across samples (rarefaction)

otu.all.norm = rrarefy(otu.all, min(rowSums(otu.all)))

# Alpha Diversity 
# calc OTU richness

allRichness = specnumber(otu.all.norm)

shapiro.test(allRichness[1:22]); shapiro.test(allRichness[23:44]); 
names(allRichness) = sites
t.test(allRichness[1:22], allRichness[23:44])
# t = 2.754, df = 41.701, p-value = 0.008684, sig difference in OTU richness

# alpha evenness measures 

allShannon = diversity(otu.all.norm, index = 'shannon')
allSimpson = diversity(otu.all.norm, index = 'simpson')

shapiro.test(allShannon[1:22]); shapiro.test(allShannon[23:44])
# w = 431, p-value < .001; summer higher
wilcox.test(allShannon[1:22], allShannon[23:44])

shapiro.test(allSimpson[1:22]); shapiro.test(allSimpson[23:44])
# w = 434, p-value < .001; summer higher 
wilcox.test(allSimpson[1:22], allSimpson[23:44])

## Beta diversity 
# produces matrix with similarities between each site 

allCommSim <- vegdist(otu.all.norm, "bray")








