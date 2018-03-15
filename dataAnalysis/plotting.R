#!/usr/bin/Rscript

# alpha diversity measures

png('boxplots.png', width = 1046, height = 411)

par(mfrow = c(1, 3))
par(mar=c(5,6,4,2)+0.1, mgp=c(4,.5,0))

boxplot(allRichness[1:22], allRichness[23:44], names = c('Summer', 'Winter'),
        ylab = 'Species Richness', 
        ylim = c(5000, 13000), cex = 2, cex.lab = 2, cex.axis = 1.5, yaxt="n")
axis(side = 2, at = seq(5000, 13000, 2000), las = 1, cex.axis = 1.5)
minor.tick(ny = 2, nx = 0)

boxplot(allShannon[1:22], allShannon[23:44], names = c('Summer', 'Winter'),
        xlab = 'Season', ylab = 'Shannon Diversity', 
        ylim = c(4, 8), cex = 2,cex.lab = 2, cex.axis = 1.5, yaxt="n")
axis(side = 2, at = seq(4, 8, 1), las = 1, cex.axis = 1.5)
minor.tick(ny = 2, nx = 0)

boxplot(allSimpson[1:22], allSimpson[23:44], names = c('Summer', 'Winter'),
        ylab = 'Simpson Diversity', 
        ylim = c(.75, 1), cex = 2, cex.lab = 2, cex.axis = 1.5, yaxt="n")
axis(side = 2, at = seq(.75, 1, .05), las = 1, cex.axis = 1.5)
minor.tick(ny = 2, nx = 0)

dev.off()

# NMDS plot

treat=c(rep("Summer", 22),rep("Winter", 22))
test = metaMDS(otu.all.norm, k=2)

envData1 = envData
names(envData1)[c(12, 39)] = c('OM', 'VD')

png('nMDS.png', width = 750, height = 500)

ef <- envfit(test, envData1[,c(9, 12, 39)], na.rm = TRUE)
colvec = c('Green', 'Red')

plot(test, display = "sites", type = 'n', 
     xlab = "NMDS axis 1", ylab = "NMDS axis 2")
# draw polygon around sites for each season
ordihull(test, groups = treat, draw="polygon", col= c("white", "grey90"), 
         border = "grey25", label=F)
# add in sites as dots 
points(test, display="sites", pch=16, cex = 1.25,
       col = c(rep("red", 22), rep("blue", 22)))
plot(ef, p.max = .1, col = 'black', cex = 1.75)

dev.off()

# Taxonomy analysis 

insert_minor <- function(major_labs, n_minor) {
  labs <- c(sapply(major_labs, function(x) c(x, '')))
  labs[1:(length(labs)-n_minor)]
}

relative = otu.all.tax
sites = names(relative[ , 2:45])
relative[ , 2:45] = data.frame(t(otu.all.norm))

# convert OTU abundances to relative abundances

relative[, 2:45] = data.frame(apply(relative[, 2:45], 2, function(x) {x/sum(x)}))

# convert to long and get total abundance for each phylum 
# need to detach Hmisc 

detach("package:Hmisc", unload=TRUE)
relative <- gather(relative, site, relativeAbundance, SQ1:WQ22, factor_key=TRUE)

relative_grouped = group_by(relative, phylum, site)
relative = summarize(relative_grouped, totAbundance = sum(relativeAbundance)) %>% data.frame

# filter out phyla that aren't at least 2% of sample

relative = relative[relative$totAbundance > .02, ]

# check for number of phyla present 

phyla = vector()
for(i in sites){
  x = unique(relative[relative$site == i, 1]) %>% length
  phyla = c(phyla, x)
}
phyla

col = c("#ba6800", "#5550d6", "#52c33b", "#c12ab0","#13e0ab", "#f01685",
       "#878000", "#ffa8ed", "#ee5824", "#6e407f", "#e2bc8c", "#a31530")

ggplot(data = relative, aes(x = site, y = totAbundance, fill = phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col, name = 'Phylum') +
  xlab('Site') + ylab("Relative Abundance (Phyla > 2%)") +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  coord_flip() + scale_y_continuous(breaks = seq(0, 1, .125), labels = insert_minor(seq(0,1, .25), 1), 
                                    expand = c(0,0), limits = c(0,1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(color = 'black'))

ggsave('phyla.png', plot = last_plot())

# number of occurances 
relative$
set.seed(1)
relative %>% 
  group_by(phylum) %>%
  summarise(no_rows = length(phylum))

# total abundance for each site 
total = relative %>% 
  group_by(site) %>%
  summarise(abund = sum(totAbundance), std = std(totAbundance)) %>%
  data.frame()
mean(total$abund); std(total$abund)

# avg per season

relative$season = sapply(relative$site, substring, first = 1, last = 1)

total = relative %>% 
  group_by(phylum, season) %>%
  summarise(abund = 100 * mean(totAbundance), std = std(totAbundance)) %>%
  data.frame() 
total[order(total$abund, decreasing = T),]
total


