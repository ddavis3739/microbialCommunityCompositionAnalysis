#!/usr/bin/Rscript

# read in deblur file 
 
deblur = read.table('reference-hit_all.table.txt')

# add in colnames

colnames(deblur) = c('sOTU', sites)

# make sOTU's rownames and transpose 

row.names(deblur) = deblur$sOTU
deblur = deblur[,-1]

deblur = as.data.frame(t(deblur))

# normalise to account for differences in read numbers across samples (rarefaction)

deblur = rrarefy(deblur, min(rowSums(deblur)))

### get alpha diveristy measures and test

deblurRichness = specnumber(deblur)
deblurShannon = diversity(deblur, index = 'shannon')
deblurSimpson = diversity(deblur, index = 'simpson')

# not sig ; t = 1.9281, df = 38.918, p-value = 0.06116)
t.test(deblurRichness[1:22], deblurRichness[23:44])

# sig; W = 339, p-value = 0.02243
wilcox.test(deblurShannon[1:22], deblurShannon[23:44])

# sig; W = 331, p-value = 0.03686
wilcox.test(deblurSimpson[1:22], deblurSimpson[23:44])


# between tests

par(mfrow = c(1, 3))
boxplot(allRichness, deblurRichness, names = c('OTU', 'Deblur'),
        ylab = 'Species Richness', 
        ylim = c(4000, 14000), cex.lab = 2, cex.lab = 1.5, yaxt="n")
axis(side = 2, at = seq(4000, 14000, 2000), las = 1)
minor.tick(ny = 2, nx = 0)

# summer higher for shannon and simpson 

boxplot(allShannon, deblurShannon, names = c('OTU', 'Deblur'),
        xlab = 'Season', ylab = 'Shannon Diversity', 
        ylim = c(3, 9), cex.lab = 2, cex.lab = 1.5, yaxt="n")
axis(side = 2, at = seq(3, 9, 1), las = 1)
minor.tick(ny = 2, nx = 0)

boxplot(allSimpson, deblurSimpson, names = c('OTU', 'Deblur'),
        ylab = 'Simpson Diversity', 
        ylim = c(.7, 1), cex.lab = 2, cex.lab = 1.5, yaxt="n")
axis(side = 2, at = seq(.7, 1, .05), las = 1)
minor.tick(ny = 2, nx = 0)

# not sig ; t = 0.8274 df = 81.22, p-value > .05)
t.test(allRichness, deblurRichness)

# sig; W = 606, p-value < .01
wilcox.test(allShannon, deblurShannon)

# sig; W = 709, p-value < 0.05
wilcox.test(allSimpson, deblurSimpson)

### get beta measures and test 

deblurCommSim <- vegdist(deblur, "bray")

permDeblur = adonis(deblurCommSim ~ Season * 
                      (pH + OM_0.10 + plantDiv +vegDen +   RB_0.30 +
                         Elevation_.m.), 
                 data = envData, permutations = 10000)
permDeblur

# test for equal dispersion  
# is homogenous; (F1, 42 = 1.71, p > 0.05)

betad = betadisper(deblurCommSim, envData$Season)
anova(betad)
plot(betad)

# create NMDS plot to visualize

treat=c(rep("Summer", 22),rep("Winter", 22))

test = metaMDS(deblur, k=2)

envData1 = envData
names(envData1)[c(12, 39)] = c('OM', 'VD')

png('nMDS_deblur.png', width = 750, height = 500)
ef <- envfit(test, envData1[,c(9, 12)], na.rm = TRUE)
colvec = c('Green', 'Red')

plot(test, display = "sites", type = 'n', 
     xlab = "NMDS axis 1", ylab = "NMDS axis 2")
# draw polygon around sites for each season
ordihull(test, groups = treat, draw="polygon", col= c("white", "grey90"), 
         border = "grey25", label = F)
# add in sites as dots 
points(test, display = "sites", pch = 16, cex = 1.25,
       col=c(rep("red", 22),rep("blue", 22)))
plot(ef, p.max = .1, col = 'black', cex = 1.75)
dev.off()

############################
# Alpha effects on Plant div
############################
# no effect of alpha diverisy on plant diversity 

x = glm(vegDen ~ deblurRichness, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

x = glm(plantDiv ~ deblurRichness, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

x = glm(vegDen ~ deblurShannon, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')
cor(envData$vegDen, allShannon)

x = glm(plantDiv ~ deblurShannon, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

x = glm(vegDen ~ deblurSimpson, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')
cor(envData$vegDen, allSimpson)

x = glm(plantDiv ~ deblurSimpson, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

### testing for alpha effects on nutrient flux

nuts = envData[c(1:12, 23:26, 28:34), c(1, 14:29)]
nutRich = deblurRichness[c(1:12, 23:26, 28:34)]
nutShannon = deblurShannon[c(1:12, 23:26, 28:34)]
nutSimpson = deblurSimpson[c(1:12, 23:26, 28:34)]

grabCoef = function(x){
  a = summary(x)
  a = a$coefficients
  a = a[2,4]
  b = anova(x, test = 'F')
  return(b$`Pr(>F)`[2])
}

alpha <- cbind(nutRich, nutShannon, nutSimpson)
models <- list();models1 = list(); models2 = list()

for (i in 2:17){
  models[[names(nuts)[i]]] = 
    glm(nuts[,i] ~ nutRich, data = nuts, family = gaussian()) 
}
for (i in 2:17){
  models1[[names(nuts)[i]]] = 
    glm(nuts[,i] ~ nutShannon, data = nuts, family = gaussian()) 
}
for (i in 2:17){
  models2[[names(nuts)[i]]] = 
    glm(nuts[,i] ~ nutSimpson, data = nuts, family = gaussian()) 
}

rich = unlist(lapply(models, grabCoef))
shan = unlist(lapply(models1, grabCoef))
simp = unlist(lapply(models2, grabCoef))

richSig = vector(); shanSig = vector(); simpSig = vector()
for (i in rich){
  if(i < .001) a = ('***')
  else if(i < .01) a = ('**')
  else if(i < .05) a = ('*')
  else a = ('')
  richSig = c(richSig, a)
}
for (i in shan){
  if(i < .001) a = ('***')
  else if(i < .01) a = ('**')
  else if(i < .05) a = ('*')
  else a = ('')
  shanSig = c(shanSig, a)
}
for (i in simp){
  if(i < .001) a = ('***')
  else if(i < .01) a = ('**')
  else if(i < .05) a = ('*')
  else a = ('')
  simpSig = c(simpSig, a)
}
alphaCors = t(cor(alpha, nuts[,-1]))

nutEffects = data.frame(alphaCors[,1], rich, richSig, alphaCors[,2], shan, shanSig, 
                        alphaCors[,3], simp, simpSig)
colnames(nutEffects)[c(1, 4, 7)] = c('richCor', 'shanCor', 'simpCor')
nutEffects
plot(nutRich, nuts$LightSilicateFlux)


# DarkNitriteFlux; (F1, 21 = 7.7191, p < 0.01) 
x = glm(nuts[,13] ~ nutRich, data = nuts, family = gaussian())
anova(x, test = 'F')

# LightNPOCFlux; shan (F1, 21 = 9.751, p < 0.01) 
x = glm(nuts[,4] ~ nutShannon, data = nuts, family = gaussian())
anova(x, test = 'F')

# LightNitriteFlux; shan (F1, 21 = 5.9338, p < 0.05)

x = glm(nuts[,12] ~ nutShannon, data = nuts, family = gaussian())
anova(x, test = 'F')



