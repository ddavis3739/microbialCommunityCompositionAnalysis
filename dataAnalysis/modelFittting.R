#!/usr/bin/Rscript

# test environmental factor for collinearity

envData$richness = allRichness
envData$shannon = as.numeric(allShannon)
envData$simpson = allSimpson

colinTest = envData[c(8:13,39, 40)]

# check correlations 
library(GGally); library(car)
ggpairs(colinTest)

library(corpcor)
mat = cov(colinTest) %>% cor2pcor 
colnames(mat) = names(colinTest)
rownames(mat) = names(colinTest)
mat

# Farrar-Glauber test for each alpha measure
library(mctest)
omcdiag(colinTest, envData$richness)
imcdiag(colinTest, envData$richness)

omcdiag(colinTest, envData$shannon)
imcdiag(colinTest, envData$shannon)

omcdiag(colinTest, envData$simpson)
imcdiag(colinTest, envData$simpson)

# compute compute the partial correlation 
# BD and SMC semi to be colinear, remove SMC
library(ppcor)
pcor(colinTest, method = "pearson")

# retest without SMC
colinTest1 = envData[c(8:10, 12:13, 39, 40)]

ggpairs(colinTest1)

mat = cov(colinTest1) %>% cor2pcor 
colnames(mat) = names(colinTest1)
rownames(mat) = names(colinTest1)
mat

# Farrar-Glauber test for each alpha measure
omcdiag(colinTest1, envData$richness)
imcdiag(colinTest1, envData$richness)

omcdiag(colinTest1, envData$shannon)
imcdiag(colinTest1, envData$shannon)

omcdiag(colinTest1, envData$simpson)
imcdiag(colinTest1, envData$simpson)

# compute compute the partial correlation 
# BD and OM semi to be colinear, remove BD
pcor(colinTest1, method = "pearson")

# retest without BD
colinTest2 = envData[c(8:9, 12:13, 39, 40)]
names(envData[c(8:9, 12:13, 39, 40)])
ggpairs(colinTest2)

mat = cov(colinTest2) %>% cor2pcor 
colnames(mat) = names(colinTest2)
rownames(mat) = names(colinTest2)
mat

# Farrar-Glauber test for each alpha measure
omcdiag(colinTest2, envData$richness)
imcdiag(colinTest2, envData$richness)

omcdiag(colinTest2, envData$shannon)
imcdiag(colinTest2, envData$shannon)

omcdiag(colinTest2, envData$simpson)
imcdiag(colinTest2, envData$simpson)

# compute compute the partial correlation 
# BD and OM semi to be colinear, remove BD
pcor(colinTest2, method = "pearson")

### quantify further using lms and related to environmental data ###

# alpha diversity measures 

# test for normality 

allShapTest = vector(mode="numeric", length=0)
for(i in c(6:42)){
  p = shapiro.test(envData[,i])
  p = as.numeric(unlist(p)[2])
  allShapTest = c(allShapTest, p)
}
names(allShapTest) = names(envData[, c(6:42)])
allShapTest[allShapTest < .05]

# richness model 

x = glm(richness ~ Season * 
          (pH + OM_0.10 + plantDiv + vegDen + RB_0.30 +
             Elevation_.m.), 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

# Shannon model 

x = glm(shannon ~ Season * 
          (pH + OM_0.10 + plantDiv + vegDen + RB_0.30 +
             Elevation_.m.), 
        data = envData, family = gaussian())
# test for overdispersion in model to confirm gaussian choice
x$deviance/x$df.residual
anova(x, test = 'F')

# Simpson model 

x = glm(simpson ~ Season * 
          (pH + OM_0.10 + plantDiv + vegDen + RB_0.30 +
             Elevation_.m.), 
        data = envData, family = gaussian())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

##### beta diveristy #########
# PERMANOVA
# prioritize biotic interactions above other measures in stepwise model

permAll = adonis(allCommSim ~ Season * 
                  (pH + OM_0.10 + plantDiv +vegDen +   RB_0.30 +
                  Elevation_.m.),
              data = envData, permutations = 10000)
permAll

# test for equal dispersion between groups 
# OTU not homogenous

betad = betadisper(allCommSim, envData$Season)
anova(betad)
plot(betad)

############################
# Alpha effects on Plant div
############################

# only veg den effected by shannon and simpson measures 
# shan F1, 42 = 7.34, p < 0.01; simp F1, 42 = 12.73, p < 0.001

x = glm(vegDen ~ richness, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

x = glm(plantDiv ~ richness, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

x = glm(vegDen ~ shannon, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')
cor(envData$vegDen, allShannon)

x = glm(plantDiv ~ shannon, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

x = glm(vegDen ~ simpson, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')
cor(envData$vegDen, allSimpson)

x = glm(plantDiv ~ simpson, 
        data = envData, 
        family = quasipoisson())
# test for overdispersion in model to confirm quasipossion choice
x$deviance/x$df.residual
anova(x, test = 'F')

#################################################################
# test to see what nutrient flux measures are effected by alpha diversity measures 
#################################################################

nuts = envData[c(1:12, 23:26, 28:34), c(1, 14:29)]
nutRich = allRichness[c(1:12, 23:26, 28:34)]
nutShannon = allShannon[c(1:12, 23:26, 28:34)]
nutSimpson = allSimpson[c(1:12, 23:26, 28:34)]

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

# LightNPOCFlux; (F1 = 9.11, p < 0.01) 
x = glm(nuts[,4] ~ nutRich, data = nuts, family = gaussian())
anova(x, test = 'F')

# LightSilicateFlux; (F1 = 4.60, p < 0.05) 
x = glm(nuts[,8] ~ nutRich, data = nuts, family = gaussian())
anova(x, test = 'F')

# DarkNitrateFlux; rich (F1 = 4.61, p < 0.05); shan (F1 = 4.96, p < 0.05)
x = glm(nuts[,15] ~ nutRich, data = nuts, family = gaussian())
anova(x, test = 'F')

x = glm(nuts[,15] ~ nutShannon, data = nuts, family = gaussian())
anova(x, test = 'F')

# DarkNOxFlux; (F1 = 4.65, p < 0.05);  shan (F1 = 4.99, p < 0.05)
x = glm(nuts[,17] ~ nutRich, data = nuts, family = gaussian())
anova(x, test = 'F')

x = glm(nuts[,17] ~ nutShannon, data = nuts, family = gaussian())
anova(x, test = 'F')


