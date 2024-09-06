
######################################################
#Try HMSC? Similar model based approach to ordination
library(Hmsc)

fit <- Hmsc(ordispe, XData = ordi_extras, XFormula = ~NH4M + NO3M + SO4SM, distr="normal")
fit2 = sampleMcmc(fit, thin = 2, samples = 100, transient = 50, nChains = 2, verbose = 50, nParallel = 2)
mpost = convertToCodaObject(fit2)
summary(mpost$Beta)
m <- fit2
preds = computePredictedValues(m) 
evaluateModelFit(hM=m, predY=preds)

postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95)

rL$nfMin=2
rL$nfMax=2
m = Hmsc(Y=Y, XData=XData, XFormula=~1, studyDesign=studyDesign, ranLevels=list(sample=rL))
etaPost=getPostEstimate(m, "Eta") 
lambdaPost=getPostEstimate(m, "Lambda")
biPlot(m, etaPost = etaPost, lambdaPost = lambdaPost, factors = c(1,2), "x2")

# Before we look at our latent variables, let's have a look at collinearity. The following
# code shows any significant collinearity between our environmental covariates.

source("http://www.sthda.com/upload/rquery_cormat.r")
colin <- rquery.cormat(ordi_extras, type = "flatten", graph = FALSE)
colin$r %>% filter(abs(cor) > 0.5)

# row    column   cor        p
# 1 DDEG0      MIND -0.89 6.9e-141
# 2 DDEG0 ELEVATION -0.99  0.0e+00
# 3  MIND ELEVATION  0.90 4.0e-148

# The high cor values and low p values indicate some serious collinearity between positive degree days, 
# elevation and moisture index. 

# Let's see if our latent variables correspond to any of these covariates.

# We define colours according to the values of covariates. The darker blue indicates a higher value
# of the relevant covariate.

par(mfrow=c(2,2))
for (i in 1:length(colnames(ordi_extras))) {
  covariate <- ordi_extras[,i]
  rbPal <- colorRampPalette(c('mediumspringgreen', 'blue'))
  Colorsph <- rbPal(20)[as.numeric(cut(covariate, breaks = 20))]
  breaks <- seq(min(covariate), max(covariate), length.out = 30)
  
  ordiplot(fit_base, main = paste0("Ordination of sites, color: ",colnames(ordi_extras)[i]),
           symbols = TRUE, s.colors = Colorsph, xlim = c(-1.2,1.2), ylim = (c(-1.2, 1.2)))
}

# We can see some quite clear gradients related to the four collinear variables we mentioned
# above. At this point let's take one of the two climate related covariates that could have a
# direct impact on vegetation (DDEG0) and SLOPE, since MIND is so collinear to DDEG0, and
# slope might have a more direct impact than elevation on our community.

# Let's use these two and use some code from last session to figure out how many latent
# variables would be appropriate. Basically we're running the same model again and again, but 
# increasing the number of latent variables each time. I'm using the lowest AICc value to determine the
# model which fits best.

fit_list <- list()
for(i in 0:3){
  fit_sub <- gllvm(ordispe, ordi_extras, family = binomial(link="probit"), Ntrials = 10, num.lv = i, sd.errors = FALSE, 
                   formula = ~ SLOPE + MIND, seed = 1234)
  fit_list[[i+1]] <- fit_sub
}

# Let's have a look at how our AICc values look.

AICc <- sapply(fit_list, function(X) {summary(X)$AICc})
data.frame(AICc, model= paste0("LV-",0:3))

# AICc      model
# 14687.36  LV-0
# 12957.51  LV-1
# 12625.36  LV-2
# 12970.21  LV-3

# We can see that the best model here uses 2 latent variables. Let's have a look at 
# how these variables compares to our remaining environmental covariates. 

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
remaining_covariates <- c("DDEG0","SOLRAD","ELEVATION","TPI")

for(i in 1:length(remaining_covariates)) {
  covariate <- X[,remaining_covariates[i]]
  rbPal <- colorRampPalette(c('mediumspringgreen', 'blue'))
  Colorsph <- rbPal(20)[as.numeric(cut(covariate, breaks = 20))]
  breaks <- seq(min(covariate), max(covariate), length.out = 30)
  ordiplot(fit_list[[3]], main = paste0("Ordination of sites, color: ",remaining_covariates[i]),
           symbols = TRUE, s.colors = Colorsph, xlim = c(-1.2,1.2), ylim = (c(-1.2, 1.2)))
}

# From this we can see that there is still a bit of variation explained by positive degree days,
# despite its collinearity with moisture index (and by elevation, but we'll focus on that tomorrow).
# Let's see what happens when we include degree days over zero.

fit_DegreeDays <- gllvm(Y, X, family = binomial(link="probit"), num.lv = 2, sd.errors = FALSE, 
                        formula = ~ NH4M + NO3M + SO4SM, seed = 1234)

summary(fit_DegreeDays)$AICc

# [1] 12743.35

# We can see the AICc values stay pretty much the same, even rising a bit. But leaving it out means 
# we may attribute variation to our latent variable that is the result of the environment.


### EXTENDED QUESTION ####
# Have a look at the coefficient effects using the basic command below, after switching sd.errors to 
# TRUE in your gllvm commands. How do the covariate effects change with the introduction of new variables?

# coefplot(fit_Elevation, cex.ylab = 0.5)

###################################################
# Breakout questions #
# What is a species association?
# If we accounted for all possible covariates producing environmental variation, 
# would we still see species associations?
###################################################

# What I've previously done here is group species together based on approximately what elevation
# their occurrence peaks at. This left us with three groups; montane, subalpine and alpine species.
# The colour plots below mean we can see each group easily in our ordination plots.

# NB: The 'ordiplot' function below is a product of a recent update. If you find that it's not 
# working, I've attached a fully coded function at the end of this script as a substitute.

colour.groups <- c("red","blue","green")[ordi_extras$survey_year]

par(mfrow=c(1,1))
ordiplot(fit_base, biplot=TRUE, main = "Ordination of sites: no covariates",
         symbols = TRUE, s.colors = "white", xlim = c(-4,4),ylim=c(-3,3))

# There are some very obvious trends here. Let's see what happens when we introduce MIND and DDEG0.

ordiplot(fit_list[[3]], biplot=TRUE, main = "Ordination of sites: two covariates",
         symbols = TRUE, s.colors = "white", xlim = c(-4,4),ylim=c(-3,3))

# And now when we introduce degree days as an extra covariate.

ordiplot(fit_DegreeDays, biplot=TRUE, main = "Ordination of species: three covariates",
         symbols = TRUE, s.colors = "white", xlim = c(-4,4),ylim=c(-3,3), spp.colors=colour.groups)


# You can see that the species group together more clearly, as the effect of the latent variable becomes
# weaker.

