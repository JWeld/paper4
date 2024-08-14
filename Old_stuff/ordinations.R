library(vegan)
library(tidyverse)
library(ggfortify)

#Ordinations####

#Arrange data
bryo_env_means_scaled <- filter(env_means_scaled, ID %in% veg_bryo$ID )

common <- intersect(bryo_env_means_scaled$ID, bryo_means$ID)
bryo.ord <- filter(bryo_means, ID %in% common)
bryo.env.ord <- filter(bryo_env_means_scaled, ID %in% common)
combo <- left_join(bryo.ord, bryo.env.ord, by ="ID")
combo$n <- NULL
combo$p <- NULL
combo$N.P <- NULL
combo$survey_year.y <- NULL
combo <- rename(combo, survey_year = survey_year.x)
extra <- dplyr::select(data5, c(ID, mean_age, tree))
extra$mean_age <- as.numeric(as.character(extra$mean_age))
combo2 <- left_join(combo, extra,  by="ID")
combo2 <- combo2 %>% filter(ID %in% dat.all.countries$ID)
#combo2 <- drop_na(combo2)

env.ord <- dplyr::select(combo2, ID, latitude, longitude, n_nh4, n_no3, tree, mean_age,
                         sum_canopy, mean_temp, mean_precip, mean_age,survey_year)
spp.ord <- dplyr::select(combo2, -c(latitude, longitude, n_nh4, n_no3,survey_year, tree, mean_age,
                                    sum_canopy, mean_temp, mean_precip,
                                    survey_year, mean_age, country, code_country))
env.ord <- column_to_rownames(env.ord, var = "ID")
spp.ord <- column_to_rownames(spp.ord, var = "ID")

spe.log <- log1p (spp.ord)  # species data are in percentage scale which is strongly 
#rightskewed, better to transform them
spp.hell <- decostand (spe.log, 'hell') 
#or not?
#spe.hell <- decostand (spp.ord, 'hell')  # we are planning to do tb-RDA, 
#this is Hellinger pre-transformation

#PCA with environmental variables over hellinger spp data
#spe.log <- log1p (spp.ord)  # species data are in percentage scale which is strongly 
#rightskewed, better to transform them
#spe.hell <- decostand (spe.log, 'hell')  # we are planning to do tb-RDA, 
PCA.1 <- rda (spp.hell)  
ef <- envfit (PCA.1 ~ ., data = env.ord, perm = 999, na.rm = TRUE)
ef
ef.adj <- ef 
pvals.adj <- p.adjust (ef$vectors$pvals, method = 'bonferroni')
ef.adj$vectors$pvals <- pvals.adj
ef.adj
summary(PCA.1)
ordiplot (PCA.1, display = 'sites')
plot (ef)

#CCA
bryoplot.cca <- cca(spp.ord ~ ., data = env.ord, na.action = na.omit)
summary(bryoplot.cca)
anova(bryoplot.cca)
plot(bryoplot.cca, display = c("sites","species"), scaling = 3)
text(bryoplot.cca, scaling = 3, display = "bp") 

#hellinger rda with L.x
spe.log <- log1p (spp.ord)  # species data are in percentage scale which is strongly 
#rightskewed, better to transform them
spe.hell <- decostand (spe.log, 'hell')  # we are planning to do tb-RDA, 
#this is Hellinger pre-transformation
tbRDA <- rda (spe.hell ~ ., data = env.ord, na.action = na.omit)  # calculate tb-RDA with all explanatory variables
tbRDA
constrained_eig <- tbRDA$CCA$eig/tbRDA$CA$tot.chi*100
unconstrained_eig <- tbRDA$CA$eig/tbRDA$CA$tot.chi*100
barplot (c(constrained_eig, unconstrained_eig), col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))), las = 2, ylab = '% variation')
ordiplot (tbRDA)
test1 <- anova(tbRDA)
test1
#test2 <- anova(tbRDA, by ='axis')
test2
ordiplot(tbRDA, display = "species")
pl <- ordiplot(tbRDA, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "species", col="blue", cex=0.9)

R2.obs <- RsquareAdj (tbRDA)$r.squared
R2.obs

# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = spe.hell2 ~ n_nh4 + n_no3 + latitude + longitude + sum_canopy + mean_temp + mean_precip + survey_year, data = bryo_env_means2_na, na.action = na.omit)
# Df Variance        F Pr(>F)    
# RDA1       1  0.07158 113.2268  0.001 ***
#   RDA2       1  0.02441  38.6100  0.001 ***
#   RDA3       1  0.00989  15.6508  0.001 ***
#   RDA4       1  0.00709  11.2124  0.001 ***
#   RDA5       1  0.00527   8.3315  0.001 ***
#   RDA6       1  0.00349   5.5213  0.001 ***
#   RDA7       1  0.00251   3.9767  0.001 ***
#   RDA8       1  0.00128   2.0213  0.007 ** 
#   Residual 993  0.62771                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
test_margin <- anova (tbRDA, by = 'margin', parallel = 12)
test_margin

anova.cca(tbRDA, first = TRUE)
anova (tbRDA, permutations = 999)

Nresponse <- dplyr::select(dat,N.y)
Met <- dplyr::select(dat, mean_temp, mean_precip)
Spacetime<- dplyr::select(dat, longitude, latitude)
Light <- dplyr::select(dat,sum_canopy)
Ndep <- dplyr::select(dat,n_nh4, n_no3)

mod <- varpart(Nresponse,  Light, Ndep, scale = TRUE)
mod
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)

env2 <- dplyr::select(dat1, mean_temp, mean_precip, longitude, latitude,survey_year, sum_canopy,n_nh4, n_no3)
mod <- varpart(dat1$N.y, dat1$n_nh4, dat1$n_no3, dat1$latitude, dat1$L.x, scale = TRUE)
mod
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5, digits = 5)




