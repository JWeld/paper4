library(vegan)

#It has three matrices: spe (with species abundances per site), 
#env (with environmental variables per site), and spatial 
#(with longitude and latitude values per site; this matrix will be only used in partial CCA).

# Import data into R:

spe <- ordispe
env <- ordienv[,-(4:5)]
spatial <- ordienv[,(4:5)]
# Apply log+1 transformation to your species occurrences data (spe matrix)
#in order to correct for possible statistical errors associated to rare or very common species:

spelog <- decostand(spe, "log") %>% as.data.frame()
env <- as.data.frame(env)
## Perform CCA. We need to specify that spe (species distribution matrix) is explained by env (environmental matrix).

ccamodel <- cca(spe~., env)
summary(ccamodel)
## To perform a Partial CCA, we need to use a third matrix (conditioning matrix). To do that, we 
#must first combine variables 
#from the environmental ("env") and conditional ("spatial") matrices. We have to do that in order to later 
#apply another function ("ordistep" function). After we have combined all variables, we will apply the
#Partial CCA using a formula where we specify the response ("spe" matrix), the constraint variables 
#(each of the variables from "env" matrix), and the conditioning variables (variables from the 
#conditioning matrix; in our case "spatial" matrix)

envspatial<-as.data.frame(ordienv)

nams <- colnames(envspatial)

# partialccamodel <- formula(paste("spe ~", paste(nams[1: (length(envspatial)-(length(spatial)) )], collapse = " + "),
#                                  "+ Condition(", paste(nams[(length(envspatial)-(length(spatial)-1) ):length(envspatial)],
#                                                        collapse ="+"),")"))

partialccamodel<-cca(spe ~ NH4M + NO3M + SO4SM + latitude + longitude + PREC + survey_year, envspatial)

# However, we also have to automatically select variables of "env" matrix that best explain "spe" matrix. 
#We can do that by using a stepwise model from "ordistep" function. Let us do that with our "ccamodel" 
#(not the partial cca).

finalmodel<- ordistep(ccamodel)

# Then, we can calculate Variance Inflation Factors (VIF) for each of the constraints (variables) 
#from the "env" matrix (environmental matrix). 
#If we find an environmental variable with VIF>10, we'll know that this variable presents colinearity with
#another or other variables. In that case, we would have to delete the variable from our initial dataset 
#and redo all the analysis. In our example, no variable is redundant with each other (all of them have VIF<10).

vif.cca(finalmodel)

# Let us now fit the stepwise model (ordistep function) using our partial cca model ("partialccamodel"). 
#Note that we will use X (longitude) and Y (latitude) variables from the previously created "envspatial" 
#object as our conditioning variables. After "vif.cca", note that "X" (spatial) variable has a value > 10,
#so one should consider to delete the variable before the analysis.

simplemodel<-cca(partialccamodel, envspatial)

finalmodelpartial<- ordistep(simplemodel, scope=formula(partialccamodel))

vif.cca(finalmodelpartial)

### Ok, let us now call "ccamodel" object (not the "partialccamodel") to see how to interpret results.

ccamodel
# You will see this:

"                        Inertia   Proportion   Rank
Total                 1.2113     1.0000     
Constrained      0.8718     0.7197          6
Unconstrained  0.3395     0.2803         11
Inertia is mean squared contingency coefficient 

Eigenvalues for constrained axes:
  CCA1   CCA2   CCA3   CCA4   CCA5   CCA6 
0.5403 0.2243 0.0735 0.0196 0.0106 0.0036 

Eigenvalues for unconstrained axes:
    CA1     CA2     CA3     CA4     CA5     CA6     CA7     CA8     CA9    CA10    CA11 
0.13990 0.05623 0.03906 0.03192 0.02422 0.01940 0.01272 0.00742 0.00377 0.00265 0.00220 
"

## Note that "Total Inertia" is the total variance in species (observations matrix) distributions. 
#"Constrained Inertia" is the variance explained by the environmental variables (gradients matrix). 
#The "Proportion" values represent the percentages of variance of species distributions explained by
#Constrained (environmental) and Unconstrained variables. Eigenvalues of constrained and unconstrained axes 
#represent the amount of variance explained by each CCA axis (graphs usually present the first two 
#constrained axes, so take a look at their values).


# This is a critical step when doing the CCA. When we have our final model, 
#we must use permutation tests to observe if our whole CCA model, the CCA terms (environmental varibles),
#and CCA axes explain more variance of "spe" (observations) matrix than expected by chance
#(tests should be significant; p values lower or equal 0.05). If the tests are not significant,
#there is no point in using the CCA. In our example, we'll see that we can continue using the CCA results.

# Testing the significance of the CCA model:
anova.cca(finalmodel)

# Testing the significance of terms (environmental variables):
anova.cca(finalmodel, by="terms")

# Testing the significance of CCA axes (at least the first two or three should present a significant p value):
anova.cca(finalmodel, by="axis")

plot(finalmodel)

### Finally, we may want to generate graphs in order to better understand results.
#To do that, we have to take a look at "xlim" (limits of x axis), "ylim" (limits of y axis), 
#and "display" (if we want to observe species, environmental gradients, and/or sites in the graph) 
#arguments. For example, to show only species scores along CCA axes, use only "sp" inside display" argument.

plot(finalmodel, xlim=c(-1.5,2), ylim=c(-1,1.5), display=c("sp"))

# If you want to show species and environmental gradients in our plot, use both "sp" and "cn" inside "display" argument,

plot(finalmodel, xlim=c(-3,3), ylim=c(-3,3), display=c("sp","cn"))




#NMDS####
#DCA to check axis length for appropriate methods
ground_dca=decorana(spelog,iweigh = TRUE)
ground_dca
plot(ground_dca)

ground_NMDS=metaMDS(ordispe,k=3,trymax=100,autotransform=TRUE,expand=FALSE, plot=FALSE, parallel=6)
plot(ground_NMDS)
ground_NMDS #gives summary of methods used to generate ordination
stressplot(ground_NMDS)

#build a data frame with NMDS coordinates and metadata
MDS1 = ground_NMDS$points[,1]
MDS2 = ground_NMDS$points[,2]
gNMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Year = test$survey_year,ID_plot= test$ID_plot)
head(gNMDS)
#plot ellipses of year/refuge
ggplot(gNMDS, aes(x=MDS1, y=MDS2, col=ID_plot)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "gNMDS Plot")

ordiplot(ground_dca,type = "n")
orditorp(ground_dca, display = "species")
ordihull(ground_dca, test$survey_year, scaling= "symmetric", label = TRUE)
ordihull(ground_dca, test$ID_plot, scaling= "symmetric", label = TRUE)

ordiplot(ground_NMDS,type="n")
orditorp(ground_NMDS,display="species",col="red",air=0.01)
orditorp(ground_NMDS,display="sites",cex=1,air=0.01)
ordihull(ground_NMDS,groups=test$survey_year,draw="polygon",col="grey90",label=T)
ordihull(ground_NMDS,groups=test$ID_plot,draw="polygon",col="grey90",label=T)

#directional change?
plot(ground_NMDS, display = "sites", type="n")

site.scores<-scores(ground_NMDS, display="sites", choices=c(1,2))

s <- seq(length(site.scores[,1])-1)  # one shorter than data
arrows(site.scores[,1][s], site.scores[,2][s], site.scores[,1][s+1], site.scores[,2][s+1], col = 1:s)

###MultiSE#### AKA is sample size adequate? https://github.com/jslefche/multSE
#Asesses precision for dissimilarity-based multivariate analysis of ecological communities
# Run multSE function on a (bray curtis in this case) distance matrix
bray <- vegdist(ordispe, method="bray")
(output = multSE(bray, group = factor(test$ID_plot)))
print(output)
print(g2)
plot(output)
#plot output
ggplot(output, aes(x = n.samp, y = means, group = group)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.2)+
  geom_point(aes(shape = group, fill = group), size = 4) + 
  scale_shape_manual(values = c(21, 24:25), name = "")+
  scale_fill_manual(values = c("black","grey50","white"), name = "")+
  coord_cartesian(ylim = c(0, 0.65)) +
  theme_bw(base_size = 18) +
  labs(x = "Sample size (n)", y = "Multivariate pseudo SE") +
  theme(legend.position = c(0.8, 0.8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(output, aes(x = n.samp, y = means)) +
  geom_errorbar(aes(ymax = upper.ci, ymin = lower.ci), width = 0.2)+
  geom_point(size = 4) + 
  # coord_cartesian(ylim = c(0, 0.65)) +
  theme_bw(base_size = 18) +
  labs(x = "Sample size (n)", y = "MultSE based on residual mean square") +
  theme(legend.position = c(0.8, 0.8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


###BETADISPR and PERMATEST####
# To test if one or more groups is more variable than the others, ANOVA of the distances 
# to group centroids can be performed and parametric theory used to interpret the significance of F. 
# An alternative is to use a permutation test. permutest.betadisper permutes model residuals to 
# generate a permutation distribution of F under the Null hypothesis of no difference in dispersion 
# between groups.
# Pairwise comparisons of group mean dispersions can be performed by setting argument pairwise to TRUE.
# A classical t test is performed on the pairwise group dispersions. This is combined with a permutation
# test based on the t statistic calculated on pairwise group dispersions. An alternative to the classical
# comparison of group dispersions, is to calculate Tukey's Honest Significant Differences between groups,
# via TukeyHSD.betadisper.

## Bray-Curtis distances between samples
gdis <- vegdist(ordispe)

## group by year and refuge
ggroups <- factor(test$ID_plot)
ggroups2 <- factor(test$survey_year)

## Calculate multivariate dispersions
gmod <- betadisper(gdis, ggroups)
gmod
gmod2 <- betadisper(gdis, ggroups2)
gmod2

## Perform test
anova(gmod)#refuges sig 0.005
anova(gmod2)#years sig 0.017

## Permutation test
permutest(gmod, pairwise = TRUE)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq     F N.Perm Pr(>F)    
# Groups      13  7.2586 0.55836 24.93    999  0.001 ***
#   Residuals 1408 31.5352 0.02240                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

permutest(gmod2, pairwise = TRUE)


## Tukey's Honest Significant Differences
(gmod.HSD <- TukeyHSD(gmod))
plot(gmod.HSD)

(gmod2.HSD <- TukeyHSD(gmod2))
plot(gmod2.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(gmod, main="multivariate dispersion, ground layer, refuge")
plot(gmod2, main="multivariate dispersion, ground layer, year")
legend(0.5, 0.3, legend=c("2006", "2011", "2016"),
       col=c("black", "red", "green"), lty=1, cex=0.8)


## Draw a boxplot of the distances to centroid for each group
boxplot(gmod)

library(labdsv)
indic <- indval(as.data.frame(ordispe), as.numeric(as.factor(test$ID_plot)))
summary(indic, short = TRUE)
summary(indic)
indic$maxcls
summary(indval(ordispe, test$ID_plot))

data(bryceveg) # returns a vegetation data.frame
data(brycesite)
clust <- cut(brycesite$elev,5,labels=FALSE)
summary(indval(bryceveg,clust))


dune <- as.data.frame(ordispe_h)
dune.env <- as.data.frame(ordienv)
### Dune data
data(dune)
data(dune.env)
mod0 <- rda(dune ~ 1, dune.env)  # Model with intercept only
mod1 <- rda(dune ~ ., dune.env)  # Model with all explanatory variables

## With scope present, the default direction is "both"
mod <- ordistep(mod0, scope = formula(mod1))
mod
## summary table of steps
mod$anova
ordiplot(mod)
ordiplot(mod, display = "sites")
## Example of ordistep, forward
ordistep(mod0, scope = formula(mod1), direction="forward")

## Example of ordiR2step (always forward)
## stops because R2 of 'mod1' exceeded
ordiR2step(mod0, mod1)
