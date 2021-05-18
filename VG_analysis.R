library(tidyverse)
library(vegan)
options(mc.cores = parallel::detectCores())


#Intensive plots (VG)####
# hell_VG_fm <- decostand(VG_fm, "hellinger") 
# hell_VG_bm <- decostand(VG_bm, "hellinger") 
ordispe_h <- decostand(ordispe, "hellinger") #hellinger transform to make usable with PCA
#DCA
#DCA_f <- decorana(VG_fm, iweigh = 0)
#DCA_b <- decorana(VG_bm, iweigh = 0)
DCA_g <- decorana(ordispe, iweigh = 0)

#pl <- ordiplot(DCA_f)
#ordispider (pl, groups = VG_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

#pl <- ordiplot(DCA_b)
#ordispider (pl, groups = VG_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_g)
ordispider (pl, groups = test$ID_plot, col = c("red","blue","green","grey","pink","yellow"),
            label = F)
ordiellipse (pl, groups = test$ID,
             label = F)
#PCA####
PCA_g <- rda(ordispe)
pl <- ordiplot(PCA_g, display = "sites")
#text(pl, "sites", col="blue", cex=0.5)
ordispider (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(PCA_g, display = "species")
text(pl, "species", col="blue", cex=0.5)


#on hell transformed data
#PCA_fh <- rda(hell_VG_fm, scale = FALSE)
#PCA_bh <- rda(hell_VG_bm, scale = FALSE)
PCA_gh <- rda(ordispe_h, scale = FALSE)
#pl <- ordiplot(PCA_gh, display = "sites")
#ordihull(pl, VG_g_wide$AreaCode, label = TRUE) 
#pl <- ordiplot(PCA_fh, display = "sites")
#ordispider (pl, groups = VG_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
#pl <- ordiplot(PCA_bh, display = "sites")
#ordispider (pl, groups = VG_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "sites")
ordispider (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "sites")
ordiellipse (pl, groups = test$ID_plot,
             label = F)
#by plot
pl <- ordiplot(PCA_gh, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = test$ID_plot,
             label = F, col = "red")
#text(pl, "sites", col="blue", cex=0.9)
#by site
pl <- ordiplot(PCA_gh, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = test$ID_site,
             label = F, col = "red")

#barplot(as.vector(PCA_gh$CA$eig)/sum(PCA_gh$CA$eig)) 
#sum((as.vector(PCA_gh$CA$eig)/sum(PCA_gh$CA$eig))[1:3])

table(VG_g_wide$ID_subplot)
table(VG_g_wide$ID_plot)

#CCA####
test.cca <- cca(ordispe ~ ordienv)
test.cca
plot(test.cca)

#RDA
test.rda <- rda(ordispe_h ~ ordienv)
test.rda
plot(test.rda)

#Intensive plots per intensive/year####
ordispe_h2 <- decostand(VG_gm2, "hellinger")

#DCA
# DCA_g2 <- decorana(VG_gm2, iweigh = 1)
# pl <- ordiplot(DCA_g2)
# ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
# ordiellipse (pl, groups = VG_g_wide2$ID,
#              label = F)
#PCA
#on hell transformed data
 PCA_gh2 <- rda(ordispe_h2)
 pl <- ordiplot(PCA_gh2, display = "sites")
 ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
 pl <- ordiplot(PCA_gh2, display = "species")
 ordiellipse (pl, groups = VG_g_wide2$ID_site,
              label = F)
 text(pl, "sites", col="blue", cex=0.6)
text(pl, "species", col="blue", cex=0.6)

#extract scores
PCA_scores <- scores(PCA_gh)
PCA_scores_sites <- as.data.frame(PCA_scores$sites)
PCA_scores_species <- as.data.frame(PCA_scores$species)

PCA_scores_sites <- PCA_scores_sites %>% rownames_to_column(.,var = "ID_fine2")
PCA_scores_species <- PCA_scores_species %>% rownames_to_column(.,var = "Species")

#And at plot level
PCA_scores2 <- scores(PCA_gh2)
PCA_scores_sites2 <- as.data.frame(PCA_scores2$sites)
PCA_scores_species2 <- as.data.frame(PCA_scores2$species)

PCA_scores_sites2 <- PCA_scores_sites2 %>% rownames_to_column(.,var = "ID_fine")
PCA_scores_species2 <- PCA_scores_species2 %>% rownames_to_column(.,var = "Species")
combo_sites <- left_join(PCA_scores_sites2, select(extras, -c(ID_fine2, ID_subplot)),
                         by = "ID_fine") %>% distinct()


extras <- test %>% 
  select(ID, ID_site,ID_fine, ID_fine2, ID_plot, ID_subplot,survey_year)
PCA_scores_sites <- left_join(PCA_scores_sites, extras, by = "ID_fine2")

#combine PCA scores with other data####
all_data <- left_join(PCA_scores_sites, IM_env, by = "ID")
#discard where deposition data is available
all_data <- filter(all_data, NH4M > 0)
#create by plot co-efficient of variation.
#use cv function in raster package, manual seems to fail with group_by?
#hogs select from dplyr though...
#library(raster)

#cv <- sd(data) / mean(data) * 100

# all_data <- all_data %>% group_by(ID_plot) %>%
#   mutate(CV1 = cv(PC1)) %>% mutate(CV2 = cv(PC2)) %>% 
#   ungroup()

#betadisper, nestedness and turnover####
#bray_curtis####
dist <- vegdist(ordispe,  method = "bray")
disps <- betadisper(dist, test$ID_fine)
disps_df <- tibble::enframe(disps$distances) %>%
  rename(ID_fine2 = name, dispersion = value)

plot(disps)
anova(disps)

#add distance to centroids to all_data
all_data <- left_join(all_data, disps_df, by = "ID_fine2")
#group mean distance to centroid
all_data <- all_data %>% group_by(ID_fine) %>% mutate(ID_fine_disp = mean(dispersion))
str(all_data)
all_data$ID_fine2 <-as.factor(all_data$ID_fine2)
all_data$ID <-as.factor(all_data$ID)
all_data$year.i <- I(dat$survey_year - 1998)
#all_data$survey_year <-as.factor(all_data$survey_year)
#all_data$year.i <- as.factor(all_data$year.i)


#BETAPART STUFF###############################################################
library(betapart)
library(tidyverse)
library(vegan)
#Use sign function to convert abundance data 
#into prescence abscence matrix that can be fed to Betapart. Can't include factors (non numeric)
test_part <- test %>% group_by(ID_plot) %>% mutate(first= min(survey_year)) %>% 
  mutate(last=max(survey_year)) %>% select(ID_fine2, first, last, everything()) %>% ungroup()
#convert to PA for betapart
test_part2 <- as.data.frame(test_part) %>% select(-c(1:9, 333:339)) %>% sign() %>% 
  cbind(test_part$ID_fine) %>% 
  rename(ID_fine = 'test_part$ID_fine') %>% select(ID_fine, everything())
#make a list of matrices, one for each plot/year combo
test_list <- test_part2 %>% 
  split(.$ID_fine) %>% lapply(., function(x) { x["ID_fine"] <- NULL; x }) %>% 
  lapply(.,function(x) {as.matrix(x)})

#apply betapart core to create a betapart object for each matrix
core_list <- lapply(test_list,function(x) betapart.core(x))

multi_list <- lapply(core_list,function(x) beta.multi(x))

multi_df <- data.frame(matrix(unlist(multi_list), nrow=85, byrow=TRUE),stringsAsFactors=FALSE) %>% 
  cbind(names(multi_list))

multi_df <- multi_df %>% rename(beta.SIM = X1, beta.SNE = X2, beta.SOR = X3, ID_fine = 'names(multi_list)')

#LT plots report only one value for each year, not possible to calculate beta 
#diversity, and betapart returns NaN

all_data <- left_join(all_data, multi_df)


## Matrix temperature
# out <- nestedtemp(ordispe)
# out
# plot(out)
# plot(out, kind="incid")
## Use oecosimu to assess the non-randomness of checker board units
#nestedchecker(ordispe)
#oecosimu(ordispe, nestedchecker, "quasiswap")
## Another Null model and standardized checkerboard score
#oecosimu(ordispe, nestedchecker, "r00", statistic = "C.score")

#library(adespatial)

library(corrplot)
all_data %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()


#simple multiple regression####
ml <- lm(dispersion  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           PREC + survey_year, data = all_data)
summary(ml)

ml <- lm(beta.SOR  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           PREC + survey_year, data = all_data)
summary(ml)

glm <- glm(dispersion  ~  NH4M + NO3M + SO4SM + TEMP + latitude + longitude +
  PREC + survey_year, family = gaussian, data = all_data)

summary(glm)
#seems more going on with PC2 as response
library(visreg)
visreg(ml, "NO3M")
visreg(ml, "NH4M")
visreg(ml, "SO4SM")
visreg(ml, "NO3M", by = "survey_year")
visreg(ml, "NH4M", by = "survey_year")
visreg(ml, "NO3M", by = "year_scaled")
visreg(ml, "NH4M", by = "year_scaled")
visreg(ml, "NH4M", by = "NO3M")
visreg(ml, "NO3M", by = "latitude")
visreg(ml, "SO4SM", by = "survey_year")
visreg(ml, "SO4SM", by = "TEMP")


#Nmle AR regressions####

#index year! ??
library(nlme)

summary(mod<-lme(ID_fine_disp~NH4M + NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_fine,
                 #correlation=corCAR1(form=~survey_year|ID_fine),
                 data=all_data))

summary(mod<-lme(beta.SOR~NH4M + NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_fine,
                 #correlation=corCAR1(form=~survey_year|ID_fine),
                 data=drop_na(all_data)))

summary(mod<-lme(ID_fine_disp~survey_year,
                 random=~1|ID_fine,
                 correlation=corCAR1(form=~survey_year|ID_fine),
                 data=all_data))
anova(mod)
fixed.effects(mod)

plot(mod, col = as.numeric(factor(all_data$survey_year,
                                  levels = unique(all_data$survey_year))),
     pch = 16, main = "resids")


#NMDS####
# First step is to calculate a distance matrix. See PCOA for more information about the distance measures
# Here we use bray-curtis distance, which is recommended for abundance data

# In this part, we define a function NMDS.scree() that automatically 
# performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
# NMDS.scree <- function(x) { #where x is the name of the data frame variable
#   plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress),
#        xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions",
#        ylab = "Stress", main = "NMDS stress plot")
#   for (i in 1:10) {
#     points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
#   }
# }

# Use the function that we just defined to choose the optimal nr of dimensions
#NMDS.scree(dist)

# Here, we perform the final analysis and check the result
NMDS1 <- metaMDS(dist, k = 2, trymax = 50, trace = T, parallel = 6)
NMDS1
ordiplot(NMDS1)
ordihull (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"),
          label = T)

NMDS2 <- metaMDS(dist, k = 4, trymax = 50, trace = T, parallel = 6)
NMDS2

mod <- NMDS2
stressplot(mod)
gof <- goodness(mod)
gof
plot(mod, display = "sites", type = "n")
points(mod, display = "sites", cex = 2*gof/mean(gof))

pl <- plot(NMDS2)
ordihull (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
ordiellipse (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

#keep only sites with observations over at least three years
n.obs.dat <-as.data.frame(table(all_data$ID_site))
n.obs.dat <-as.data.frame(table(all_data$ID_plot))

# drop.list <- filter(n.obs.dat, Freq <3)
# dat <- filter(all_data2, ID_site %!in% drop.list$Var1)
# dat2 <- dat

#plots- cup and ball models####
ggplot(all_data, aes(PC1, PC2)) + geom_point() + facet_wrap(facets = "ID_plot")

ggplot(all_data, aes(PC1, PC2)) + geom_point() + geom_path() +
  facet_wrap(facets = "ID_plot")

ggplot(all_data, aes(PC1, PC2, colour = ID_plot, label = ID_plot)) + geom_point() + 
  geom_path() + geom_text()

ggplot(all_data, aes(PC1, PC2,colour = ID_plot)) + geom_point()

ggplot(all_data, aes(PC1, PC2,colour = ID_plot)) + geom_point()

ggplot(all_data, aes(PC1, PC2, colour = ID_plot)) + geom_point() + 
  geom_path(arrow = arrow())
#facet_wrap(facets = "ID_plot")

ggplot(drop_na(combo_sites), aes(PC1, PC2, colour = ID_plot)) + geom_point()+ 
  geom_path() 


#GAM
#try same model with bam
# N3_final_b <- bam(Ellenberg_N ~ 
#                    #s(n_nh4, bs="tp")+
#                    #s(n_no3, bs="tp")+
#                    #ti(n_nh4, n_no3, bs="tp")+
#                    te(NH4, NO3, bs = "tp") +
#                    s(canopy,bs="tp") +
#                    #s(L.x, bs= "tp")+
#                    #s(herb_layer.x, bs="tp")+
#                    #s(temperature,bs= "tp")+
#                    s(precipitation,bs="tp")+
#                    s(age, bs="tp", k=5)+
#                    #as.factor(age)+ 
#                    as.factor(tree)+
#                    s(survey_year,ID_site, bs="re")+
#                    s(ID_site, bs="re"),
#                    #te(longitude,latitude,survey_year,bs="tp",d=c(2,1)),argGam = list(select=T),
#                  data = dat_names)
# 
# b <- getViz(N3_final_b)
# summary(b)
# draw(b, parametric = FALSE, residuals = TRUE)
# appraise(b)
# acf(residuals(b))
# #show sig areas only
# plot(sm(b, 1))+l_fitRaster(pTrans = function(.p).p<=0.05)
  