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

summary(PCA_gh)
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
             label = T, col = "red")
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

#RDA####
test.rda <- rda(ordispe_h ~ ordienv)
test.rda
plot(test.rda)

mod <- ordistep(test.rda)
mod0 <- rda(ordispe_h ~ 1, ordienv)  # Model with intercept only
mod1 <- rda(ordispe_h ~ ordienv)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
## summary table of steps
mod$anova

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


extras <- test %>% 
  select(ID, ID_site,ID_fine, ID_fine2, ID_plot, ID_subplot,survey_year)
PCA_scores_sites <- left_join(PCA_scores_sites, extras, by = "ID_fine2")

combo_sites <- left_join(PCA_scores_sites2, select(extras, -c(ID_fine2, ID_subplot)),
                         by = "ID_fine") %>% distinct()
#year as index from first year
combo_sites$year.i <- I(combo_sites$survey_year-1997)

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
all_data$year.i <- I(all_data$survey_year -1997)
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

#X Dispersion, betapart reponses X####
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
# NMDS1 <- metaMDS(dist, k = 2, trymax = 50, trace = T, parallel = 6)
# NMDS1
# ordiplot(NMDS1)
# ordihull (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"),
#           label = T)
# 
# NMDS2 <- metaMDS(dist, k = 4, trymax = 50, trace = T, parallel = 6)
# NMDS2

# mod <- NMDS2
# stressplot(mod)
# gof <- goodness(mod)
# gof
# plot(mod, display = "sites", type = "n")
# points(mod, display = "sites", cex = 2*gof/mean(gof))
# 
# pl <- plot(NMDS2)
# ordihull (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
# ordiellipse (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

#keep only sites with observations over at least three years
#n.obs.dat <-as.data.frame(table(all_data$ID_site))
#n.obs.dat <-as.data.frame(table(all_data$ID_plot))

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

ggplot(drop_na(combo_sites), aes(PC1, PC2, colour = ID_plot, label = ID_fine)) + geom_point()+ 
  geom_path() + geom_label()

ggplot(drop_na(combo_sites), aes(PC1, PC2, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow())


#Add distances between sequential observations and distance to base-line for
#all to allow plots like in Lamothe 2019

combosites2 <- combo_sites %>% ungroup() %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(pc1diff = PC1 - lag(PC1)) %>% mutate(pc2diff = PC2 - lag(PC2)) %>% 
  mutate(pc_dist = sqrt((pc1diff)^2 + (pc2diff)^2))

ggplot(combosites2 %>% drop_na(pc_dist) %>% drop_na(ID_plot), aes(ID_plot, pc_dist, colour = ID_plot)) +
  geom_boxplot() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")


temp <- combosites2 %>% ungroup() %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  filter(row_number()==1) %>% select(ID_plot, PC1, PC2) %>% rename(basePC1 = PC1, basePC2 = PC2)

combosites3 <- left_join(combosites2, temp, by = "ID_plot")
combosites4 <- combosites3 %>% ungroup() %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(pc1diffbase = PC1 - basePC1) %>% mutate(pc2diffbase = PC2 - basePC2) %>% 
  mutate(pc_dist_base = sqrt((pc1diffbase)^2 + (pc2diffbase)^2))


ggplot(combosites4 %>% drop_na(pc_dist_base) %>%
         drop_na(ID_plot), aes(survey_year, pc_dist_base, colour = ID_plot)) +
#  geom_point()+
#  geom_smooth(method = "lm")+
  geom_line(size=1) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
                               )
#Add distance to baseline as response variable to all_data
combosites5 <- left_join(select(combosites4, ID_fine, ID, ID_site, ID_plot, survey_year,
                                pc_dist, pc_dist_base),
          IM_env, by = "ID") %>% ungroup()
#Add betapart data
t1 <- select(combosites5, ID,pc_dist, pc_dist_base) %>% drop_na(ID)
all_data2 <- left_join(all_data, t1) %>% ungroup()
all_data2 <- all_data2 %>% mutate(NTOT = NH4M + NO3M)

combosites5 <- left_join(combosites5, multi_df, by = "ID_fine")

combosites5 <- combosites5 %>% drop_na(ID) %>% ungroup()
combosites5 <- combosites5 %>% mutate(NTOT = NH4M + NO3M)

#create lag terms for atmospheric pollutants
combosites6 <- combosites5 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(NTOT.lag = lag(NTOT))
combosites6 <- combosites6 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(NH4M.lag = lag(NH4M))
combosites6 <- combosites6 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(NO3M.lag = lag(NO3M))
combosites6 <- combosites6 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(SO4SM.lag = lag(SO4SM))
combosites6 <- ungroup(combosites6)
combosites6 <- drop_na(combosites6, pc_dist)


#subplotlevel PCA movements
combosubs <- all_data2 %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  mutate(pc1diff = PC1 - lag(PC1)) %>% mutate(pc2diff = PC2 - lag(PC2)) %>% 
  mutate(pc_dist = sqrt((pc1diff)^2 + (pc2diff)^2))

ggplot(combosubs %>% drop_na(pc_dist) %>% drop_na(ID_plot), aes(ID_plot, pc_dist, colour = ID_plot)) +
  geom_boxplot() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")


temp2 <- combosubs %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  filter(row_number()==1) %>% select(ID_subplot, PC1, PC2) %>% rename(basePC1 = PC1, basePC2 = PC2)

combosubs2 <- left_join(combosubs, temp2, by = "ID_subplot")
combosubs3 <- combosubs2 %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  mutate(pc1diffbase = PC1 - basePC1) %>% mutate(pc2diffbase = PC2 - basePC2) %>% 
  mutate(pc_dist_base = sqrt((pc1diffbase)^2 + (pc2diffbase)^2))


ggplot(combosubs3 %>% drop_na(pc_dist_base) %>%
         drop_na(ID_plot), aes(survey_year, pc_dist_base, colour = ID_plot)) +
   geom_point()+
    geom_smooth(method = "lm", se=FALSE)
  #geom_line(size=1) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#X PCA distance to base reponses X####
#combine PCA base distance and env variables and model
library(corrplot)
combosites5 %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()
combosites5 %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() 

#simple multiple regression####
ml <- lm(pc_dist_base  ~  NTOT + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           #beta.SOR +
           PREC +
           survey_year, data = combosites5)
summary(ml)
library(car)
vif(ml)
#attributes(alias(ml)$Complete)$dimnames[[1]] #beta.SOR is linearly dependant on SOR + SNE

ml <- lm(pc_dist_base  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           #beta.SOR +
           PREC + survey_year, data = all_data2)
summary(ml)

ml <- lm(pc_dist_base  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           beta.SOR+
           PREC + survey_year, data = combosites5)
summary(ml)


glm <- glm(pc_dist_base  ~  NH4M + NO3M + SO4SM + TEMP + latitude + longitude +
             PREC + survey_year, family = gaussian, data = combosites5)

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

#index year?
library(nlme)
library(lme4)
combosites5$year.i <- I(combosites5$survey_year-1997)
combosites5$year.i <- as.numeric(combosites5$year.i)

barplot(table(combosites5$ID_plot))
barplot(table(combosites5$latitude))
barplot(combosites5$NTOT)
barplot(combosites5$year.i)
barplot(table(combosites5$year.i))
hist(combosites5$pc_dist_base)#right skewed
hist(log(combosites5$pc_dist_base))#log transf?

summary(nmod1 <- lme(pc_dist_base~
                       NH4M +
                       NO3M +
                       SO4SM +
                       latitude + longitude,
                 random=~1|ID_plot,
                 correlation=corCAR1(form=~survey_year|ID_plot),
                 data=combosites5))

summary(nmod2 <- lme(pc_dist_base~
                   #NH4M * NO3M +
                   NTOT + 
                   SO4SM +
                   latitude + longitude +
                   survey_year,
                 random=~1|ID_plot,
                 correlation=corCAR1(form=~survey_year|ID_plot),
                 data=combosites5))

summary(nmod3 <- lme(log(pc_dist_base+1)~
                       #NH4M * NO3M +
                       NTOT + 
                       SO4SM +
                       latitude + longitude +
                       survey_year,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~survey_year|ID_plot),
                     data=combosites5))
#basic LM
summary(mod<-lm(pc_dist_base~
                   #NH4M * NO3M +
                   NTOT + 
                   SO4SM +
                   latitude + longitude +
                   survey_year,
                 data=combosites5))

par(mfrow = c(2,2))
plot(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))
vif(mod)#car

#GLM null model
GLM <- gls(log(pc_dist_base+1)~
             #NH4M * NO3M +
             NTOT + 
             SO4SM +
             latitude + longitude +
             survey_year,
           data=combosites5,
           method = "ML")
summary(GLM)
plot(GLM)

#LME mixed
lmm0 <- lmer(pc_dist_base~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year+
               (1|ID_plot),
             data=combosites5)

lmm1 <- lmer(log(pc_dist_base+1)~
              #NH4M * NO3M +
              NTOT + 
              SO4SM +
              latitude +
              longitude +
              survey_year+
    (1|ID_plot),
    data=combosites5)

lmm1.5 <- lmer(pc_dist_base~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year+
               (1|ID_plot),
             data=combosites5)

lmm2 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year+
               (1|ID),
             data=combosites5)

lmm3 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               #year.i+
               (1|ID),
             data=combosites5)

lmm4 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year +
               (1|ID_plot)+
               (1|ID),
             data=combosites5)

AIC(GLM, lmm1, lmm1.5, lmm2,lmm3)

plot(ranef(lmm4)) 
# QQ plots (drawn to the same scale!)
par(mfrow = c(1,2))
lims <- c(-3.5,3.5)
qqnorm(resid(GLM, type = "pearson"),
       xlim = lims, ylim = lims,main = "GLM")
abline(0,1, col = "red", lty = 2)
qqnorm(resid(lmm1, type = "pearson"),
       xlim = lims, ylim = lims, main = "lmm1")
abline(0,1, col = "red", lty = 2)

mod <- nmod2

summary(lmm1)
#percentage variance explained by random factor ID_plot
0.0007465/(0.0007465 + 0.0026378)#22%

summary(lmm1.5)
#percentage variance explained by random factor ID_plot
0.0009634/(0.0009634 + 0.0035281)#22%

summary(lmm2)
#percentage variance explained by random factor ID
0.002088/(0.002088 + 0.001063)#66%

anova(lmm1, lmm2)

qqnorm(resid(mod))#2016_SE14_1 is an outlier on 0.20+
qqline(resid(mod))

#Facet plots of model predictions####
#comboplots
(mm_plot <- ggplot(combosites5, aes(x = survey_year, y = pc_dist_base, colour = ID_plot)) +
    facet_wrap(~ID_plot, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(combosites5, pred = predict(gam1)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

#all_data2
(mm_plot <- ggplot(all_data2, aes(x = survey_year, y = pc_dist_base, colour = ID_plot)) +
    facet_wrap(~ID_plot, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(all_data2, pred = predict(gam1)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)


summary(mod<-lme(pc_dist_base~NH4M * NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_subplot,
                 correlation=corCAR1(form=~survey_year|ID_subplot),
                 data=drop_na(all_data2)))


summary(mod<-lme(Sensitivity~time+sp,random=~1|Plot/Tree,
                 correlation=corCAR1(form=~time|Plot/Tree), data=data_Gd))
anova(Gd.mod3.s)
#GAM
#try same model with bam
k=5
bs="ts"
gam1 <- bam(pc_dist_base ~ #s(NH4M) +
             #s(NO3M)+
             #s(SO4SM,k=k, bs=bs)+
             #s(NTOT,k=k, bs=bs)+
             #s(PREC)+
             s(latitude,k=k, bs=bs)+
             s(longitude,k=k, bs=bs)+
             #te(latitude, longitude,k=k, bs=bs)+
             s(survey_year,k=k, bs=bs)+
             s(ID_plot, bs="re"),
             nthreads=12, 
             family = gaussian,
             data = all_data2)

b <- getViz(gam1)
summary(b)
draw(b, parametric = FALSE, residuals = TRUE)
appraise(b)
acf(residuals(gam))
concurvity(b)

library(dsm)#vis.concurvity

vis.concurvity(b, type = "estimate")
vis.concurvity(b, type = "worst")



#X PCA sequential distance reponses X####

#lm
summary(lm1 <- lm(pc_dist~
                    NTOT.lag +
                       #NH4M +
                       #NO3M +
                       SO4SM +
                       latitude +
                       longitude +
                       survey_year, 
                     data=combosites6))

#nmle
hist(combosites6$pc_dist)
hist(log(combosites6$pc_dist))#log transf?


summary(nmod1 <- lme(pc_dist~
                       NTOT.lag +
                       #NH4M +
                       #NO3M +
                       SO4SM +
                       latitude + longitude,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~survey_year|ID_plot),
                     data=combosites6))



#MAPS#####

library(ggplot2)  # ggplot() fortify()
library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)  # readOGR() spTransform()
library(raster)  # intersect()
library(ggsn)  # north2() scalebar()
library(rworldmap)
world <- getMap(resolution = "low")

(site_check <- ggplot(all_data, mapping = aes(x = longitude, y = latitude)) + 
    geom_point(alpha = 0.5))

clipper_europe <- as(extent(3, 32, 50, 72), "SpatialPolygons")

proj4string(clipper_europe) <- CRS(proj4string(world))

world_clip <- raster::intersect(world, clipper_europe)

world_clip_f <- fortify(world_clip)

(site_map <- ggplot() + 
    geom_polygon(data = world_clip_f, 
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "black") + 
    geom_point(colour = "black", size = 2.5, shape = 21, fill = "grey70",
               stroke = 1,
               aes(x = longitude, y = latitude),
               data = all_data) +
    # labs(shape="Forest type")+
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") + 
    coord_quickmap())

(site_map <- ggplot() + 
    geom_polygon(data = world_clip_f, 
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "black") + 
    geom_point(colour = "black", size = 1.5,
               aes(x = longitude, y = latitude),
               data = all_data) +
    # labs(shape="Forest type")+
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") + 
    coord_quickmap())


#rotation forest
# library(rotationForest)
# data(iris)
# y <- as.factor(ifelse(iris$Species[1:100]=="setosa",0,1))
# x <- iris[1:100,-5]
# rF <- rotationForest(x,y)
# predict(object=rF,newdata=x)
# summary(rF)
# plot(rF$loadings)
