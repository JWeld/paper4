library(tidyverse)
library(vegan)
options(mc.cores = parallel::detectCores())

#Intensive plots (VG)####
#Data#
ordispe
ordispe_sub #matrix of subplots
ordi_extras #factors columns removed from subplots matrix
ordispe_pl #matrix of plots
ordispe_pl_extras #factors columns removed from plots matrix
ordienv #factor columns removed from subplots matrix plus env variables
#ordispe_h <- decostand(ordispe_sub, "hellinger") #hellinger transform to make usable with PCA
ordispe_h <- decostand(ordispe, "hellinger") #hellinger transform to make usable with PCA

ordispe_h2 <-decostand(ordispe_pl, "hellinger") 

#Ordinations####

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
PCA_gh <- rda(ordispe_h, scale = F)

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
             label = T, col = "red")

table(VG_g_wide$ID_subplot)
table(VG_g_wide$ID_plot)


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
ordispe_h2 <- decostand(ordispe_pl, "hellinger")

#PCA
#on hell transformed data
PCA_gh2 <- rda(ordispe_h2)
pl <- ordiplot(PCA_gh2, display = "sites")
ordispider (pl, groups = ordispe_pl_extras$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh2, display = "species")
points(pl, "sites", pch=21, col="grey60", bg="grey80")

text(pl, "sites", col="blue", cex=0.6)
text(pl, "species", col="blue", cex=0.6)

#extract scores
PCA_scores <-  vegan::scores(PCA_gh)
PCA_scores_sites <- as.data.frame(PCA_scores$sites)
PCA_scores_species <- as.data.frame(PCA_scores$species)

PCA_scores_sites <- PCA_scores_sites %>% rownames_to_column(.,var = "ID_fine2")
PCA_scores_species <- PCA_scores_species %>% rownames_to_column(.,var = "Species")

#And at plot level
PCA_scores2 <-  vegan::scores(PCA_gh2)
PCA_scores_sites2 <- as.data.frame(PCA_scores2$sites)
PCA_scores_species2 <- as.data.frame(PCA_scores2$species)

PCA_scores_sites2 <- PCA_scores_sites2 %>% rownames_to_column(.,var = "ID_fine")
PCA_scores_species2 <- PCA_scores_species2 %>% rownames_to_column(.,var = "Species")


PCA_scores_sites <- left_join(PCA_scores_sites, extras, by = "ID_fine2")
#combine PCA scores with other data####

combo_sites <- left_join(PCA_scores_sites2, select(extras, -c(ID_fine2, ID_subplot)),
                         by = "ID_fine") %>% distinct()
#year as index from first year
combo_sites$year.i <- I(combo_sites$survey_year-1997)

all_data <- left_join(PCA_scores_sites, IM_env, by = "ID")
#discard where deposition data is available
all_data <- filter(all_data, NH4M > 0)

dist <- vegdist(ordispe_sub,  method = "bray")

library(corrplot)
all_data %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()

#Add PC movement from base and sequential distances to ####
#both subplot and plot level data
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


#Add distance to baseline as response variable to all_data
combosites5 <- left_join(select(combosites4, ID_fine, ID, ID_site, ID_plot, PC1, PC2, survey_year,
                                pc_dist, pc_dist_base),
                         IM_env, by = "ID") %>% ungroup()

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
#combosites6 <- drop_na(combosites6, pc_dist)

#Plot level final data #####
plot.level.dat <- combosites6
plot.level.dat$year.i <- I(plot.level.dat$survey_year -1997)


#subplotlevel PCA movements
all_data2 <- all_data %>% mutate(NTOT = NH4M + NO3M)

combosubs <- all_data2 %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  mutate(pc1diff = PC1 - lag(PC1)) %>% mutate(pc2diff = PC2 - lag(PC2)) %>% 
  mutate(pc_dist = sqrt((pc1diff)^2 + (pc2diff)^2))

ggplot(combosubs %>% drop_na(pc_dist) %>% drop_na(ID_plot), aes(ID_plot, pc_dist, colour = ID_plot)) +
  geom_boxplot() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")


temp2 <- combosubs %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  filter(row_number()==1) %>% select(ID_subplot, PC1, PC2) %>%
  rename(basePC1 = PC1, basePC2 = PC2)

combosubs2 <- left_join(combosubs, temp2, by = "ID_subplot")
combosubs3 <- combosubs2 %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  mutate(pc1diffbase = PC1 - basePC1) %>% mutate(pc2diffbase = PC2 - basePC2) %>% 
  mutate(pc_dist_base = sqrt((pc1diffbase)^2 + (pc2diffbase)^2))

#create lag terms for atmospheric pollutants
combosubs3 <- combosubs3 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(NTOT.lag = lag(NTOT))
combosubs3 <- combosubs3 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(NH4M.lag = lag(NH4M))
combosubs3 <- combosubs3 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(NO3M.lag = lag(NO3M))
combosubs3 <- combosubs3 %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(SO4SM.lag = lag(SO4SM))
combosubs3 <- ungroup(combosubs3)
combosubs3 <- drop_na(combosubs3, pc_dist)


ggplot(combosubs3 %>%
         drop_na(ID_plot), aes(SO4SM.lag, pc_dist_base, colour = ID_plot)) +
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)
#geom_line(size=1) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#Subplot level final data #####
subplot.level.dat <- combosubs3

#PLOTS- cup and ball models####
ggplot(subplot.level.dat, aes(PC1, PC2, colour = ID_plot)) + geom_point() +
  facet_wrap(facets = "ID_plot")

ggplot(plot.level.dat, aes(PC1, PC2, colour = ID_plot)) + geom_point() +
  facet_wrap(facets = "ID_plot")

ggplot(subplot.level.dat, aes(PC1, PC2,colour = ID_plot)) + geom_point()

ggplot(plot.level.dat, aes(PC1, PC2,colour = ID_plot)) + geom_point()

#Figure####
ggplot(drop_na(plot.level.dat), aes(PC1, PC2, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow()) +
  theme_light()


ggplot((plot.level.dat), aes(PC1, PC2, colour = ID_plot)) + geom_point()+ 
  geom_path(linewidth=1.25, arrow=arrow()) 

test1 <- subplot.level.dat %>% group_by(ID_fine) %>% mutate(PC1m = mean(PC1)) %>% 
  mutate(PC2m = mean(PC2)) %>% ungroup()
ggplot(test1, aes(PC1m, PC2m, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow())
ggplot(subplot.level.dat, aes(PC1, PC2, colour = ID_subplot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow())+ theme(legend.position="none")

#Figure - cumulative distance moved from baseline
ggplot(plot.level.dat %>% drop_na(pc_dist_base) %>%
         drop_na(ID_plot), aes(survey_year, pc_dist_base, colour = ID_plot)) +
  #  geom_point()+
  #  geom_smooth(method = "lm")+
  geom_line(size=1) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
  )

