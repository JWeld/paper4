library(tidyverse)
library(vegan)
options(mc.cores = parallel::detectCores())
library(gvlma)#checks assumptions for models
library(sjPlot)
library(mgcv)
library(mgcViz)
library(gratia)
#Intensive plots (VG)####
#Data#
ordispe_sub #matrix of subplots
ordispe_pl #matrix of plots
ordispe_pl_extras #factors columns removed from plots matrix
ordienv #factor columns removed from subplots matrix plus env variables
ordispe_h <- decostand(ordispe_sub, "hellinger") #hellinger transform to make usable with PCA
#ordispe_h2 <-decostand(ordispe_pl) 

#Exclude sites with 
#Ordinations####
#NMDS####
#subplots
nmds_sub <- metaMDS(ordispe_sub, distance = "bray", k=2, parallel = 8, try =20, trymax = 25)
#ordiplot(nmds_sub)
#stressplot(nmds_sub)
pl <- ordiplot(nmds_sub, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = ordi_extras$ID_plot,
             label = T, col = "red")
#plots
nmds_pl <- metaMDS(ordispe_pl, distance = "bray", k=2, trymax = 100)
ordiplot(nmds_pl)
ordiellipse (nmds_pl, groups = ordispe_pl_extras$ID_plot)
ordihull (nmds_pl, groups = ordispe_pl_extras$ID_plot, label = T)
stressplot(nmds_pl)
pl <- ordiplot(nmds_pl, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = ordispe_pl_extras$ID_plot,
             label = T, col = "red")

NMDS_sub_scores <- vegan::scores(nmds_sub)
NMDS_pl_scores <- vegan::scores(nmds_pl)

NMDS_sub_sites_scores <- NMDS_sub_scores$sites
NMDS_sub_sites_scores <- NMDS_sub_sites_scores %>% as.data.frame() %>% 
  rownames_to_column(var="ID_fine2")
NMDS_sub_sites_scores <- left_join(NMDS_sub_sites_scores, extras, by = "ID_fine2")

NMDS_sub_sp_scores <- NMDS_sub_scores$species

NMDS_pl_sites_scores <- NMDS_pl_scores$sites
NMDS_pl_sites_scores <- NMDS_pl_sites_scores %>% as.data.frame() %>% 
  rownames_to_column(var="ID_fine") %>% distinct()
# NMDS_pl_sites_scores <- left_join(NMDS_pl_sites_scores, select(extras, -c(ID_fine2, ID_subplot)), by = "ID_fine") %>% 
#   distinct() %>%
#   filter(!ID_site %in% c("AT01", "CZ02", "IT01","IT02",
#                          "IT03","IT04","IT05","IT06",
#                          "IT07","IT08","IT09","IT10",
#                          "IT11","IT12", "LT03", "ES02",
#                          "DE02", "RU04", "RU13"))

NMDS_pl_sp_scores <- NMDS_pl_scores$species

ggplot(NMDS_sub_sites_scores, aes(NMDS1, NMDS2, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow()) +
  theme_light()

df_labels <- aggregate(cbind(NMDS1, NMDS2) ~ ID_plot, data = NMDS_pl_sites_scores, FUN = mean)

# Create the plot
ggplot(filter(NMDS_sub_sites_scores, NMDS1 > -10),aes(x = NMDS1, y = NMDS2, group = ID_plot, color = ID_plot)) +
  #geom_path(linewidth=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Connect points within each group
  geom_point() + # Plot points
  geom_text(data = df_labels, aes(label = ID_plot, x = NMDS1, y = NMDS2), vjust = 0,
            check_overlap = T) + # Add group labels
  theme_bw() +  # Optional: Use a minimal theme for the plot
  theme(legend.position = "none")
  #labs(color = "Group") # Label the legend

# Create the plot
ggplot(NMDS_pl_sites_scores, aes(x = NMDS1, y = NMDS2, group = ID_plot, color = ID_plot)) +
  geom_path(linewidth=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Connect points within each group
  geom_point() + # Plot points
  geom_text(data = df_labels, aes(label = ID_plot, x = NMDS1, y = NMDS2), vjust = -1,
            check_overlap = T) + # Add group labels
  theme_bw() +  # Optional: Use a minimal theme for the plot
  theme(legend.position = "none")
#labs(color = "Group") # Label the legend

#Highlight

# Add a new column for highlight color
highlight_group <- c("SE14_0001") # Specify the group you want to highlight
NMDS_pl_sites_scores$colour <- ifelse(NMDS_pl_sites_scores$ID_plot == highlight_group, "red", "grey") # Change "red" to your preferred highlight color
NMDS_sub_sites_scores$colour <- ifelse(NMDS_sub_sites_scores$ID_plot == highlight_group, "red", "grey") # Change "red" to your preferred highlight color

# Calculate representative points for labels, considering only the highlight color for the label color
df_labels_hi <- aggregate(cbind(NMDS1, NMDS2) ~ ID_plot + colour, data = NMDS_pl_sites_scores, FUN = mean)
df_labels_sub_hi <- aggregate(cbind(NMDS1, NMDS2) ~ ID_plot + colour, data = NMDS_sub_sites_scores, FUN = mean)

ggplot(NMDS_pl_sites_scores, aes(x = NMDS1, y = NMDS2, group = ID_plot)) +
  geom_path(aes(colour = colour),size=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Use the new color column for line colors
  geom_point(aes(colour = colour)) + # Use the new color column for point colors
  scale_color_identity() + # Tell ggplot2 to use the colors as is
  geom_text(data = df_labels_hi, aes(label = ID_plot, x = NMDS1, y = NMDS2, colour = colour), vjust = -1, show.legend = FALSE) + # Add group labels
  theme_bw() +  # Optional: Use a minimal theme for the plot
  labs(color = "Group") # Adjust if you want a legend title, or use guides(color=FALSE) to remove the legend

ggplot(filter(NMDS_sub_sites_scores, NMDS1 > -10), aes(x = NMDS1, y = NMDS2, group = ID_plot)) +
  #geom_path(aes(colour = colour),size=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Use the new color column for line colors
  geom_point(aes(colour = colour)) + # Use the new color column for point colors
  scale_color_identity() + # Tell ggplot2 to use the colors as is
  geom_text(data = df_labels_sub_hi, aes(label = ID_plot, x = NMDS1, y = NMDS2, colour = colour), vjust = -1, show.legend = FALSE) + # Add group labels
  theme_bw() +  # Optional: Use a minimal theme for the plot
  labs(color = "Group") # Adjust if you want a legend title, or use guides(color=FALSE) to remove the legend

#PCA####
PCA_g <- rda(ordispe_sub)
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

#barplot(as.vector(PCA_gh$CA$eig)/sum(PCA_gh$CA$eig)) 
#sum((as.vector(PCA_gh$CA$eig)/sum(PCA_gh$CA$eig))[1:3])

table(VG_g_wide$ID_subplot)
table(VG_g_wide$ID_plot)

#CCA####
test.cca <- cca(ordispe_sub ~ ordienv)
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
ordispe_h2 <- decostand(ordispe_pl, "hellinger")

#DCA
DCA_g2 <- decorana(ordispe_pl, iweigh = 0)
ordiplot(DCA_g2)
pl <- ordiplot(DCA_g2, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = ordispe_pl_extras$ID_plot,
             label = T, col = "red")
#PCA
#on hell transformed data
PCA_gh2 <- rda(ordispe_h2)
pl <- ordiplot(PCA_gh2, display = "sites")
ordispider (pl, groups = ordispe_pl_extras$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh2, display = "species")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = ordispe_pl_extras$ID_plot,
             label = T, col = "red")
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
dist <- vegdist(ordispe_sub,  method = "bray")
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
test_part2 <- as.data.frame(test_part) %>% select( -c(survey_year,first, last,ID,ID_site,ID_plot,ID_subplot,ID_fine, ID_fine2,
                                                         TEMP,NH4M,NO3M,SO4SM,PREC,TEMP,latitude,longitude)) %>%
  sign() %>% 
  cbind(test_part$ID_fine) %>% 
  rename(ID_fine = 'test_part$ID_fine') %>% select(ID_fine, everything())
#make a list of matrices, one for each plot/year combo
test_list <- test_part2 %>% 
  split(.$ID_fine) %>% lapply(., function(x) { x["ID_fine"] <- NULL; x }) %>% 
  lapply(.,function(x) {as.matrix(x)})


# generate species sums these will be NA is any are missing
csum <- colSums(test_part2[,-c(1)])
# check if any are missing
any(is.na(csum))
# yes, some missing, so which ones?
which(is.na(csum))

#apply betapart core to create a betapart object for each matrix
core_list <- lapply(test_list,function(x) betapart.core(x))

multi_list <- lapply(core_list,function(x) beta.multi(x))

multi_df <- data.frame(matrix(unlist(multi_list), nrow=88, byrow=TRUE),stringsAsFactors=FALSE) %>% 
  cbind(names(multi_list))

multi_df <- multi_df %>% rename(beta.SIM = X1, beta.SNE = X2, beta.SOR = X3, ID_fine = 'names(multi_list)')

#LT plots report only one value for each year, not possible to calculate beta 
#diversity, and betapart returns NaN

all_data <- left_join(all_data, multi_df)


## Matrix temperature
# out <- nestedtemp(ordispe_sub)
# out
# plot(out)
# plot(out, kind="incid")
## Use oecosimu to assess the non-randomness of checker board units
#nestedchecker(ordispe_sub)
#oecosimu(ordispe_sub, nestedchecker, "quasiswap")
## Another Null model and standardized checkerboard score
#oecosimu(ordispe_sub, nestedchecker, "r00", statistic = "C.score")

#library(adespatial)

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
#combosites6 <- drop_na(combosites6, pc_dist)

#Plot level final data #####
plot.level.dat <- combosites6
plot.level.dat$year.i <- I(plot.level.dat$survey_year -1997)


#subplotlevel PCA movements
#Add betapart data
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

#Figure
ggplot(plot.level.dat %>% drop_na(pc_dist_base) %>%
         drop_na(ID_plot), aes(survey_year, pc_dist_base, colour = ID_plot)) +
  #  geom_point()+
  #  geom_smooth(method = "lm")+
  geom_line(size=1) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
  )

#MODELS with env variables####

#X Dispersion, betapart reponses X####
#simple multiple regression####
ml <- lm(dispersion  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           PREC + survey_year, data = subplot.level.dat)
summary(ml)

ml <- lm(beta.SOR  ~  NH4M * NO3M + SO4SM  + latitude + longitude +
           PREC + survey_year, data = subplot.level.dat)
summary(ml)

ml <- lm(beta.SNE  ~  NTOT + SO4SM  + latitude + longitude +
           PREC + survey_year, data = plot.level.dat)

summary(lm)
#seems more going on with PC2 as response
library(visreg)
visreg(ml, "NO3M")
visreg(ml, "NH4M")
visreg(ml, "NTOT")
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

summary(mod<-lme(pc_dist~NTOT + SO4SM + latitude + longitude,
                 random=~1|ID_fine,
                 correlation=corCAR1(form=~year.i|ID_fine),
                 data=plot.level.dat))

summary(mod<-lme(beta.SOR~NH4M + NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_fine,
                 #correlation=corCAR1(form=~year.i|ID_fine),
                 data=drop_na(subplot.level.dat)))

summary(mod<-lme(ID_fine_disp~year.i,
                 random=~1|ID_fine,
                 correlation=corCAR1(form=~year.i|ID_fine),
                 data=subplot.level.dat))
anova(mod)
fixed.effects(mod)

plot(mod, col = as.numeric(factor(subplot.level.dat$year.i,
                                  levels = unique(subplot.level.dat$year.i))),
     pch = 16, main = "resids")


#X PCA distance to base reponses X####
#combine PCA base distance and env variables and model
library(corrplot)
combosites5 %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()
combosites5 %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() 
subplot.level.dat %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()

#simple multiple regression####
ml <- lm(pc_dist  ~  NTOT+ SO4SM  + latitude + longitude +
           #beta.SNE +
           #beta.SIM +
           #beta.SOR +
           PREC +
           year.i, data = plot.level.dat)
summary(ml)
library(car)
vif(ml)
#attributes(alias(ml)$Complete)$dimnames[[1]] #beta.SOR is linearly dependant on SOR + SNE

ml <- lm(pc_dist_base  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           #beta.SOR +
           PREC + year.i, data = subplot.level.dat)
summary(ml)

ml <- lm(pc_dist_base  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           beta.SOR+
           PREC + year.i, data = plot.level.dat)
summary(ml)


glm <- glm(pc_dist_base  ~  NH4M + NO3M + SO4SM + TEMP + latitude + longitude +
             PREC + year.i, family = gaussian, data = plot.level.dat)

summary(glm)
#seems more going on with PC2 as response
library(visreg)
visreg(ml, "NO3M")
visreg(ml, "NH4M")
visreg(ml, "SO4SM")
visreg(ml, "NO3M", by = "year.i")
visreg(ml, "NH4M", by = "year.i")
visreg(ml, "NO3M", by = "year_scaled")
visreg(ml, "NH4M", by = "year_scaled")
visreg(ml, "NH4M", by = "NO3M")
visreg(ml, "NO3M", by = "latitude")
visreg(ml, "SO4SM", by = "year.i")
visreg(ml, "SO4SM", by = "TEMP")


#Nmle AR regressions####

#index year?
library(nlme)
library(lme4)
plot.level.dat$year.i <- I(plot.level.dat$survey_year-1997)
plot.level.dat$year.i <- as.numeric(plot.level.dat$year.i)

barplot(table(plot.level.dat$ID_plot))
barplot(table(plot.level.dat$latitude))
barplot(plot.level.dat$NTOT)
barplot(plot.level.dat$year.i)
barplot(table(plot.level.dat$year.i))
hist(plot.level.dat$pc_dist_base)#right skewed
hist(log(plot.level.dat$pc_dist_base))#log transf?

summary(nmod1 <- lme(pc_dist~
                       NH4M +
                       NO3M +
                       SO4SM +
                       latitude + longitude,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=plot.level.dat))

summary(nmod2 <- lme(pc_dist_base~
                       #NH4M * NO3M +
                       NTOT + 
                       SO4SM +
                       latitude + longitude +
                       year.i,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=plot.level.dat))

summary(nmod3 <- lme(log(pc_dist_base+1)~
                       #NH4M * NO3M +
                       NTOT + 
                       SO4SM +
                       latitude + longitude +
                       year.i,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=plot.level.dat))
AIC(nmod1,nmod2,nmod3)

#basic LM
summary(mod<-lm(pc_dist_base~
                  #NH4M * NO3M +
                  NTOT + 
                  SO4SM +
                  latitude + longitude +
                  year.i,
                data=plot.level.dat))

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
             year.i,
           data=subplot.level.dat,
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
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1.5 <- lmer(log(pc_dist_base)~
                 #NH4M * NO3M +
                 NTOT + 
                 SO4SM +
                 latitude +
                 longitude +
                 year.i+
                 (1|ID_subplot),
               data=subplot.level.dat)

lmm2 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm3 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               #year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm4 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i +
               (1|ID_fine),
             
             data=subplot.level.dat)

AIC(GLM, lmm1, lmm1.5, lmm2,lmm3,lmm4)
plot(lmm1.5)
plot(ranef(lmm1.5)) 
# QQ plots (drawn to the same scale!)
par(mfrow = c(1,2))
lims <- c(-3.5,3.5)
qqnorm(resid(GLM, type = "pearson"),
       xlim = lims, ylim = lims,main = "GLM")
abline(0,1, col = "red", lty = 2)
qqnorm(resid(lmm1.5, type = "pearson"),
       xlim = lims, ylim = lims, main = "lmm1")
abline(0,1, col = "red", lty = 2)

mod <- gam1

summary(lmm1)
#percentage variance explained by random factor ID_plot
0.0007465/(0.0007465 + 0.0026378)#22%

summary(lmm1.5)
#percentage variance explained by random factor ID_plot
0.0009634/(0.0009634 + 0.0035281)#22%

summary(lmm2)
#percentage variance explained by random factor ID
0.002088/(0.002088 + 0.001063)#66%

anova(lmm1, lmm1.5)

qqnorm(resid(mod))#2016_SE14_1 is an outlier on 0.20+
qqline(resid(mod))
#sjPlot
plot_model(lmm1.5, type = "re", show.values = TRUE)
plot_model(lmm1.5, show.values = TRUE)


#Facet plots of model predictions####
#comboplots
(mm_plot <- ggplot(plot.level.dat, aes(x = year.i, y = pc_dist_base, colour = ID_plot)) +
   facet_wrap(~ID_plot, nrow=2) +   # a panel for each mountain range
   geom_point(alpha = 0.5) +
   theme_classic() +
   geom_line(data = cbind(plot.level.dat, pred = predict(mod)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
   theme(legend.position = "none",
         panel.spacing = unit(2, "lines"))  # adding space between panels
)

#subplot.level.dat
(mm_plot <- ggplot(subplot.level.dat, aes(x = year.i, y = pc_dist_base, colour = ID_plot)) +
    facet_wrap(~ID_plot, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(subplot.level.dat, pred = predict(mod)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)


summary(mod<-lme(pc_dist_base~NH4M * NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_subplot,
                 correlation=corCAR1(form=~year.i|ID_subplot),
                 data=drop_na(subplot.level.dat)))


#GAM####
#try same model with bam
k=10
bs="ts"
gam0 <- bam(round(pc_dist*100) ~ s(NH4M) +
              s(NO3M)+
              #s(SO4SM.lag,k=k, bs=bs)+
              s(NTOT.lag,k=k, bs=bs)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              s(year.i, bs="re")+
              #s(year.i,k=k, bs=bs)+
              #s(ID_site, bs= "re")+
              s(ID_plot, bs= "re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = poisson,
            data = subplot.level.dat)
#data = subplot.level.dat)

gam0f <- bam(round(pc_dist*100) ~ #s(NH4M) +
               #s(NO3M)+
               s(SO4SM.lag,k=k, bs=bs)+
               s(NTOT.lag,k=k, bs=bs)+
               s(PREC)+
               #s(latitude,k=k, bs=bs)+
               #s(longitude,k=k, bs=bs)+
               te(latitude, longitude,k=k, bs=bs)+
               s(year.i,k=k, bs=bs)+
               s(ID_subplot, bs="re"),
             nthreads=8, 
             family = poisson,
             data = subplot.level.dat)
#data = subplot.level.dat)

gam1 <- bam(round(pc_dist*100) ~ #s(NH4M) +
              #s(NO3M)+
              s(SO4SM.lag,k=k, bs=bs)+
              s(NTOT.lag,k=k, bs=bs)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              s(year.i,k=k, bs=bs)+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = poisson,
            data = subplot.level.dat)
#data = subplot.level.dat)

gam1.5 <- bam(round(pc_dist*100) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT.lag,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(year.i,k=k, bs=bs)+
                s(ID_plot, bs="re")+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = poisson,
              data = subplot.level.dat)
#data = subplot.level.dat)

gam2 <- bam(log(pc_dist) ~ #s(NH4M) +
              #s(NO3M)+
              #s(SO4SM,k=k, bs=bs)+
              s(NTOT,k=k, bs=bs)+
              #s(beta.SOR)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              #year.i+
              s(year.i, bs=bs)+
              s(ID_plot,bs="re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = gaussian,
            method = "REML",
            #data = all_data_dca)
            data = subplot.level.dat.no)
#data = plot.level.dat)

gam2.1 <- bam(round(pc_dist*100) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT.lag,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(year.i,k=k, bs=bs)+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = poisson,
              data = subplot.level.dat)


gam2.2 <- bam(sqrt(pc_dist) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT.lag,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(year.i, bs=bs)+
                s(ID_plot,bs="re")+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = betar,
              data = subplot.level.dat)

gam2.3 <- bam(log(pc_dist) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(year.i,k=k, bs=bs)+
                s(ID_plot, bs="re")+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = gaussian,
              data = subplot.level.dat)
#data = plot.level.dat)

gam3 <- bam(log(pc_dist_base) ~ s(NH4M) +
              s(NO3M)+
              #s(NTOT)+
              s(SO4SM)+
              ti(SO4SM,NO3M,NH4M)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              s(year.i,k=15, bs=bs)+
              s(ID_plot, bs= "re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            select = TRUE,
            family = gaussian,
            data = subplot.level.dat)

gam4 <- bam(log(pc_dist) ~ #s(NH4M) +
              #s(NO3M)+
              #s(SO4SM,k=k, bs=bs)+
              #s(NTOT,k=k, bs=bs)+
              te(NTOT, SO4SM)+
              #s(beta.SOR)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              #year.i+
              s(year.i, bs=bs)+
              s(ID_plot,bs="re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = gaussian,
            method = "REML",
            data = subplot.level.dat.no)

# qdat <- select(subplot.level.dat, NTOT, SO4SM, year.i, ID_plot, ID_subplot, latitude, longitude,
#                pc_dist) %>% drop_na()

b <- getViz(gam3)
summary(b)
draw(b, parametric = FALSE, residuals = TRUE)
appraise(b)
acf(residuals(b))
concurvity(b, full = FALSE)
library(sjPlot)
tab_model(b, show.ci = FALSE, show.stat = TRUE, show.r2 = FALSE,
          show.obs = FALSE, string.stat = "Chi.sq", string.est = "edf", file ="GAM1.html")#on dat 1493.023
library(dsm)#vis.concurvity

vis_concurvity(b, type = "estimate")
vis_concurvity(b, type = "worst")

AIC(gam0, gam0f, gam1,gam1.5,gam2.2,gam2.3, gam3)

#X PCA sequential distance reponses X####

#lm
summary(lm1 <- lm(log(pc_dist)~
                    NTOT.lag +
                    #NH4M +
                    #NO3M +
                    SO4SM.lag +
                    #TEMP +
                    PREC +
                    #SO4SM.lag +
                    latitude +
                    longitude +
                    year.i, 
                  data=subplot.level.dat))
# data=plot.level.dat))
gvlma(lm1)
acf(resid(lm1))
mean(lm1$residuals)
#nmle
hist(subplot.level.dat$pc_dist)
hist(log(subplot.level.dat$pc_dist))#log transf?


summary(nmod1 <- lme(pc_dist~
                       NTOT.lag +
                       year.i ,
                     #NH4M +
                     #NO3M +
                     #SO4SM.lag ,
                     #latitude + longitude,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=drop_na(select(plot.level.dat, - TEMP))))

summary(nmod2 <- lme(pc_dist~
                       NTOT.lag +
                       year.i + 
                       #NH4M +
                       #NO3M +
                       SO4SM.lag +
                       latitude + longitude,
                     random=~1|ID_plot/ID_subplot,
                     correlation=corCAR1(form=~year.i|ID_plot/ID_subplot),
                     data=drop_na(select(subplot.level.dat, - TEMP))))


summary(nmod3 <- lme(log10(pc_dist)~
                       NTOT.lag +
                       #year.i +
                       #NH4M * NO3M +
                       SO4SM.lag +
                       latitude + longitude,
                     random=~1|ID_plot/ID_subplot,
                     correlation=corCAR1(form=~year.i|ID_plot/ID_subplot),
                     data=drop_na(select(subplot.level.dat, - TEMP))))

summary(nmod4 <- lme(sqrt(pc_dist)~
                       NTOT.lag +
                       #year.i +
                       #NH4M * NO3M +
                       SO4SM.lag +
                       latitude + longitude,
                     random=~1|ID_plot/ID_subplot,
                     correlation=corCAR1(form=~year.i|ID_plot/ID_subplot),
                     data=drop_na(select(subplot.level.dat, - TEMP))))

AIC(nmod1,nmod2,nmod3,nmod4)
plot(nmod1)
anova(nmod1)
fixed.effects(nmod1)
acf(resid(nmod1))
nmod1
plot(nmod1, col = c(1:nlevels(as.factor(subplot.level.dat$ID_plot))), pch = 16)
vif(nmod1)

#Lmer####
#GLM null model
GLM <- gls(log(pc_dist)~
             #NH4M * NO3M +
             NTOT + 
             SO4SM +
             latitude + longitude +
             year.i,
           data=subplot.level.dat,
           # family = "poisson",
           method = "ML")
summary(GLM)
plot(GLM)

#GLM null model
GLM2 <- glm(pc_dist~
              #NH4M * NO3M +
              NTOT + 
              SO4SM +
              latitude +
              longitude +
              year.i,
            data=subplot.level.dat,
            family = "poisson")
summary(GLM2)
plot(GLM2)
plot_model(GLM2, show.values = TRUE) #sjplot

#LME mixed
lmm0 <- lmer(pc_dist~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1.5 <- lmer(log(pc_dist)~
                 #NH4M * NO3M +
                 NTOT + 
                 #SO4SM +
                 latitude +
                 longitude +
                 year.i+
                 (1|ID_subplot),
               data=subplot.level.dat)

lmm2 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT.lag + 
               #SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm3 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               #year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm4 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               (1|year.i) +
               (1|ID_fine),
             
             data=subplot.level.dat)

lmm5 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               #SO4SM +
               latitude +
               longitude +
               (1|ID_plot),
             #(1|year.i),
             #(1|ID_fine),
             
             data=subplot.level.dat)

AIC(GLM, lmm1, lmm1.5, lmm2,lmm3,lmm4,lmm5)

mod <- lmm5

plot(mod)
plot(ranef(mod)) 
# QQ plots (drawn to the same scale!)
par(mfrow = c(1,2))
lims <- c(-3.5,3.5)
qqnorm(resid(GLM, type = "pearson"),
       xlim = lims, ylim = lims,main = "GLM")
abline(0,1, col = "red", lty = 2)
qqnorm(resid(mod, type = "pearson"),
       xlim = lims, ylim = lims, main = "lmm1")
abline(0,1, col = "red", lty = 2)
summary(mod)
#percentage variance explained by random factor ID_plot
0.0007465/(0.0007465 + 0.0026378)#22%

summary(lmm1.5)
#percentage variance explained by random factor ID_plot
0.0009634/(0.0009634 + 0.0035281)#22%

summary(lmm2)
#percentage variance explained by random factor ID
0.002088/(0.002088 + 0.001063)#66%

anova(lmm1, lmm1.5)

qqnorm(resid(mod))#2016_SE14_1 is an outlier on 0.20+
qqline(resid(mod))
#sjPlot
plot_model(mod, type = "re", show.values = TRUE)
plot_model(mod, show.values = TRUE)

acf(resid(mod))

library(ggeffects)
ggpredict(mod, terms = c("NTOT")) %>% plot()
ggpredict(mod, terms = c("NTOT", "ID_plot"), type = "re") %>% plot() +
  theme(legend.position = "bottom")

#MCMCglmm####
library(MCMCglmm)
library(MCMCvis)
mod <- MCMCglmm(round(pc_dist*100) ~ NH4M.lag + year.i, random = ~ID_subplot,
                family = "poisson", data = subplot.level.dat)

summary(mod)
plot(mod$Sol)
plot(mod$VCV)
MCMCplot(mod$Sol)
MCMCplot(mod$VCV)


#glmmTMB####
library(glmmTMB)
library(broom)
library(DHARMa)
mod0 <- glmmTMB(log(pc_dist) ~ NTOT + (1|year.i) + #(1|year.i/ID_subplot) +
                  (1|ID_site/ID_plot/ID_subplot),
                #(1|ID_fine),
                family = "gaussian", data = subplot.level.dat)

mod0 <- glmmTMB(log(pc_dist) ~ NTOT + SO4SM +(1|year.i) + #(1|year.i/ID_subplot) +
                  (1|ID_site)+ (1|ID_plot)+ (1|ID_subplot),
                #(1|ID_fine),
                family = "gaussian", data = subplot.level.dat)

mod1 <- glmmTMB(round(pc_dist*100)~ NTOT +
                  #SO4SM +
                  year.i+
                  (1|year.i) + #(1|year.i/ID_subplot) +
                  (1|ID_site)+
                  (1|ID_plot)+
                  (1|ID_subplot),
                #(1|ID_fine),
                family = "poisson", data = subplot.level.dat)
mod <- mod1
summary(mod)
res = simulateResiduals(mod)
plot(res, rank = T)
acf(resid(mod))

mod1 <- glmmTMB(round(pc_dist*100) ~ NTOT + year.i +(1|ID_subplot) +
                  (1|ID_plot),
                family = "poisson", data = subplot.level.dat)


mod2 <- glmmTMB(log(pc_dist) ~ NTOT + year.i +
                  #beta.SIM +
                  (1|ID_subplot) +
                  (1|ID_plot),
                family = "gaussian", data = subplot.level.dat.no)

mod3 <- glmmTMB(log(pc_dist) ~ NTOT +(1|year.i) +(1|ID_subplot) +
                  (1|ID_plot),
                family = "gaussian", data = subplot.level.dat.no)

mod4 <- glmmTMB(sqrt(pc_dist) ~ NTOT + (1|year.i) +(1|ID_subplot) +
                  (1|ID_plot),
                family = "gaussian", data = subplot.level.dat)

mod5 <- glmmTMB(pc_dist ~ NTOT + year.i +(1|ID_subplot) +
                  (1|ID_plot), beta_family(),
                data = subplot.level.dat)

mod6 <- glmmTMB(pc_dist ~ NTOT +(1|ID_subplot/year.i) +
                  (1|ID_plot/year.i), beta_family(),
                data = subplot.level.dat)

AIC(mod1, mod2,mod3,mod4,mod5,mod6)

BIC(mod1, mod2,mod3,mod4,mod5,mod6,GLM, lmm1, lmm1.5,
    lmm2,lmm3,lmm4,lmm5,gam0, gam0f, gam1,gam1.5, gam2,gam2.2,
    gam2.3, gam3) %>% arrange(BIC)

mod <- mod2
summary(mod)
res = simulateResiduals(mod)
#res = simulateResiduals(mod, form = subplot.level.dat$year.i)
plot(res, rank = T)

out <- boxplot.stats(sqrt(subplot.level.dat$pc_dist))$out
out_ind <- which(sqrt(subplot.level.dat$pc_dist) %in% c(out))
out_ind

out <- boxplot.stats(log(subplot.level.dat$pc_dist))$out
out_ind <- which(log(subplot.level.dat$pc_dist) %in% c(out))
out_ind
outliers <- subplot.level.dat[out_ind, ]

subplot.level.dat.no <- subplot.level.dat[-out_ind, ]

acf(resid(mod))

testDispersion(res)
testUniformity(res)
testOutliers(res)
testZeroInflation(res)
testTemporalAutocorrelation(res.t)
testSpatialAutocorrelation(res, x= subplot.level.dat$longitude, y = subplot.level.dat$latitude)

res.t <- recalculateResiduals(res, group = subplot.level.dat$year.i)
plot(res, quantreg = FALSE)



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


