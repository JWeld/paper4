library(tidyverse)
library(vegan)

#Intensive plots (VG)####
hell_VG_fm <- decostand(VG_fm, "hellinger") 
hell_VG_bm <- decostand(VG_bm, "hellinger") 
hell_VG_gm <- decostand(VG_gm, "hellinger") 

#DCA
DCA_f <- decorana(VG_fm, iweigh = 0)
DCA_b <- decorana(VG_bm, iweigh = 0)
DCA_g <- decorana(VG_gm, iweigh = 0)

pl <- ordiplot(DCA_f)
ordispider (pl, groups = VG_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_b)
ordispider (pl, groups = VG_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_g)
ordispider (pl, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
ordiellipse (pl, groups = VG_g_wide$ID,
             label = F)
#PCA
PCA_g <- rda(VG_gm)
pl <- ordiplot(PCA_g, display = "sites")
#text(pl, "sites", col="blue", cex=0.5)
ordispider (pl, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(PCA_g, display = "species")
text(pl, "species", col="blue", cex=0.5)


#on hell transformed data
PCA_fh <- rda(hell_VG_fm)
PCA_bh <- rda(hell_VG_bm)
PCA_gh <- rda(hell_VG_gm)
#pl <- ordiplot(PCA_gh, display = "sites")
#ordihull(pl, VG_g_wide$AreaCode, label = TRUE) 
pl <- ordiplot(PCA_fh, display = "sites")
ordispider (pl, groups = VG_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_bh, display = "sites")
ordispider (pl, groups = VG_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "sites")
ordispider (pl, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "sites")
ordiellipse (pl, groups = VG_g_wide$ID_subplot,
             label = F)
#text(pl, "sites", col="blue", cex=0.9)

table(VG_g_wide$ID_subplot)
#Intensive plots per intensive/year####

# hell_VG_gm2 <- decostand(VG_gm2, "hellinger") 
# 
# #DCA
# DCA_g2 <- decorana(VG_gm2, iweigh = 1)
# pl <- ordiplot(DCA_g2)
# ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
# ordiellipse (pl, groups = VG_g_wide2$ID,
#              label = F)
# #PCA
# PCA_g <- rda(VG_gm2)
# pl <- ordiplot(PCA_g, display = "sites")
# ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
# #on hell transformed data
# PCA_gh <- rda(hell_VG_gm2)
# pl <- ordiplot(PCA_gh, display = "sites")
# ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
# pl <- ordiplot(PCA_gh, display = "species")
# ordiellipse (pl, groups = VG_g_wide2$ID,
#              label = F)
# text(pl, "sites", col="blue", cex=0.6)
# #text(pl, "species", col="blue", cex=0.6)

#extract scores
PCA_scores <- scores(PCA_gh)
PCA_scores_sites <- as.data.frame(PCA_scores$sites)
PCA_scores_species <- as.data.frame(PCA_scores$species)

PCA_scores_sites <- PCA_scores_sites %>% rownames_to_column(.,var = "ID_fine2")
PCA_scores_species <- PCA_scores_species %>% rownames_to_column(.,var = "Species")

extras <- VG_g_wide %>% 
  select(ID, ID_site,ID_fine, ID_fine2, survey_year)
PCA_scores_sites <- left_join(PCA_scores_sites, extras, by = "ID_fine2")

#combine PCA scores with other data####
all_data <- left_join(PCA_scores_sites, IM_env, by = "ID")
#discard where deposition data is available
all_data <- filter(all_data, NH4M > 0)

#simple multiple regression####
ml <- lm(PC1  ~  NH4M * NO3M + SO4SM + TEMP + 
           PREC + survey_year, data = all_data)
summary(ml)

#seems more going on with PC2 as response
library(visreg)
visreg(ml, "NO3M", by = "survey_year")
visreg(ml, "NH4M", by = "survey_year")
visreg(ml, "NO3M", by = "year_scaled")
visreg(ml, "NH4M", by = "year_scaled")
visreg(ml, "NH4M", by = "NO3M")
visreg(ml, "NO3M", by = "latitude")
visreg(ml, "SO4SM", by = "survey_year")
visreg(ml, "SO4SM", by = "TEMP")


#keep only sites with observations over at least three years
# n.obs.dat <-as.data.frame(table(all_data2$ID_site))
# drop.list <- filter(n.obs.dat, Freq <3)
# dat <- filter(all_data2, ID_site %!in% drop.list$Var1)
# dat2 <- dat
