library(tidyverse)
library(vegan)
library(parallel)

#data we're using
#VG_f_wide - field layer
#VG_b_wide - bottom layer
#VG_g_wide -combo of above, ground layer

# VG_fm <- select(VG_f_wide, -c(1:10))
# VG_bm <- select(VG_b_wide, -c(1:10))
# VG_gm <- select(VG_g_wide, -c(1:10))

#And matrix versions of the above
# VG_fm <- select(VG_f_wide, -c(1:10))
# VG_bm <- select(VG_b_wide, -c(1:10))
# VG_gm <- select(VG_g_wide, -c(1:10))



#Circular plots (VS)####

#data we're using
#VG_f_wide - field layer
#VG_b_wide - bottom layer
#VG_g_wide -combo of above, ground layer

#PCA
hell_VS_fm <- decostand(VS_fm, "hellinger") 
hell_VS_bm <- decostand(VS_bm, "hellinger") 
hell_VS_gm <- decostand(VS_gm, "hellinger") 

#DCA
DCA_f <- decorana(VS_fm, iweigh = 1)
DCA_b <- decorana(VS_bm, iweigh = 1)
DCA_g <- decorana(VS_gm, iweigh = 1)

pl <- ordiplot(DCA_f, display = "sites")
ordispider (pl, groups = VS_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_b)
ordispider (pl, groups = VS_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_g)
ordispider (pl, groups = VS_g_wide$ID, col = c("red","blue","green","grey","pink","yellow"), label = T)
ordihull (pl, groups = VS_g_wide$ID_site, col = c("red","blue","green","grey70","pink","yellow","grey30","purple"),
          draw= "polygon", alpha = 100, 
          label = T)
ordiellipse (pl, groups = VS_g_wide$ID_site, col = c("red","blue","green","grey70","pink","yellow","grey30","purple"),
             draw= "polygon", alpha = 100, 
             label = T)
ordiellipse (pl, groups = VS_g_wide$ID, col = c("red","blue","green","grey70","pink","yellow","grey30","purple"),
             draw= "polygon", alpha = 100, 
             label = F)
#PCA
PCA_g <- rda(VS_gm)
pl <- ordiplot(PCA_g, display = "sites")
ordispider (pl, groups = VS_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
#on hell transformed data
PCA_fh <- rda(hell_VS_fm)
PCA_bh <- rda(hell_VS_bm)
PCA_gh <- rda(hell_VS_gm)
#pl <- ordiplot(PCA_gh, display = "sites")
#ordihull(pl, VG_g_wide$AreaCode, label = TRUE) 
pl <- ordiplot(PCA_fh, display = "sites")
ordispider (pl, groups = VS_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_bh, display = "sites")
ordispider (pl, groups = VS_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "sites")
ordispider (pl, groups = VS_gr$ID, col = c("red","blue","green","grey","pink","yellow"), label = T)
ordiellipse (pl, groups = VS_g_wide$ID, col = c("red","blue","green","grey70","pink","yellow","grey30","purple"),
             draw= "polygon", alpha = 0.2, 
             label = F)

#NMDS
#nmds_f <- metaMDS(VG_fm, k=2, try = 20, parallel = getOption("mc.cores"))
#nmds_b <- metaMDS(VG_bm, k=2, try = 20, parallel = getOption("mc.cores")) #no convergence
#nmds_g <- metaMDS(VG_gm, k=2, try = 20, parallel = getOption("mc.cores"))

ordiplot(nmds_g, type="n")
ordihull(nmds_g, groups=VG_g_wide$ID_site ,draw="polygon",col="grey90",label=F)
orditorp(nmds_g, display="sites",col="red",air=0.01)
ordihull(nmds_g, groups=VG_g_wide$ID_site ,draw="polygon",col="grey90",label=F)
ordispider (nmds_f, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)












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
ordiellipse (pl, groups = VG_g_wide$ID,
             label = F)
#text(pl, "sites", col="blue", cex=0.9)


#Intensive plots per intensive/year####

hell_VG_gm2 <- decostand(VG_gm2, "hellinger") 

#DCA
DCA_g2 <- decorana(VG_gm2, iweigh = 1)
pl <- ordiplot(DCA_g2)
ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
ordiellipse (pl, groups = VG_g_wide2$ID,
             label = F)
#PCA
PCA_g <- rda(VG_gm2)
pl <- ordiplot(PCA_g, display = "sites")
ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
#on hell transformed data
PCA_gh <- rda(hell_VG_gm2)
pl <- ordiplot(PCA_gh, display = "sites")
ordispider (pl, groups = VG_g_wide2$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "species")
ordiellipse (pl, groups = VG_g_wide2$ID,
             label = F)
text(pl, "sites", col="blue", cex=0.6)
#text(pl, "species", col="blue", cex=0.6)

#CA
CA_g <- cca(VG_gm2)
pl <- ordiplot(CA_g, display = "species")
text(pl, "species", col="blue", cex=0.5)

#ICP Forests data####
env_means <- readRDS("~/Documents/r/paper4/Forests_data/env_means2.RDS")
env_means_scaled <- readRDS("~/Documents/r/paper4/Forests_data/env_means_sc2.RDS")
#For_V <-  veg_vasc %>% group_by(ID, survey_number, species_name) %>% summarise(mean_cover = mean(cover))

hell_VF_vasc <- decostand(VF_vasc_m, "hellinger") 
hell_VF <- decostand(VF_m, "hellinger") 

#PCA
# PCA_VF_v <- rda(VF_vasc_m)
# pl <- ordiplot(PCA_VF_v, display = "sites")
# ordispider (pl, groups = VF_vasc$ID_site, col = c("red","blue","green","grey",
#                                                   "pink","yellow"), label = F)
# pl <- ordiplot(PCA_VF_v, display = "species")
# text(pl, "species", col="blue", cex=0.5)
# 
# 
# PCA_VF <- rda(VF_m)
# pl <- ordiplot(PCA_VF, display = "sites")
# ordispider (pl, groups = VF$ID_site, col = c("red","blue","green","grey",
#                                                   "pink","yellow"), label = F)
# 
# pl <- ordiplot(PCA_VF, display = "species")
# text(pl, "species", col="blue", cex=0.5)


#on hell transformed data
PCA_VF_v_h <- rda(hell_VF_vasc)
PCA_VF_h <- rda(hell_VF)
#pl <- ordiplot(PCA_gh, display = "sites")
#ordihull(pl, VG_g_wide$AreaCode, label = TRUE) 
#vasc only
pl <- ordiplot(PCA_VF_v_h, display = "sites")
ordispider (pl, groups = VF_vasc$ID_site, col = c("red","blue","green","grey",
                                                  "pink","yellow"), label = F)
#all
pl <- ordiplot(PCA_VF_h, display = "sites")
ordispider (pl, groups = VF$ID_site, col = c("red","blue","green","grey",
                                             "pink","yellow"), label = F)

ordiellipse (pl, groups = VF$ID_site,
             label = F)
#text(pl, "sites", col="blue", cex=0.9)

#extract scores
PCA_scores <- scores(PCA_VF_h)
PCA_scores_sites <- as.data.frame(PCA_scores$sites)
PCA_scores_species <- as.data.frame(PCA_scores$species)

PCA_scores_sites <- PCA_scores_sites %>% rownames_to_column(.,var = "ID")
extras <- VF %>% 
  select(ID, ID_site, country, survey_year)
PCA_scores_sites <- left_join(PCA_scores_sites, extras, by = "ID")
all_data_scaled <- left_join(PCA_scores_sites, env_means_scaled, by = "ID")
all_data <- left_join(PCA_scores_sites, env_means, by = "ID")

library(VIM)
matrixplot(all_data)
all_data %>% matrixplot()
all_data2 <- all_data %>% drop_na()
all_data2 <- rename(all_data2, survey_year = survey_year.x, year_scaled = survey_year.y)

all_data_scaled %>% matrixplot()
all_data_scaled_2 <- all_data_scaled %>% drop_na()
all_data_scaled_2 <- rename(all_data_scaled_2, survey_year = survey_year.x,
                            year_scaled = survey_year.y)

ggplot(all_data2, aes(x=n_no3, y = PC2, colour = country))+
  geom_point() + geom_smooth(method = "lm", se=FALSE)

#simple multiple regression####
ml <- lm(PC2  ~  n_nh4 * n_no3 + latitude + longitude + sum_canopy + mean_temp + 
           mean_precip + year_scaled, data = all_data_scaled_2)
summary(ml)

#seems more going on with PC2 as response
library(visreg)
visreg(ml, "n_no3", by = "survey_year")
visreg(ml, "n_nh4", by = "survey_year")
visreg(ml, "n_no3", by = "year_scaled")
visreg(ml, "n_nh4", by = "year_scaled")
visreg(ml, "n_nh4", by = "n_no3")
visreg(ml, "n_no3", by = "latitude")

#keep only sites with observations over at least three years
# n.obs.dat <-as.data.frame(table(all_data2$ID_site))
# drop.list <- filter(n.obs.dat, Freq <3)
# dat <- filter(all_data2, ID_site %!in% drop.list$Var1)
# dat2 <- dat

#save data
saveRDS(all_data2, file="dat.RDS")
saveRDS(all_data_scaled_2, file="dat_scaled.RDS")

