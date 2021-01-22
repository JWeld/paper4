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

hell_VG_fm <- decostand(VG_fm, "hellinger") 
hell_VG_bm <- decostand(VG_bm, "hellinger") 
hell_VG_gm <- decostand(VG_gm, "hellinger") 

DCA_f <- decorana(VG_fm, iweigh = 0)
DCA_b <- decorana(VG_bm, iweigh = 0)
DCA_g <- decorana(VG_gm, iweigh = 0)

pl <- ordiplot(DCA_f)
ordispider (pl, groups = VG_f_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_b)
ordispider (pl, groups = VG_b_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

pl <- ordiplot(DCA_g)
ordispider (pl, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

#PCA
PCA_g <- rda(VG_gm)
pl <- ordiplot(PCA_g, display = "sites")
ordispider (pl, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
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

#text(pl, "sites", col="blue", cex=0.9)

#CA
CA_g <- cca(VG_gm)
pl <- ordiplot(CA_g, display = "sites")


VS_h <- decostand(VS_m, "hellinger")
PCA_VS <- rda(VS_h)
pl <- ordiplot(PCA_VS, display = "sites")
ordispider (pl, groups = VS_under_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

#NMDS
#nmds_f <- metaMDS(VG_fm, k=2, try = 20, parallel = getOption("mc.cores"))
#nmds_b <- metaMDS(VG_bm, k=2, try = 20, parallel = getOption("mc.cores")) #no convergence
#nmds_g <- metaMDS(VG_gm, k=2, try = 20, parallel = getOption("mc.cores"))

ordiplot(nmds_g, type="n")
ordihull(nmds_g, groups=VG_g_wide$ID_site ,draw="polygon",col="grey90",label=F)
orditorp(nmds_g, display="sites",col="red",air=0.01)
ordihull(nmds_g, groups=VG_g_wide$ID_site ,draw="polygon",col="grey90",label=F)
ordispider (nmds_f, groups = VG_g_wide$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)

