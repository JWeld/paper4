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
#Intensive plots (VG)####

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


