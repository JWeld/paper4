#run ordinations on VG only, no explanatory env variables, just exploratory

library(tidyverse)
library(vegan)
options(mc.cores = parallel::detectCores())
#library(gvlma)#checks assumptions for models

#Intensive plots (VG)####
#Data#
ordispe #matrix of subplots
ordispe_pl #matrix of plots
ordispe_pl_extras #factors columns removed from plots matrix
ordienv #factor columns removed from subplots matrix plus env variables
#hellinger transformed
ordispe_h <- decostand(ordispe, "hellinger") #hellinger transform to make usable with PCA
ordispe_h2 <- decostand(ordispe_pl, "hellinger")

#Ordinations####
#NMDS####
#subplots
nmds_sub <- metaMDS(ordispe, distance = "bray", k=2, parallel = 8, try =20, trymax = 25)
#ordiplot(nmds_sub)
stressplot(nmds_sub)
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
NMDS_pl_sites_scores <- left_join(NMDS_pl_sites_scores, select(extras, -c(ID_fine2, ID_subplot)), by = "ID_fine") %>% 
  distinct() %>% filter(!ID_site %in% c("AT01", "CZ02", "IT01","IT02","IT03","IT04","IT05","IT06",
                                        "IT07","IT08","IT09","IT10","IT11","IT12", "LT03", "ES02", "DE02"))

NMDS_pl_sp_scores <- NMDS_pl_scores$species

ggplot(NMDS_sub_sites_scores, aes(NMDS1, NMDS2, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow()) +
  theme_light()

df_labels <- aggregate(cbind(NMDS1, NMDS2) ~ ID_plot, data = NMDS_pl_sites_scores, FUN = mean)

# Create the plot
ggplot(NMDS_pl_sites_scores, aes(x = NMDS1, y = NMDS2, group = ID_plot, color = ID_plot)) +
  geom_path(size=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Connect points within each group
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

# Calculate representative points for labels, considering only the highlight color for the label color
df_labels_hi <- aggregate(cbind(NMDS1, NMDS2) ~ ID_plot + colour, data = NMDS_pl_sites_scores, FUN = mean)

ggplot(NMDS_pl_sites_scores, aes(x = NMDS1, y = NMDS2, group = ID_plot)) +
  geom_path(aes(colour = colour),size=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Use the new color column for line colors
  geom_point(aes(colour = colour)) + # Use the new color column for point colors
  scale_color_identity() + # Tell ggplot2 to use the colors as is
  geom_text(data = df_labels_hi, aes(label = ID_plot, x = NMDS1, y = NMDS2, colour = colour), vjust = -1, show.legend = FALSE) + # Add group labels
  theme_bw() +  
  labs(color = "Group") # Adjust legend title, or use guides(color=FALSE) to remove the legend

#PCA####
#on hell transformed data

#By subplots
PCA_gh <- rda(ordispe_h, scale = F)

summary(PCA_gh)
pl <- ordiplot(PCA_gh, display = "sites")
ordispider (pl, groups = test$ID_site, col = c("red","blue","green","grey","pink","yellow"), label = T)
pl <- ordiplot(PCA_gh, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = test$ID_plot,
             label = T, col = "red")

pl <- ordiplot(PCA_gh, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = test$ID_plot,
             label = T, col = "red")

pl <- ordiplot(PCA_gh, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = test$ID_site,
             label = T, col = "red")


# By plots ##
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


test.df <- PCA_scores_sites2 #CHANGE HERE AS NEEDED
test.df <- left_join(test.df, extras)

#Add PC movement from base and sequential distances to ####
#both subplot and plot level data
#Add distances between sequential observations and distance to base-line for
#all to allow plots like in Lamothe 2019

test.df <- test.df %>% ungroup() %>% group_by(ID_plot) %>% arrange(survey_year) %>% 
  mutate(pc1diff = PC1 - lag(PC1)) %>% mutate(pc2diff = PC2 - lag(PC2)) %>% 
  mutate(pc_dist = sqrt((pc1diff)^2 + (pc2diff)^2))

temp <- test.df %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  filter(row_number()==1) %>% select(ID_subplot, PC1, PC2) %>%
  rename(basePC1 = PC1, basePC2 = PC2)

test.df <- left_join(test.df, temp, by = "ID_subplot")
test.df <- test.df %>% ungroup() %>% group_by(ID_subplot) %>% arrange(survey_year) %>% 
  mutate(pc1diffbase = PC1 - basePC1) %>% mutate(pc2diffbase = PC2 - basePC2) %>% 
  mutate(pc_dist_base = sqrt((pc1diffbase)^2 + (pc2diffbase)^2))

ggplot(test.df %>% drop_na(pc_dist) %>% drop_na(ID_plot), aes(ID_plot, pc_dist, colour = ID_plot)) +
  geom_boxplot() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")


#PLOTS- cup and ball models####
ggplot(test.df, aes(PC1, PC2, colour = ID_plot)) + geom_point() +
  facet_wrap(facets = "ID_plot")

ggplot(test.df, aes(PC1, PC2, colour = ID_plot)) + geom_point() +
  facet_wrap(facets = "ID_site")

ggplot(test.df, aes(PC1, PC2,colour = ID_plot)) + geom_point()

#Figure####
ggplot(drop_na(test.df), aes(PC1, PC2, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow()) +
  theme_light()

ggplot(test.df %>% drop_na(pc_dist_base) %>%
         drop_na(ID_plot), aes(survey_year, pc_dist_base, colour = ID_plot)) +
  #geom_point()+
  geom_smooth(method = "loess", se = F)+
  #geom_line(size=1) +
  theme_light()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  theme(legend.position = "none")

# Create the plot
ggplot(test.df, aes(x = PC1, y = PC2, group = ID_plot, color = ID_plot)) +
  geom_path(size=1.25, arrow=arrow(length= unit(4, "mm"), type = "open")) +  
  geom_point() + # Plot points
 # geom_text(data = df_labels, aes(label = ID_plot, x = NMDS1, y = NMDS2), vjust = -1,
  #          check_overlap = T) + # Add group labels
  theme_bw() + 
  theme(legend.position = "none")
#labs(color = "Group") # Label the legend

## Model based ###
# fitnb <- gllvm(y = ordispe_pl, family = "negative.binomial") 
# fitg <- gllvm(y = ordispe_pl, family = "gaussian") 
# 
# AIC(fitnb)
# AIC(fitg)
# 
# ordiplot.gllvm(
#   fitnb,
#   main = NULL,
#   biplot = T,
#   s.colors = "black",
#   spp.colors = "blue",
#   alpha = 0.45,
#   cex.spp = 0.7,
#   jitter = F,
#   s.cex = 0.6
# )
# 
# abline(h = 0, v = 0, lty = 2)




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


