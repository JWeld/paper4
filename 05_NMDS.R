library(tidyverse)
library(vegan)
options(mc.cores = parallel::detectCores())
library(gvlma)#checks assumptions for models

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
#NMDS####
#subplots #ordispe_sub
nmds_sub <- metaMDS(ordispe, distance = "bray", k=2, parallel = 8, try =20, trymax = 25)
#ordiplot(nmds_sub)
#stressplot(nmds_sub)
pl <- ordiplot(nmds_sub, display = "sites", type="none")
text(pl, "sites", cex=0.5, col="grey60", bg="grey80")
ordiellipse (pl, groups = ordi_extras$ID_plot,
             label = T, col = "red")

#Weird outliers for 2016_SE14_0001_24? and 2016_SE14_0001_5?
check <- rownames_to_column(ordispe_sub, var="ID")
check2 <- filter(check, ID !="2016_SE14_0001_24")
check2 <- filter(check2, ID !="2016_SE14_0001_5")
check2 <- column_to_rownames(check2, "ID")

ordi_extras_check <- filter(ordi_extras.df, ID_fine2 !="2016_SE14_0001_24")
ordi_extras_check <- filter(ordi_extras_check, ID_fine2 !="2016_SE14_0001_5")

nmds_sub <- metaMDS(check2, distance = "bray", k=2, parallel = 8, try =20, trymax = 25)
#ordiplot(nmds_sub)
#stressplot(nmds_sub)
pl <- ordiplot(nmds_sub, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordiellipse (pl, groups = ordi_extras_check$ID_plot,
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
NMDS_pl_sites_scores <- left_join(NMDS_pl_sites_scores, extras, by = "ID_fine")

# NMDS_pl_sites_scores <- left_join(NMDS_pl_sites_scores, select(extras, -c(ID_fine2, ID_subplot)), by = "ID_fine") %>% 
#   distinct() %>%
#   filter(!ID_site %in% c("AT01", "CZ02", "IT01","IT02",
#                          "IT03","IT04","IT05","IT06",
#                          "IT07","IT08","IT09","IT10",
#                          "IT11","IT12", "LT03", "ES02",
#                          "DE02", "RU04", "RU13"))

NMDS_pl_sp_scores <- NMDS_pl_scores$species

ggplot(NMDS_sub_sites_scores, aes(NMDS1, NMDS2, colour = ID_plot)) + geom_point()+ 
  geom_path(linewidth=0.5, arrow=arrow()) +
  theme_light()

df_labels <- aggregate(cbind(NMDS1, NMDS2) ~ ID_plot, data = NMDS_pl_sites_scores, FUN = mean)

# Create the plot
ggplot(filter(NMDS_sub_sites_scores, NMDS1 > -10),aes(x = NMDS1, y = NMDS2, group = ID_subplot, color = ID_subplot)) +
  geom_path(linewidth=0.5, arrow=arrow(length= unit(4, "mm"), type = "open")) +  # Connect points within each group
  geom_point() + # Plot points
  #geom_text(data = df_labels, aes(label = ID_plot, x = NMDS1, y = NMDS2), vjust = 0,
  #         check_overlap = T) + # Add group labels
  theme_bw() +  # Optional: Use a minimal theme for the plot
  theme(legend.position = "none")+
  labs(color = "Group") # Label the legend

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