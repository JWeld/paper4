#DCA
#by subplot
DCA_sub <- decorana(ordispe, iweigh = 0)
pl <- ordiplot(DCA_sub, display = "sites", type="none")
points(pl, "sites", pch=21, col="grey60", bg="grey80")
ordispider (pl, groups = ordi_extras$ID_plot, col = c("red","blue","green","grey","pink","yellow"),
            label = F)
#extract scores
DCA_scores_sub <- vegan::scores(DCA_sub, choices = c(1,2))
DCA_scores_sub <- as.data.frame(DCA_scores_sub)
DCA_scores_sub <- DCA_scores_sub %>% rownames_to_column(.,var = "ID_fine2")

DCA_scores_sub.df <- left_join(DCA_scores_sub, extras, by = "ID_fine2")

#by plot
DCA_pl <- decorana(ordispe_pl, iweigh = 0)
ordiplot(DCA_pl)
pl <- ordiplot(DCA_pl)
ordiellipse (pl, groups = ordispe_pl_extras$ID_site,
             label = T, col = "red")
ordispider (pl, groups = ordispe_pl_extras$ID_plot, col = c("red","blue","green","grey","pink","yellow"),
            label = F)

#extract  vegan::scores
DCA_scores_pl <-  vegan::scores(DCA_pl, choices = c(1,2))
DCA_scores_pl <- as.data.frame(DCA_scores_pl)
DCA_scores_pl <- DCA_scores_pl %>% rownames_to_column(.,var = "ID_fine")

#combine DCA scores with other data####

combo_sites_dca <- left_join(DCA_scores_sub.df, select(extras, -c(ID_subplot)),
                             by = "ID_fine2") %>% distinct()

combo_sites_dca <- combo_sites_dca %>% select(ID_fine, ID.x, ID_site.x, ID_plot.x,
                                              survey_year.x, DCA1, DCA2) %>% 
  rename(ID=ID.x, ID_site = ID_site.x, ID_plot = ID_plot.x, survey_year=survey_year.x) %>% 
  distinct()

combo_sites_dca <- drop_na(combo_sites_dca) #CHECK THAT MORE CANT BE KEPT!
#year as index from first year
combo_sites_dca$year.i <- I(combo_sites_dca$survey_year-1997)

all_data_dca <- left_join(DCA_scores_sites2, IM_env, by = "ID")
#discard where deposition data is available
all_data_dca <- filter(all_data_dca, NH4M > 0)


#PLOTS- cup and ball models####
ggplot(DCA_scores_sub.df, aes(DCA1, DCA2, colour = ID_site)) + geom_point() +
  facet_wrap(facets = "ID_site")

ggplot(all_data_dca, aes(DCA1, DCA2, colour = ID_plot)) + geom_point() +
  facet_wrap(facets = "ID_plot")

ggplot(all_data_dca, aes(DCA1, DCA2,colour = ID_plot)) + geom_point()

ggplot(all_data_dca, aes(DCA1, DCA2,colour = ID_plot)) + geom_point()

#Figure
ggplot(drop_na(DCA_scores_sub.df), aes(DCA1, DCA2, colour = ID_site)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow())

test2 <- all_data_dca %>% group_by(ID_fine) %>% mutate(DCA1m = mean(DCA1)) %>% 
  mutate(DCA2m = mean(DCA2)) %>% ungroup()
ggplot(test2, aes(DCA1m, DCA2m, colour = ID_plot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow())
ggplot(all_data_dca, aes(DCA1, DCA2, colour = ID_subplot)) + geom_point()+ 
  geom_path(size=1.25, arrow=arrow())+ theme(legend.position="none")

#Figure####
ggplot(plot.level.dat %>% drop_na(pc_dist_base) %>%
         drop_na(ID_plot), aes(survey_year, pc_dist_base, colour = ID_plot)) +
  #  geom_point()+
  #  geom_smooth(method = "lm")+
  geom_line(size=1) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
  )





#ICP Forests data###
# env_means <- readRDS("~/Documents/r/paper4/Forests_data/env_means_2.RDS")
# veg <- readRDS("~/Documents/r/paper4/Forests_data/veg.RDS")
# veg_vasc <- readRDS("~/Documents/r/paper4/Forests_data/veg_vasc.RDS")
# 
# #Scale numeric variables####
# non_nums <- select_if(env_means, negate(is.numeric))
# nums <- as.data.frame(select_if(env_means, is.numeric) %>% scale())
# env_means_scaled <- bind_cols(non_nums, nums)
# 
# #create plot ID_site without year
# VF_vasc <- veg_vasc
# #VF_vasc <- mutate(VF_vasc, ID_site = substr(VF_vasc$ID, 6, 12))
# ID_site = substr(VF_vasc$ID, 6, 12)
# VF_vasc$ID_site <- ID_site
# 
# VF <- veg
# #VF_vasc <- mutate(VF_vasc, ID_site = substr(VF_vasc$ID, 6, 12))
# ID_site = substr(VF$ID, 6, 12)
# VF$ID_site <- ID_site
# 
# VF_vasc <- VF_vasc %>% ungroup() %>% 
#   group_by(survey_year, country, ID_site, ID, species_name) %>%
#   summarise(avg = mean(cover)) %>%
#   spread(species_name, avg) %>% ungroup()
# 
# VF <- VF %>% ungroup() %>% 
#   group_by(survey_year, country, ID_site, ID, species_name) %>%
#   summarise(avg = mean(cover)) %>%
#   spread(species_name, avg) %>% ungroup()
# 
# #we need the plots that have a time series...
# consistency <- as.data.frame(table(VF$ID_site))
# constant <- filter(consistency, Freq > 2)
# cons_list <- constant$Var1
# VF <- filter(VF, ID_site %in% cons_list)
# 
# consistency <- as.data.frame(table(VF_vasc$ID_site))
# constant <- filter(consistency, Freq > 2)
# cons_list <- constant$Var1
# VF_vasc <- filter(VF_vasc, ID_site %in% cons_list)
# 
# VF_vasc <- VF_vasc %>% replace(is.na(.), 0)
# VF_vasc <-  VF_vasc[rowSums(VF_vasc[, -c(1:4)]) != 0, ] 
# 
# VF <- VF %>% replace(is.na(.), 0)
# VF <-  VF[rowSums(VF[, -c(1:4)]) != 0, ] 
# 
# VF_vasc_m <- VF_vasc %>% select(.,-c(1:3))
# VF_m <- VF %>% select(.,-c(1:3))
# 
# VF_vasc_m <- column_to_rownames(VF_vasc_m, var="ID")
# VF_m <- column_to_rownames(VF_m, var="ID")
# 
# #environmental variables
# dat <- readRDS("Forests_data/dat.RDS")
# 




#-------------------------------------------------#
# 
# #Import VS (circular plots)####
# VS <- read_csv("raw_data/ICPIM_VS.csv",
#                        col_types = cols(DataProcessingDate = col_skip(),
#                                         DataQualityFlag = col_skip(),
#                                         Institute = col_skip(),
#                                         PFlag = col_skip(),
#                                         PrameterList = col_skip(),
#                                         SUBPROG = col_skip(),
#                                         StatusFlag = col_skip(),
#                                         TreeOrQuarterNumber = col_skip(),
#                                         Unit = col_skip()),
#     trim_ws = TRUE
#   )
# VS$Name <- NULL
# VS$SIZE <- NULL
# as.data.frame(table(VS$SpeciesName))
# as.data.frame(table(VS$StationCode))
# #drop 9999s for StationCode, not plot level
# VS <- filter(VS, StationCode != "9999")
# 
# #Don't need information on month just year. Add column for year only (but also keep yearmonth in case...)
# VS$survey_year <- as.integer(substr(VS$YearMonth, 0, 4))
# 
# 
# #add country
# VS <-mutate(VS, country = str_sub(VS$AreaCode, 1, 2)) #keep first two letters of country code and use to match
# VS$country <-
#   VS$country %>% str_replace_all(
#     c(
#       "AT" = "Austria",
#       "DE" = "Germany",
#       "EE" = "Estonia",
#       "FI" = "Finland",
#       "IT" = "Italy",
#       "LT" = "Lithuania",
#       "LV" = "Latvia",
#       "SE" = "Sweden",
#       "CZ" = "Czech",
#       "ES" = "Spain",
#       "IE" = "Ireland"
#     )
#   )
# #create variable "ID"
# VS <- transform(VS,ID = paste0(VS$survey_year, sep = "_",  VS$AreaCode))
# #create variable "ID_fine"
# VS <- transform(VS,ID_fine = paste0(VS$ID, sep = "_",  VS$StationCode))
# #create ID_site
# VS <- VS %>% rename(ID_site = AreaCode)
# #add lat long
# VS <- left_join(VS, IM_locations, by = "ID_site")
# 
# VS <- filter(VS, SpeciesName != "NULL")
# VS$Class <- NULL
# VS$Sample_ID <- NULL
# VS$SpatialPool <- NULL
# VS$SpeciesName <- NULL
# VS$SpeciesList <- NULL
# VS$StationName <- NULL
# VS$Species_ID <- NULL
# VS$MEDIUM <- NULL
# VS$description <- NULL
# #filter for only parameters of interest to reduce size and keep canopy cover seperate
# #VS_canopy <-
# #  VS %>% filter(VS$Parameter %in% c("COVE_T", "COVE_T1", "COVE_T2"))
# 
# VS_ground <-
#   VS %>% filter(VS$Parameter %in% c("COVE_F", "COVE_B"))
# VS_ground <- VS_ground %>% rename(Species = SpeciesCode, Cover = Value)
# VS_ground <- VS_ground %>% select(ID_site,ID,ID_fine, survey_year, country, Species, Cover,
#                                   latitude, longitude, everything())
# as.data.frame(table(VS_ground$ID_fine))
# as.data.frame(table(VS_ground$survey_year))
# 
# 
# #xXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# VS_2 <-
#   read_csv(
#     "raw_data/From_paper2/Veg_IM/data_from_Thomas/VS_Assessment_2012_AT_DE_ES_SE.csv",
#     col_types = cols(
#       CLASS = col_skip(),
#       Comment = col_skip(),
#       FLAGQUA = col_skip(),
#       FLAGSTA = col_skip(),
#       SCODE = col_character(),
#       INST = col_skip(),
#       LISTSPE = col_skip(),
#       MEDIUM = col_skip(),
#       PARLIST = col_skip(),
#       PFLAG = col_skip(),
#       SUBPROG = col_skip(),
#       TREE = col_skip(),
#       source_datafile = col_skip()
#     ),
#     trim_ws = TRUE
#   )
# IM_veg <- VS_2
# as.data.frame(table(IM_veg$SPECIES))
# 
# #Don't need information on month just year. Add column for year only (and keep month just in case...
# IM_veg$survey_year <- as.integer(substr(IM_veg$YYYYMM, 0, 4))
# IM_veg$LISTMED <- NULL
# 
# #add country
# IM_veg <-
#   mutate(IM_veg, country = str_sub(IM_veg$AREA, 1, 2)) #keep first two letters of country code and use to match
# IM_veg$country <-
#   IM_veg$country %>% str_replace_all(
#     c(
#       "AT" = "Austria",
#       "DE" = "Germany",
#       "EE" = "Estonia",
#       "FI" = "Finland",
#       "IT" = "Italy",
#       "LT" = "Lithuania",
#       "LV" = "Latvia",
#       "SE" = "Sweden",
#       "CZ" = "Czech",
#       "ES" = "Spain",
#       "IE" = "Ireland"
#     )
#   )
# #total covers have NA for species, separate out these
# IM_veg_totals <- filter(IM_veg,!is.na(SPECIES))
# IM_veg <- filter(IM_veg,!is.na(SPECIES))
# 
# #make new column ID (to year/site level)
# IM_veg <-
#   transform(IM_veg, ID = paste0(IM_veg$survey_year, sep = "_",  IM_veg$AREA))
# 
# #make new column ID_fine (to circle plot level)
# IM_veg <- IM_veg %>% transform(IM_veg, ID_fine = paste0(IM_veg$ID, sep = "_", IM_veg$SCODE))
# 
# #filter table for only parameters of interest to reduce size and keep canopy cover seperate
# IM_veg_cover <-
#   IM_veg %>% filter(IM_veg$PARAM %in% c("COVE_T", "COVE_T1", "COVE_T2", "COVE_T3"))
# IM_veg <- IM_veg %>% filter(IM_veg$PARAM %in% c("COVE_F", "COVE_B"))
# 
# IM_veg <- IM_veg %>% rename(Species = SPECIES, Cover = VALUE, ID_site = AREA)
# 
# #add location
# IM_veg <- left_join(IM_veg, IM_locations, by = "ID_site")
# 
# IM_veg <- IM_veg %>% select(ID_site, ID, ID_fine,PARAM, survey_year, country, Species, Cover, latitude, longitude)
# IM_veg <- rename(IM_veg, Parameter = PARAM)
# #Add location
# IM_veg <- left_join(IM_veg, IM_locations, by = "ID_site")
# 
# #XXXXXXXXXXXXXXXXXXXXXXX
# 
# 
# #double entries for each plot as combo of bottom and field, reduce to one per plot
# # VS_ground <- VS_ground %>% group_by(survey_year,country,ID,ID_site,ID_fine,Species) %>% summarise(Cover = sum(Cover)) %>% 
# #   ungroup()
# # 
# # VS_field <- VS_field %>% group_by(survey_year,country,ID,ID_site,ID_fine,Species) %>% summarise(Cover = sum(Cover)) %>% 
# #   ungroup()
# 
# #add lat long
# # VS_ground <- left_join(VS_ground, IM_locations, by = "ID_site")
# # VS_field <- left_join(VS_field, IM_locations, by = "ID_site")
# # VS_bot <- left_join(VS_bot, IM_locations, by = "ID_site")
# 
# #Combine the two ground cover dataframes
# VS_ground2 <-  full_join(VS_ground, IM_veg)
# VS_ground2 <- select(VS_ground2, survey_year, country, ID, ID_fine, ID_site,
#                      latitude, longitude, Species, Cover, Parameter)
# #double entries for each plot as combo of bottom and field, reduce to one per plot.Only a 
# #few species are found in both levels, but these are summed
#  VS_ground3 <- VS_ground2 %>% group_by(survey_year,country,ID,ID_site,ID_fine,Species) %>%
#    summarise(Cover = sum(Cover)) %>% 
#    ungroup()
# 
#  
#  #Field layer
#  VS_field <- filter(VS_ground2, Parameter == "COVE_F")
#  VS_field <- VS_field %>% select(ID, ID_site, ID_fine, country, survey_year, Species, Cover)
#  VS_field <- VS_field %>% group_by(survey_year,country,ID,ID_site,ID_fine,Species) %>%
#    summarise(Cover = sum(Cover)) %>% 
#    ungroup()
#  
#  #Bottom layer
#  VS_bot <- filter(VS_ground2, Parameter == "COVE_B")
#  VS_bot <- VS_bot %>% select(ID, ID_site, ID_fine, country, survey_year, Species, Cover)
#  VS_bot <- VS_bot %>% group_by(survey_year,country,ID,ID_site,ID_fine,Species) %>% summarise(Cover = sum(Cover)) %>% 
#    ungroup()
# # test <- VS_ground2 %>% group_by(ID_fine, Species) %>% summarise(Cover = mean(Cover))
# # test2 <- select(VS_ground2, survey_year, country, ID, ID_fine, ID_site, latitude, longitude) %>% distinct()
# # test3 <- left_join(test, test2, by = "ID_fine")
# # 
# # VS_ground3 <- test3
# # rm(test, test2, test3)
# 
# #Wide versions and matrices for vegan
# VS_g_wide <- VS_ground3 %>%
#   group_by(ID_fine) %>%
#   tidyr::pivot_wider(names_from = Species, values_from = Cover)
# VS_g_wide <- VS_g_wide %>% replace(is.na(.), 0) %>% ungroup()
# VS_g_wide <- VS_g_wide[rowSums(VS_g_wide[,-c(1:7)]) != 0, ]
# 
# #1961 dates for DE and lack of lat long for DE  XXXXXX
# 
# VS_b_wide <- VS_bot %>%
#   group_by(ID,ID_site,ID_fine) %>%
#   tidyr::pivot_wider(names_from = Species, values_from = Cover)
# VS_b_wide <- VS_b_wide %>% replace(is.na(.), 0) %>% ungroup()
# VS_b_wide <- VS_b_wide[rowSums(VS_b_wide[,-c(1:5)]) != 0, ]
# 
# VS_f_wide <- VS_field %>%
#   group_by(ID,ID_site,ID_fine) %>% 
#   tidyr::pivot_wider(names_from = Species, values_from = Cover) 
# VS_f_wide <- VS_f_wide %>% replace(is.na(.), 0) %>% ungroup()
# VS_f_wide <- VS_f_wide[rowSums(VS_f_wide[,-c(1:5)]) != 0, ]
# 
# #make matrix
# VS_gm <- select(VS_g_wide, -c(1:4))
# VS_gm <- column_to_rownames(VS_gm, var = "ID_fine")
# 
# VS_fm <- select(VS_f_wide, -c(1:4))
# VS_fm <- column_to_rownames(VS_fm, var = "ID_fine")
# 
# VS_bm <- select(VS_b_wide, -c(1:4))
# VS_bm <- column_to_rownames(VS_bm, var = "ID_fine")





#-------------------------------------------#