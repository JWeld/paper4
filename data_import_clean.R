library(tidyverse)
library(lubridate)
library(vegdata)
library(FD)
library(labdsv)
#library(vegan)
#library(ggfortify)
#library(dave)
#library(TR8)
#library(magrittr)

`%!in%` = Negate(`%in%`)

# Import location data ####
#Import IM locations (in decimal degrees)
IM_locations <- read_csv("raw_data/From_paper2/Veg_IM/IM_locations.csv")
IM_locations <-
  dplyr::select(IM_locations, Code, Latitude, Longitude) %>%
  rename(ID_site = Code,
         latitude = Latitude,
         longitude = Longitude)

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

#Import VG FD (Intensive plots and Forest damage) #1990 - 2018 ####
FD_VG <- read_csv("raw_data/ICPIM_FD_VG.csv", col_types = cols(Institute = col_skip()))
#break down to FD and VG
FD <- filter(FD_VG, Subprog == "FD")
VG <- filter(FD_VG, Subprog == "VG")
VG$Damage <- NULL
VG$StatusFlag <- NULL
VG$ParameterList <- NULL
VG$Subprog <- NULL
#Don't need information on month just year. Add column for year only (but also keep yearmonth in case...)
VG$survey_year <- as.integer(substr(VG$YearMonth, 0, 4))

#add country
VG <-mutate(VG, country = str_sub(VG$AreaCode, 1, 2)) #keep first two letters of country code and use to match
VG$country <-
  VG$country %>% str_replace_all(
    c(
      "AT" = "Austria",
      "DE" = "Germany",
      "EE" = "Estonia",
      "FI" = "Finland",
      "IT" = "Italy",
      "LT" = "Lithuania",
      "LV" = "Latvia",
      "SE" = "Sweden",
      "CZ" = "Czech",
      "ES" = "Spain",
      "IE" = "Ireland",
      "NO" = "Norway",
      "RU" = "Russia"
      
    )
  )
#drop Russia, only have data for 1 year so no change possible
VG <- filter(VG, country != "Russia")
#create variable "ID"
VG <- transform(VG,ID = paste0(VG$survey_year, sep = "_",  VG$AreaCode))
#create variable "ID_fine"
VG <- transform(VG,ID_fine = paste0(VG$ID, sep = "_",  VG$StationCode))
#create variable "ID_fine2"
VG <- transform(VG,ID_fine2 = paste0(VG$ID_fine, sep = "_",  VG$TreeOrQuarterNumber))

VG <- rename(VG, Quarter = TreeOrQuarterNumber)
VG <- rename(VG, ID_site = AreaCode)

#create variable "ID_subplot"
VG <- transform(VG,ID_subplot = paste0(VG$ID_site, sep = "_",  VG$StationCode,
                                       sep = "_",  VG$Quarter))
#create variable "ID_plot"
VG <- transform(VG,ID_plot = paste0(VG$ID_site, sep = "_",  VG$StationCode))

#One observation is very very nearly empty in field layer and is a crazy outlier on ordinations, remove
VG <- VG %>% filter(!ID_fine2 == "2000_DE01_60_220")

#there is a species named "NULL", remove
VG <- VG %>% filter(!Medium == "NULL")

#filter to field and bottom layers
VG_ground <- VG %>% filter(VG$Parameter %in% c("COVE_F", "COVE_B")) %>%
  rename(Param = Parameter)
VG_ground <- select(VG_ground, -c("AreaName", "StationName", "Descrption", "country", "SpeciesName", "Class"))
VG_ground$Class <- NULL
VG_ground$ListMedium <- NULL
VG_ground$Unit <- NULL
VG_ground$Quarter <- NULL

VG_field <- filter(VG_ground, Param == "COVE_F")
VG_field <- VG_field %>% rename(Species = Medium, Cover = Value)

VG_bot <- filter(VG_ground, Param == "COVE_B")
VG_bot <- VG_bot %>% rename(Species = Medium, Cover = Value)

VG_ground$Param <- NULL
#some species occur in both F and B levels, sum these
VG_ground <- VG_ground %>% group_by(ID_fine2, Medium) %>%
  summarise(Cover = sum(Value),
            survey_year = unique(survey_year),
            ID =unique(ID),
            ID_site=unique(ID_site),
            ID_plot=unique(ID_plot),
            ID_subplot=unique(ID_subplot),
            ID_fine=unique(ID_fine)) %>% rename(Species = Medium) %>% ungroup()

VG_b_wide <- VG_bot %>%
  group_by(ID,ID_site,ID_fine2) %>%
  #mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Cover)# %>%
  #select(-row)
VG_b_wide <- VG_b_wide %>% replace(is.na(.), 0) %>% ungroup()

VG_f_wide <- VG_field %>%
  group_by(ID,ID_fine,ID_fine2) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Cover) 
VG_f_wide <- VG_f_wide %>% replace(is.na(.), 0) %>% ungroup()

VG_g_wide <- VG_ground %>% group_by(ID,ID_fine,ID_fine2) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Cover)
VG_g_wide <- VG_g_wide %>% replace(is.na(.), 0) %>% ungroup()

VG_g_wide$ID_plot <- as.factor(VG_g_wide$ID_plot)
VG_g_wide$ID_subplot <- as.factor(VG_g_wide$ID_subplot)
# VG_g_wide2 <- left_join(VG_g_wide, select(VG_ground,ID_fine2,ID_plot,ID_subplot),by = "ID_fine2") %>% 
#   select(., ID_plot, ID_subplot, everything())

#Try grouping by intensive plot/year combo rather than individual small squares####
VG_ground2 <- VG_ground %>% group_by(survey_year,ID,ID_site,ID_fine, Species) %>% 
  summarise(Cover = mean(Cover)) %>% ungroup()

VG_g_wide2 <- VG_ground2 %>% group_by(ID,ID_fine) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Cover)
VG_g_wide2 <- VG_g_wide2 %>% replace(is.na(.), 0) %>% ungroup()

#find sites that have no data for subplots (are they reporting overall mean??)
#few_f <- VG_f_wide %>% group_by(ID_fine) %>% summarise(distinctID2 = n_distinct(ID_fine2)) %>% filter(distinctID2 < 2)
#few_b <- VG_b_wide %>% group_by(ID_fine) %>% summarise(distinctID2 = n_distinct(ID_fine2)) %>% filter(distinctID2 < 2)
#few_g <- VG_g_wide %>% group_by(ID_fine) %>% summarise(distinctID2 = n_distinct(ID_fine2)) %>% filter(distinctID2 < 2)

#remove, can't calculate DCA from them
#VG_f_wide <- filter(VG_f_wide, ID_fine %!in% few_f$ID_fine)
#VG_b_wide <- filter(VG_b_wide, ID_fine %!in% few_b$ID_fine)
#VG_g_wide <- filter(VG_g_wide, ID_fine %!in% few_g$ID_fine)

#make factor?
VG_f_wide$ID <- as_factor(VG_f_wide$ID)
VG_f_wide$ID_site <- as_factor(VG_f_wide$ID_site)
VG_f_wide$ID_fine <- as_factor(VG_f_wide$ID_fine)
VG_f_wide$ID_fine2 <- as_factor(VG_f_wide$ID_fine2)

VG_b_wide$ID <- as_factor(VG_b_wide$ID)
VG_b_wide$ID_site <- as_factor(VG_b_wide$ID_site)
VG_b_wide$ID_fine <- as_factor(VG_b_wide$ID_fine)
VG_b_wide$ID_fine2 <- as_factor(VG_b_wide$ID_fine2)

VG_g_wide$ID <- as_factor(VG_g_wide$ID)
VG_g_wide$ID_site <- as_factor(VG_g_wide$ID_site)
VG_g_wide$ID_fine <- as_factor(VG_g_wide$ID_fine)
VG_g_wide$ID_fine2 <- as_factor(VG_g_wide$ID_fine2)

#later functions complain about zero sum rows, check and remove
VG_f_wide <- VG_f_wide[rowSums(VG_f_wide[,-c(1:11)]) != 0, ]
VG_b_wide <- VG_b_wide[rowSums(VG_b_wide[,-c(1:11)]) != 0, ]
VG_g_wide <- VG_g_wide[rowSums(VG_g_wide[,-c(1:7)]) != 0, ]

VG_f_wide$StationCode <- NULL
VG_f_wide$SpatialPool <- NULL
VG_f_wide$Param <- NULL

VG_b_wide$StationCode <- NULL
VG_b_wide$SpatialPool <- NULL
VG_b_wide$Param <- NULL

VG_g_wide$StationCode <- NULL
VG_g_wide$SpatialPool <- NULL
VG_g_wide$Param <- NULL
# DCA_func <- function(x){y <- select(x,-c(1:10))
#   decorana(veg = y)
#   
# }

VG_fm <- column_to_rownames(VG_f_wide, var = "ID_fine2")
VG_fm <- select(VG_fm, -c(1:7))

VG_bm <- column_to_rownames(VG_b_wide, var = "ID_fine2")
VG_bm <- select(VG_bm, -c(1:7))

VG_gm <- column_to_rownames(VG_g_wide, var = "ID_fine2")
VG_gm <- select(VG_gm, -c(1:6))

#And the same for the per intensive plot version
#make factor
VG_g_wide2$ID <- as_factor(VG_g_wide2$ID)
VG_g_wide2$ID_site <- as_factor(VG_g_wide2$ID_site)
VG_g_wide2$ID_fine <- as_factor(VG_g_wide2$ID_fine)

#later functions complain about zero sum rows, check and remove
VG_g_wide2 <- VG_g_wide2[rowSums(VG_g_wide2[,-c(1:4)]) != 0, ]


VG_gm2 <- select(VG_g_wide2, -c(1:3))
VG_gm2 <- column_to_rownames(VG_gm2, var = "ID_fine")






# split.df <- VG_f_wide %>% column_to_rownames(.,var = "ID_fine2") %>% 
#   group_by(ID_fine) %>% .[rowSums(.[,-c(1:10)]) != 0, ] %>%
#   .[complete.cases(.), ] %>% group_split() 

#run by site ordinations, probably not what we want
#split.df <- split(VG_f_wide, VG_f_wide$ID_fine)

#z <-lapply(split.df, function(x){y <- select(x,-c(1:10)) #%>% select_if(colSums(.) != 0)
#d <- decorana(veg = y)
#}
#)

#names(z)
#z$`1991_DE01_50`$evals #eigenvalues
#z$`1991_DE01_50`$evals.decorana #decorona values

#Environmental Variables#####
ICPIM_AC_AM_TF <- read_csv("raw_data/ICPIM_AC_AM_TF.csv", 
                           col_types = cols(Day = col_skip(), DataQualityFlag = col_skip(), 
                                            StatusFlag = col_skip(), NeedleAgeCode = col_skip(), 
                                            ParameterInfo = col_skip(), PretreInfo = col_skip(), 
                                            DeterInfo = col_skip()))

#Add arbitrary day to yearmonth to create date data
ICPIM_AC_AM_TF$YearMonth <- paste(ICPIM_AC_AM_TF$YearMonth,"15", sep="")
#change to date
ICPIM_AC_AM_TF$YearMonth <- ymd(ICPIM_AC_AM_TF$YearMonth)
#Add column for year only
ICPIM_AC_AM_TF$survey_year <- year(ICPIM_AC_AM_TF$YearMonth) 
#rename
ICPIM_AC_AM_TF <- ICPIM_AC_AM_TF %>% rename(ID_site = AreaCode)
#add country
ICPIM_AC_AM_TF <-
  mutate(ICPIM_AC_AM_TF, country = str_sub(ICPIM_AC_AM_TF$ID_site, 1, 2)) #keep first two letters of country code and use to match
ICPIM_AC_AM_TF$country <-
  ICPIM_AC_AM_TF$country %>% str_replace_all(
    c(
      "AT" = "Austria",
      "DE" = "Germany",
      "DK" = "Denmark",
      "EE" = "Estonia",
      "FI" = "Finland",
      "IT" = "Italy",
      "LT" = "Lithuania",
      "LV" = "Latvia",
      "NO" = "Norway",
      "PL" = "Poland",
      "RU" = "Russia",
      "SE" = "Sweden",
      "CZ" = "Czech",
      "ES" = "Spain",
      "IE" = "Ireland"
    )
  )

#Make column ID
ICPIM_AC_AM_TF <-
  transform(
    ICPIM_AC_AM_TF,
    ID = paste0(
      ICPIM_AC_AM_TF$survey_year,
      sep = "_",
      ICPIM_AC_AM_TF$ID_site
    )
  )

#Air Chemistry
IM_AC <- filter(ICPIM_AC_AM_TF, Subprog == "AC")
#Meterology- temperature
IM_AM <- filter(ICPIM_AC_AM_TF, Subprog == "AM") %>%
  filter(ParameterCode == "TEMP")
#Throughfall deposition N, S, precip
IM_TF <- filter(ICPIM_AC_AM_TF, Subprog == "TF") %>% 
  filter(ParameterCode == "NH4N" | ParameterCode == "NO3N" |
           ParameterCode == "SO4S"| ParameterCode == "PREC")

IM_TF <- select(IM_TF,ID, ID_site,survey_year, country, StationCode,
                YearMonth, ParameterCode, Value)
#take monthly mean for variables with more than one measurement in a month recorded
IM_TF <- IM_TF %>% group_by(ID,country,ID_site,survey_year, YearMonth,StationCode, ParameterCode) %>%
  summarise_at(vars(Value), funs(mean), na.rm = TRUE) 

#IM data is long, make wide version
IM_TF <- pivot_wider(IM_TF, names_from = ParameterCode, values_from = Value)
#monthly precip x TF conc
IM_dep <- IM_TF %>% 
  mutate(NH4M = (NH4N*PREC)*0.01) %>% 
  mutate(NO3M = (NO3N*PREC)*0.01) %>% 
  mutate(SO4SM = (SO4S*PREC)*0.01)
# average per site/year combo. First find annual total for each station and then take 
# mean of station values to get a whole site value.
IM_dep_means <- IM_dep %>% group_by(ID,ID_site,StationCode,survey_year) %>%
  summarise_at(vars(NH4M, NO3M, SO4SM, PREC), funs(sum), na.rm = TRUE) 

IM_dep_means <- IM_dep_means %>% group_by(ID,ID_site,survey_year) %>%
  summarise_at(vars(NH4M, NO3M, SO4SM, PREC), funs(mean), na.rm = TRUE) 

IM_dep2 <- drop_na(IM_dep_means)
IM_dep2$survey_year <- NULL

#Temp
IM_AM <- select(IM_AM,ID, ID_site,survey_year, country, StationCode,
                YearMonth, ParameterCode, Value)
#take monthly mean for variables with more than one measurement in a month recorded
IM_AM <- IM_AM %>% group_by(ID,country,ID_site,survey_year, YearMonth,StationCode,
                            ParameterCode) %>%
  summarise_at(vars(Value), funs(mean), na.rm = TRUE) 
IM_AM <- pivot_wider(IM_AM, names_from = ParameterCode, values_from = Value)
#IM_AM is monthly if you need eg summer vs winter avergae temps
#make annual mean
IM_meteo <- IM_AM %>% group_by(ID,ID_site,survey_year) %>%
  summarise_at(vars(TEMP), funs(mean), na.rm = TRUE) %>% ungroup()

#ENV combo####
IM_env <- left_join(IM_dep2, select(IM_meteo, ID, TEMP), by = "ID") %>% ungroup() %>% 
  left_join(., IM_locations) %>% 
  select(., -c(ID_site))


#create data - species and environemntal matrixes####
#There is species data back to 90 but only env data from 98 on...
test <- left_join(VG_g_wide, IM_env, by = "ID") %>% filter(NH4M > 0) %>% ungroup()

#date range for plots
last <- test %>% group_by(ID_plot) %>% summarise(last=max(survey_year))
first <- test %>% group_by(ID_plot) %>% summarise(first=min(survey_year))
range <- left_join(first, last)
range <- range %>% mutate(diff= last-first)

#German and Italian have only one year of data, drop them
droplist <- c()
test <- filter(test, ID_plot %!in% c("DE01_50", "DE01_60", "IT03_1"))

#make matrix for ordinations
ordienv <-  select(ungroup(test), ID_fine2, NH4M, NO3M, SO4SM, TEMP, latitude, longitude, 
                   PREC, survey_year) %>% column_to_rownames(var="ID_fine2")
ordispe <- select(test, -c(2:7, 331:337)) %>% column_to_rownames(var="ID_fine2")
ordienv <-  select(test, ID_fine2, NH4M, NO3M, SO4SM, latitude, longitude, 
                   PREC, survey_year) %>% column_to_rownames(var="ID_fine2")
ordispe <- as.matrix(ordispe)
ordienv <- as.matrix(ordienv)


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
