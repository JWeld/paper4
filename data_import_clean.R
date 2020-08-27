library(tidyverse)
library(lubridate)
library(vegdata)
library(FD)
library(labdsv)
#library(vegan)
#library(ggfortify)
#library(dave)
library(TR8)
library(magrittr)


#Import VS (circular plots)####
VS <- read_csv("raw_data/ICPIM_VS.csv", 
                       col_types = cols(DataProcessingDate = col_skip(), 
                                        DataQualityFlag = col_skip(),
                                        Institute = col_skip(), 
                                        PFlag = col_skip(),
                                        PrameterList = col_skip(), 
                                        SUBPROG = col_skip(),
                                        StatusFlag = col_skip(),
                                        TreeOrQuarterNumber = col_skip(),
                                        Unit = col_skip()),
    trim_ws = TRUE
  )
VS$Name <- NULL
VS$SIZE <- NULL
as.data.frame(table(VS$SpeciesName))

#Don't need information on month just year. Add column for year only (but also keep yearmonth in case...)
VS$survey_year <- as.integer(substr(VS$YearMonth, 0, 4))


#add country
VS <-mutate(VS, country = str_sub(VS$AreaCode, 1, 2)) #keep first two letters of country code and use to match
VS$country <-
  VS$country %>% str_replace_all(
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
      "IE" = "Ireland"
    )
  )
#create variable "ID"
VS <- transform(VS,ID = paste0(VS$survey_year, sep = "_",  VS$AreaCode))

#filter for only parameters of interest to reduce size and keep canopy cover seperate
VS_canopy <-
  VS %>% filter(VS$Parameter %in% c("COVE_T", "COVE_T1", "COVE_T2"))

VS_under <- 
  VS %>% filter(VS$Parameter %in% c("COVE_F", "COVE_B"))

VS_shrub <- 
  VS %>% filter(VS$Parameter == "COVE_S")

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

#create variable "ID"
VG <- transform(VG,ID = paste0(VG$survey_year, sep = "_",  VG$AreaCode))
#create variable "ID_fine"
VG <- transform(VG,ID_fine = paste0(VG$ID, sep = "_",  VG$StationCode))

# Import location data ####
#Import IM locations (in decimal degrees)
IM_locations <- read_csv("raw_data/From_paper2/Veg_IM/IM_locations.csv")
IM_locations <-
  dplyr::select(IM_locations, Code, Latitude, Longitude) %>%
  rename(ID_site = Code,
         latitude = Latitude,
         longitude = Longitude)

#Import VS veg data from paper 2 to compare- are we missing anything in the latest data####
#(this file includes the extra data from Thomas but only up to 2012
VS_Assessment_2012 <-
  read_csv(
    "raw_data/From_paper2/Veg_IM/data_from_Thomas/VS_Assessment_2012_AT_DE_ES_SE.csv",
    col_types = cols(
      CLASS = col_skip(),
      Comment = col_skip(),
      FLAGQUA = col_skip(),
      FLAGSTA = col_skip(),
      INST = col_skip(),
      LISTSPE = col_skip(),
      MEDIUM = col_skip(),
      PARLIST = col_skip(),
      PFLAG = col_skip(),
      SCODE = col_character(),
      SUBPROG = col_skip(),
      TREE = col_skip(),
      source_datafile = col_skip()
    ),
    trim_ws = TRUE
  )


IM_VS12 <- VS_Assessment_2012
rm(VS_Assessment_2012)
as.data.frame(table(IM_VS12$SPECIES))
length(complete.cases(IM_VS12$SPECIES) == TRUE)

#Don't need information on month just year. Add column for year only (and keep month just in case...
IM_VS12$survey_year <- as.integer(substr(IM_VS12$YYYYMM, 0, 4))
IM_VS12$LISTMED <- NULL

#add country
IM_VS12 <-
  mutate(IM_VS12, country = str_sub(IM_VS12$AREA, 1, 2)) #keep first two letters of country code and use to match
IM_VS12$country <-
  IM_VS12$country %>% str_replace_all(
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
      "IE" = "Ireland"
    )
  )

#filter table for only parameters of interest to reduce size and keep canopy cover seperate
IM_VS12_cover <-
  IM_VS12 %>% filter(IM_VS12$PARAM %in% c("COVE_T", "COVE_T1", "COVE_T2", "COVE_T3"))
IM_VS12 <- IM_VS12 %>% filter(IM_VS12$PARAM %in% c("COVE_F", "COVE_B"))

#add location
IM_locations <- IM_locations %>% mutate(AREA = ID_site)
IM_VS12 <- left_join(IM_VS12, IM_locations, by = "AREA")
#create variable "ID"
IM_VS12 <- transform(IM_VS12,ID = paste0(IM_VS12$survey_year, sep = "_",  IM_VS12$AREA))
#create variable "ID_fine"
IM_VS12 <- transform(IM_VS12,ID_fine = paste0(IM_VS12$ID, sep = "_",  IM_VS12$SCODE))

spp_means <- IM_VS12 %>% 
  group_by(survey_year, ID_fine, SPECIES) %>%
  summarise(avg = mean(VALUE)) %>%
  spread(SPECIES, avg)


#Import VG data from paper 2####
#Import/clean paper 2 data version of VG (intensive plots) # 1982 -1999 only
IM_VG_data <- read_csv("raw_data/From_paper2/Veg_IM/IM_VG_data.csv", 
                       col_types = cols(Inst = col_skip(), Method = col_skip(), 
                                        ParamList = col_skip(), QualityFlag = col_skip(), 
                                        Spool = col_double(), StatusFlag = col_skip(), 
                                        Subprog = col_skip(), level = col_skip()))

#Thomas' data # 1982-2012
VG_Assessment_2012 <- read_csv("raw_data/From_paper2/Veg_IM/data_from_Thomas/VG_Assessment_2012.csv", 
                               col_types = cols(Comment = col_skip(), 
                                                INST = col_skip(), PARLIST = col_skip(), FLAGSTA = col_skip(),
                                                SUBPROG = col_skip(), source_datafile = col_skip()))
#fix colnames
VG_Assessment_2012 <- rename(VG_Assessment_2012, AreaCode=AREA, StationCode=SCODE, Medium=MEDIUM, ListMedium=LISTMED, Spool=SPOOL,
       Param=PARAM, value=VALUE, Param=PARAM, Value=VALUE, Unit=UNIT, Class=CLASS)
IM_VG_data <- rename(IM_VG_data, AreaCode=Area, Medium=Species)

#create variable "ID"
IM_VG_data$survey_year <- as.integer(substr(IM_VG_data$YYYYMM, 0, 4))
IM_VG_data <- transform(IM_VG_data,ID = paste0(IM_VG_data$survey_year, sep = "_",  IM_VG_data$AreaCode))
VG_Assessment_2012$survey_year <- as.integer(substr(VG_Assessment_2012$YYYYMM, 0, 4))
VG_Assessment_2012 <- transform(VG_Assessment_2012,ID = paste0(VG_Assessment_2012$survey_year, sep = "_",
                                                               VG_Assessment_2012$AreaCode))

#filter to field and bottom layers
VG_P2_ground <- IM_VG_data %>% filter(IM_VG_data$Param %in% c("COVE_F", "COVE_B"))
VG_P2T_ground <- VG_Assessment_2012 %>% filter(VG_Assessment_2012$Param %in% c("COVE_F", "COVE_B"))
VG_ground <- VG %>% filter(VG$Parameter %in% c("COVE_F", "COVE_B")) %>% rename(Param = Parameter)



spp_means <- VG_ground %>% 
  group_by(survey_year, ID, Medium) %>%
  summarise(avg = mean(Value)) %>%
  spread(Medium, avg)
