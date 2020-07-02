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


#Import VS
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

#Import VG FD
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

#filter table for only parameters of interest to reduce size and keep canopy cover seperate
VS_canopy <-
  VS %>% filter(VS$Parameter %in% c("COVE_T", "COVE_T1", "COVE_T2"))

VS_under <- 
  VS %>% filter(VS$Parameter %in% c("COVE_F", "COVE_B"))

VS_shrub <- 
  VS %>% filter(VS$Parameter == "COVE_S")

