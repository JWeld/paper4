`%!in%` = Negate(`%in%`)

library(tidyverse)
library(lubridate)
library(vegdata)
library(FD)
library(labdsv)
library(vegan)
library(ggfortify)
library(dave)
library(TR8)
library(magrittr)

#Tree species metadata
FOR_tree_species <-
  read.csv("raw_data/FOR_tree_species.csv") %>%  #tree species codes
  rename(code_tree_species = code, species_name = description) %>% #standardise column names
  dplyr::select(-4, -5) #keep only relevant columns
FOR_tree_species$species_name  <-
  gsub("[[:punct:]]", "", FOR_tree_species$species_name) #remove asterisks
FOR_tree_species$species_name  <-
  word(FOR_tree_species$species_name, 1, 2) #cut to genus and species

#Country metadata
FOR_country <-
  read.csv("raw_data/FOR_country.csv") %>%  #country codes
  rename(code_country = code, country = lib_country) %>% #standardise column names
  dplyr::select(code_country, country)#keep only relevant columns

#Import location data ----------------------------------------------------
library(magrittr)
#Import Forest Veg locations (in degrees minutes seconds)
FOR_veg_locations <-
  read_csv("raw_data/Veg_Forests/FOR_veg_locations.csv")
#remove canary islands
FOR_veg_locations <- filter(FOR_veg_locations, longitude >= -99999)
temp <- FOR_veg_locations
temp$longitude <- as.character(temp$longitude)
temp$latitude <- as.character(temp$latitude)
temp$longitude <- gsub("-", "W 0", FOR_veg_locations$longitude)#add leading zero for minus numbers
#and W instead of minus. After Canaries removal all are single digit degrees W.
temp1 <- filter(temp, str_detect(temp$longitude,"W", negate = TRUE)) 
temp2 <- filter(temp, str_detect(temp$longitude,"W", negate = FALSE)) 
temp1$longitude <- str_pad(temp1$longitude, 6, 'left', '0')
temp1$longitude <- paste("E ", temp1$longitude, sep="")
temp <- full_join(temp1,temp2)
temp$longitude <- str_pad(temp$longitude, 8, 'right', '0')
test <- tidyr::separate(temp, longitude, c('EW','deg','min','sec'), sep = c(2,4,6), remove = FALSE)
test$EW <- str_trim(test$EW)
test <- unite(test, "long2", deg, min, sec, EW, sep="/")

#and insert divider to latitudes. Easier as all positive and in same range
test <- tidyr::separate(test, latitude, c('deg','min','sec'), sep = c(2,4), remove = FALSE)
test <- unite(test, "lat2", deg, min, sec, sep="/")
test$lat2 <-  paste0(test$lat2,"/N")
test$latitude <- NULL
test$longitude <- NULL
test <- test %>% rename(longitude=long2, latitude=lat2)
#create plot ID_site without year
test <- mutate(test, ID_site = substr(test$ID, 6, 12))
test <- mutate(test, year = substr(test$ID, 0, 4))
test <-
  test %>% group_by(ID_site) %>% summarise_at(c("latitude", "longitude"), max)

dms2dec <- function(dms, separators = c("/", "/", "/'", "/")) {
  
  # dms: a vector (or column) of latitude or longitude in degrees-minutes-seconds-hemisfere,
  #e.g. 41° 34' 10.956" N (with or without spaces)
  # separators: the characters that are separating degrees, minutes and seconds in dms
  
  dms <- as.character(dms)
  dms <- gsub(pattern = " ", replacement = "", x = dms)
  for (s in separators) dms <- gsub(pattern = s, replacement = "_splitHere_", x = dms)
  
  splits <- strsplit(dms, split = "_splitHere_")
  n <- length(dms)
  deg <- min <- sec <- hem <- vector("character", n)
  
  for (i in 1:n) {
    deg[i] <- splits[[i]][1]
    min[i] <- splits[[i]][2]
    sec[i] <- splits[[i]][3]
    hem[i] <- splits[[i]][4]
  }
  
  dec <- as.numeric(deg) + (as.numeric(min) / 60) + (as.numeric(sec) / 3600)
  sign <- ifelse (hem %in% c("N", "E"), 1, -1)
  dec <- sign * dec
  return(dec)
}  # end dms2dec function 

test$longitude <-  dms2dec(test$longitude)
test$latitude <-  dms2dec(test$latitude)

FOR_veg_locations <- test

#Import IM locations (in decimal degrees)
IM_locations <- read_csv("raw_data/Veg_IM/IM_locations.csv")
IM_locations <-
  dplyr::select(IM_locations, Code, Latitude, Longitude) %>%
  rename(ID_site = Code,
         latitude = Latitude,
         longitude = Longitude)

location <- full_join(FOR_veg_locations, IM_locations)

#Meterological data####

#library(raster) #claims select from dplyr!
#library(sp)
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

#higher resolution
#r <- getData("worldclim",var="bio",res=2.5) #uncomment for first run!!!XXXXXX
#r <- r[[c(1,12)]]
 #names(r) <- c("mean_temp", "mean_precip")
 # points <- cbind.data.frame(location$ID_site, location$latitude, location$longitude)
 # points <- points[!duplicated(points), ]
 #  values <- extract(r,points[,2:3])
 #  meteo <- cbind.data.frame(points,values)
 #  colnames(meteo) <-  c("ID_site", "latitude", "longitude", "temp", "precip")

# Create a data.frame with sample site coordinates
 # #site <- veg$ID
 # site <- location$ID_site
 # #lon <- veg$longitude
 # #lat <- veg$latitude
 # lon <- location$longitude
 # lat <- location$latitude
 # samples <- data.frame(site, lon, lat)
 # points2 <- samples[!duplicated(samples), ]
 # rownames(points2) <- points2$site
 # points2$site <- NULL
 # values2 <- extract(r,points2)
 # meteo2 <- cbind.data.frame(points2,values2)
 # meteo2$ID_site <- rownames(meteo2)
 # rownames(meteo2) <- NULL
 # colnames(meteo2) <-  c("longitude","latitude", "mean_temp", "mean_precip","ID_site")
 # meteo2$mean_temp <- meteo2$mean_temp/10
 # saveRDS(meteo2, file = "meteo2.RDS")
meteo2 <-  readRDS("output_data/meteo2.RDS")
meteo2$mean_temp <- meteo2$mean_temp/10 #check this if reloading, see line 555

#import ICP Forests meteo data
mm_mem <- read_delim("raw_data/311_mm_20181120123705/mm_mem.csv", 
                     ";", escape_double = FALSE, col_types = cols(date_observation = col_date(format = "%Y-%m-%d")), 
                     trim_ws = TRUE)

mm2 <- filter(mm_mem, code_variable %in% c("PR", "AT"))
rm(mm_mem)
mm <- dplyr::select(mm2, survey_year, code_country, code_variable, date_observation,
                    daily_mean, code_plot, line_nr)

mm <-
  transform(
    mm,
    ID = paste0(
      mm$survey_year,
      sep = "_",
      mm$code_country,
      sep = "_",
      mm$code_plot
    )
  )

mm <-
  transform(
    mm,
    ID_site = paste0(
      mm$code_country,
      sep = "_",
      mm$code_plot
    )
  )
rm(mm2)

library(lubridate)
mm$year <- floor_date(mm$date_observation, "year")

AT <- dplyr::filter(mm, code_variable == "AT") %>% dplyr::select(ID,ID_site, daily_mean) %>% drop_na()
PR <- dplyr::filter(mm, code_variable == "PR") %>% dplyr::select(ID, ID_site, date_observation, daily_mean) #%>% drop_na()

AT <- AT %>% filter(daily_mean >-50) %>% filter(daily_mean<50)
PR <- PR %>% filter(daily_mean >=0) %>% filter(daily_mean<600)

#filter out plots with only few temperature observations in a year to avoid bias and use
#mean annual from area instead
remove <- AT %>% count(ID)  %>% filter(n<300)
remove <- remove$ID
AT <- filter(AT, ID %!in% remove)
AT <- AT %>%
  group_by(ID) %>%
  summarise(temp = mean(daily_mean)) #%>% summary(temp) #%>% filter(temp >-1)

#monthly precip data
PR_m <- PR %>% group_by(ID,ID_site,month1=floor_date(date_observation, "month")) %>%
  summarize(monthly_total=sum(daily_mean))
PR_m$month <- month(PR_m$month1)
PR_m <-
  transform(
    PR_m,
    IDm = paste0(
      PR_m$ID,
      sep = "_",
      PR_m$month
    )
  )
PR_m <- filter(PR_m, ID != "2002_8_11") #wrong date, 1960

#annual precip data
PR_ann <- PR %>%
  group_by(ID,ID_site) %>%
  add_tally(daily_mean) %>% dplyr::select(ID, n) %>% distinct() %>% 
  filter(n<8000) %>% rename(precip=n) %>% filter(precip>0)
# filter(n<3000) %>% rename(precip=n) %>% filter(precip>200)

meteo_ann <- full_join(AT,PR_ann)
length(unique(meteo2$ID))
test2 <- full_join(meteo_ann, meteo2)
#check deviations from the mean are plausible
test2 <- test2 %>% mutate(tdiff = mean_temp - temp)
test2 <- test2 %>% mutate(pdiff = mean_precip - precip)
test2 <- test2 %>% mutate(pdiffp = precip/mean_precip)
#quantile(test2$pdiffp, .99, na.rm = TRUE)#2.79
test2 <- test2 %>% mutate(temp = replace(temp, tdiff > 5 | tdiff < -5, NA))
test2 <- test2 %>% mutate(precip = replace(precip, pdiffp > 5 | pdiffp < 0.03, NA))
test2$precip_c <- coalesce(test2$precip,test2$mean_precip)
test2$temp_c <- coalesce(test2$temp,test2$mean_temp)
#test2$temp_c <- coalesce(test2$temp,test2$mean_temp/10)
#combine, using direct data from plots where possible but filling gaps with values for 
#that grid square
meteo3 <- test2 %>% dplyr::select(ID,ID_site,precip_c,temp_c) %>%
  rename(mean_temp=temp_c) %>% 
  rename(mean_precip=precip_c) #%>% drop_na()

saveRDS(meteo3,file = "meteo3.RDS")

# Import/clean IM and FOR DEP ---------------------------------------------

#Forests
dp_dem <- read_csv("raw_data/Forests_deposition/dp_dem.csv", 
                             col_types = cols(date_start = col_date(format = "%Y-%m-%d"), 
                                              date_end = col_date(format = "%Y-%m-%d"),
                                              sampler_id = col_double()))
FOR_dep <- dp_dem %>% dplyr::select(1:8, 14, 16, 19, 40, 41, 49)
rm(dp_dem)
#filter for throughfall deposition
FOR_dep <- FOR_dep %>% filter(code_sampler == 1)
FOR_dep$code_sampler <- NULL
FOR_dep <- distinct(FOR_dep)
FOR_dep_bak <- FOR_dep
#-1 marks missing data replace with NA #
FOR_dep <-
  FOR_dep %>% mutate(n_nh4 = replace(n_nh4, which(n_nh4 < 0), NA)) %>%
  mutate(n_no3 = replace(n_no3, which(n_no3 < 0), NA)) %>% 
  mutate(n_total = replace(n_total, which(n_total < 0), NA))


# DEP-Deposition dat ------------------------------------------------------
#make line ID
FOR_dep <-
  transform(
    FOR_dep,
    ID = paste0(
      FOR_dep$survey_year,
      sep = "_",
      FOR_dep$code_country,
      sep = "_",
      FOR_dep$code_plot
    )
  )

#make line ID_site
FOR_dep <-
  transform(
    FOR_dep,
    ID_site = paste0(
      FOR_dep$code_country,
      sep = "_",
      FOR_dep$code_plot
    )
  )
#add country
FOR_dep <- left_join(FOR_dep, FOR_country, by = "code_country")
FOR_dep$ID <- as.character(FOR_dep$ID)
FOR_dep$ID_site <- as.character(FOR_dep$ID_site)
FOR_dep$country <- as.character(FOR_dep$country)
#create ID including period number
FOR_dep <-
  transform(
    FOR_dep,
    IDp = paste0(
      FOR_dep$ID,
      sep = "_",
      FOR_dep$period_number
    )
  )

#no period time data for early period- seperate and use annual means
FOR_dep_t <- filter(FOR_dep, !is.na(date_start))
FOR_dep_nt <- filter(FOR_dep, is.na(date_start))
#sort by month
FOR_dep_t$monthly <- round_date(FOR_dep_t$date_start, unit = "month")
#keep just month of date
FOR_dep_t$month <- month(FOR_dep_t$monthly)
#create ID with month appended to match precip data
FOR_dep_t <-
  transform(
    FOR_dep_t,
    IDm = paste0(
      FOR_dep_t$ID,
      sep = "_",
      FOR_dep_t$month
    )
  )
#Fill missing months with value from month before in sequence
#FOR_dep_t <- FOR_dep_t %>% group_by(ID) %>% arrange(month) %>% tidyr::fill(n_nh4, n_no3, .direction="downup") %>% 
#  ungroup()

#average per plot/year for checks
FOR_dep_means <- FOR_dep %>% group_by(ID) %>%
  summarise_at(vars(n_nh4, n_no3), funs(mean), na.rm = TRUE) 

#TF conc per month x precipitation summed to year
FOR_dep_m <- FOR_dep_t %>% group_by(ID,ID_site,IDm) %>% 
  summarise_at(vars(n_nh4, n_no3), funs(mean), na.rm = TRUE) 
FOR_dep_mm <- left_join(FOR_dep_m, PR_m, by = "IDm")
FOR_dep_mm <- FOR_dep_mm %>% rename(ID_site=ID_site.x)
#Where concentrations are missing and precip is < 0.1, set deposition to zero
#otherwise data is lost to calculations with NaNs
FOR_dep_mm <- FOR_dep_mm %>% mutate(n_nh4 = replace(n_nh4, which(is.nan(n_nh4) & monthly_total<=0.1), 0))
FOR_dep_mm <- FOR_dep_mm %>% mutate(n_no3 = replace(n_no3, which(is.nan(n_no3) & monthly_total<=0.1), 0))

FOR_dep_mm <- FOR_dep_mm %>% rename(n_nh4_c = n_nh4, n_no3_c = n_no3) %>% 
  mutate(n_nh4 = (n_nh4_c*monthly_total)*0.01) %>% 
  mutate(n_no3 = (n_no3_c*monthly_total)*0.01) %>%
  group_by(ID.x,ID_site) %>% 
  summarise_at(vars(n_nh4, n_no3), list(sum), na.rm = TRUE) %>% 
  rename(ID=ID.x) %>% 
  #where plot precipitation is not available we get 0, replace with NA and so
  #these will be given as mean annual deposition*mean annual regional precipitation
  mutate(n_nh4 = replace(n_nh4, which(n_nh4 < 0.1), NA)) %>% 
  mutate(n_no3 = replace(n_no3, which(n_no3 < 0.1), NA))

#there are some with period data but no matching monthly precipitation data- move to the 
#set with no period data
noperiod <- filter(FOR_dep_mm, is.na(n_nh4))
noperiods <- noperiod$ID
noperiod.df <- filter(FOR_dep_m, ID %in% noperiods)
FOR_dep_mm <- filter(FOR_dep_mm, ID %!in% noperiods)

#average per plot/year where period time data are not available to match with monthly precip
FOR_dep_nt <- dplyr::select(FOR_dep_nt, ID, ID_site, n_nh4, n_no3)
FOR_dep_nt2 <- full_join(FOR_dep_nt, noperiod.df)
FOR_dep_nt2$IDm <- NULL
FOR_dep_nt_a <- FOR_dep_nt2 %>% group_by(ID,ID_site) %>%
  summarise_at(vars(n_nh4, n_no3), list(mean), na.rm = TRUE) 
FOR_dep_nt_a <- left_join(FOR_dep_nt_a, meteo_ann, by="ID")
FOR_dep_nt_a <- FOR_dep_nt_a %>% 
  mutate(n_nh4 = (n_nh4*precip)*0.01) %>% 
  mutate(n_no3 = (n_no3*precip)*0.01)
FOR_dep_nt_a <- FOR_dep_nt_a %>% rename(ID_site = ID_site.x) %>% dplyr::select(ID,ID_site,n_nh4,n_no3)
#combine
FOR_dep2 <- full_join(FOR_dep_nt_a, FOR_dep_mm) %>% drop_na()

#IM
IM_dep_long <- read_csv(
    "raw_data/IM deposition/ICPIM_TF_data.csv",
    col_types = cols(
    DataProcessingDate = col_skip(),
    DataQualityFlag = col_skip(),
    Day = col_skip(),
    DeterInfo = col_skip(),
    InstCode = col_skip(),
    PretreCode = col_skip(),
    PretreInfo = col_skip(),
    SamplingLevel = col_skip(),
    SpatialPool = col_skip(),
    ParameterList = col_skip(),
    DeterCode = col_skip(),
    ParameterInfo = col_skip(),
    NeedleAgeCode = col_skip(),
    StatusFlag = col_skip(),
    Unit = col_skip(),
    ListMedium = col_skip(),
    Subprog = col_skip(),
    YearMonth = col_date(format = "%Y%m")
  )
)

#add year column
IM_dep_long$survey_year <-  year(IM_dep_long$YearMonth)
#make line ID
IM_dep_long <-
  transform(IM_dep_long,
            ID = paste0(IM_dep_long$survey_year, sep = "_", IM_dep_long$AreaCode))
IM_dep_long$ID <- as.character(IM_dep_long$ID)

#make line ID_site
IM_dep_long <-
  transform(IM_dep_long,
            ID_site = paste0(IM_dep_long$AreaCode))
IM_dep_long$ID_site <- as.character(IM_dep_long$ID_site)

#ws_check function shows some columns need white space trimming out
IM_dep_long$StationCode <- trimws(IM_dep_long$StationCode, "both")
IM_dep_long$Value <- trimws(IM_dep_long$Value, "both")
IM_dep_long$Sample_ID <- trimws(IM_dep_long$Sample_ID, "both")
IM_dep_long$Value <- as.numeric(IM_dep_long$Value)

#add country
IM_dep_long <-
  mutate(IM_dep_long, country = str_sub(IM_dep_long$AreaCode, 1, 2)) #keep first two letters of country code and use to match
IM_dep_long$country <-
  IM_dep_long$country %>% str_replace_all(
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

#IM data is long, Forests is wide, spread IM to match
parameters <- c("NH4N", "NO3N", "PH", "PREC")
IM_dep <- filter(IM_dep_long, ParameterCode %in% parameters) %>%
  spread(., ParameterCode, Value)
IM_dep <- IM_dep %>% filter(StationCode != "281")
# rename columns to match FOR_dep
IM_dep <- IM_dep %>% rename(n_nh4 = NH4N, n_no3 = NO3N, ph = PH, precip = PREC)
#monthly precip x TF conc
IM_dep <- IM_dep %>% 
  mutate(n_nh4_m = (n_nh4*precip)*0.01) %>% 
  mutate(n_no3_m = (n_no3*precip)*0.01)
# average per site/year combo
IM_dep_means <- IM_dep %>% group_by(ID,ID_site,StationCode,survey_year) %>%
  summarise_at(vars(n_nh4, n_no3), funs(sum), na.rm = TRUE) 

IM_dep_means <- IM_dep_means %>% group_by(ID,ID_site,survey_year) %>%
  summarise_at(vars(n_nh4, n_no3), funs(mean), na.rm = TRUE) 

IM_dep2 <- drop_na(IM_dep_means)
IM_dep2$survey_year <- NULL

# Combine IM and FOR dep --------------------------------------------------
dep <- full_join(IM_dep2, FOR_dep2)
depbak <- dep

#remove some outliers- there are a very few plots with crazy high levels (an order of magnitude)
#probably
#near a point source (or error in recording). Remove those observations outside the
#upper limit with the value of 98th %ile

q <- quantile(dep$n_nh4, .99, na.rm = TRUE)
dep <- dep %>% filter(n_nh4 < q)

q <- quantile(dep$n_no3, .99, na.rm = TRUE)
dep <- dep %>% filter(n_no3 < q)

# VEG- Vegetation Data ----------------------------------------------------
# Import/clean FOR VEG ----------------------------------------------------

#veg data ICP Forests
FOR_veg <-   read_delim("raw_data/Veg_Forests/gv_vem.csv",
             ";", col_types = cols(change_date = col_skip(), 
                                   code_certainty = col_skip(),
                                   code_substrate = col_character(), 
                                   other_obs = col_skip(), q_flag = col_skip(), 
                                   sample_id = col_character()),
             trim_ws = TRUE)
FOR_veg_site_info <-
  read_delim("raw_data/Veg_Forests/gv_plv.csv",
             ";", col_types = cols(bare_soil_cover = col_double(), 
                                   change_date = col_skip(), code_survey_type = col_character(), 
                                   herb_layer_cover = col_double(), 
                                   litter_cover = col_double(), mosses_cover = col_double(), 
                                   no_members = col_skip(), q_flag = col_skip(), 
                                   sample_id = col_character(), shrub_layer_cover = col_double(), 
                                   team_id = col_skip(), tree_layer_cover = col_double()),
             trim_ws = TRUE)
#raw_veg_FOR <- read_delim("~/Documents/R/Paper2/raw_data/Veg_Forests/gv_vem.csv", ";", trim_ws = TRUE)
#Understorey species metadata
FOR_veg_species <-
  read_delim("raw_data/FOR_species.csv",
             ";", col_types = cols(cmnt = col_skip(), valid_to_survey_year = col_double()),
             trim_ws = TRUE) #understorey species codes

FOR_veg_species <-
  FOR_veg_species %>% dplyr::select(code, family, genus, species) %>%
  mutate(species_name = paste(genus, species, sep = "_")) %>%
  dplyr::select(code, family, species_name) %>%
  rename(code_species = code)

FOR_veg_species$species_name <-
  gsub(" agg.", "", FOR_veg_species$species_name)


FOR_veg$code_species <-
  gsub("[^0-9.]", "", FOR_veg$code_species) #strip letters from species codes to allow matching with species list
FOR_veg$code_species <-
  gsub('^\\.|\\.$', '', FOR_veg$code_species) #regex: the two backslashes, \\, in the regular expression escape the dot,
#., which would actually mean any character. The caret, ^, marks the beginning of the string, the dollar, $,
#the end of the string. The, |, is a logical "or". So in essence the regular expression matches a dot at the
#beginning of the string or a dot at the end of the string and replaces it with an empty string.

#use species/country codes to add species/country names
FOR_veg <- left_join(FOR_veg, FOR_veg_species, by = "code_species")
FOR_veg <- left_join(FOR_veg, FOR_country, by = "code_country")

#make new column ID
FOR_veg <- transform(
  FOR_veg,
  ID = paste0(
    FOR_veg$survey_year,
    sep = "_",
    FOR_veg$code_country,
    sep = "_",
    FOR_veg$code_plot
  )
)

#from 2011 new ID system introduced: "For each sampling unit (inside outside fence, CSA,
#etc) use an unique ID,
#which must not change over time. If you have already a unique number for each subplot
#and you have submitted 
#this until the monitoring year 2010 as ‘Survey number’ use this number as new Sample ID". Sample ID is not found
#in data prior to 2011, column is all NAs.
#For post 2010 plots use Sample ID as most granular level
#a <- FOR_veg %>% filter(survey_year>=2011)
#a$survey_number <- NULL
#a$ID <- NULL
#make new column ID
# a <- transform(
#   a,
#   ID = paste0(
#     a$survey_year,
#     sep = "_",
#     a$code_country,
#     sep = "_",
#     a$sample_id
#   )
# )
# aa <- a %>% group_by(ID, species_name, code_layer_surface) %>%
#   mutate(species_cover2 = mean(species_cover))
# aaa <- ungroup(aa)
# aaa$survey_number <- NULL
# aaa$species_cover <- NULL
# 
# aaa <- rename(aaa, survey_number = sample_id)
# aaa <- rename(aaa, species_cover = species_cover2)
# aaa$survey_number <- as.numeric(aaa$survey_number)
# 
# b <- FOR_veg %>% filter(survey_year<=2010)
# FOR_veg2 <- full_join(aaa, b)
# rm(a,aa,aaa,b)

#create site ID column in site info dataframe
FOR_veg_site_info <-
  transform(
    FOR_veg_site_info,
    ID = paste0(
      FOR_veg_site_info$survey_year,
      sep = "_",
      FOR_veg_site_info$code_country,
      sep = "_",
      FOR_veg_site_info$code_plot
    )
  )

For_veg_bak <- FOR_veg

#create new column for vascular/non vascular plants/Lichens (also added to IM data)
FOR_veg <-
  mutate(
    FOR_veg,
    vascular = if_else(
      FOR_veg$code_layer_surface == "4" &
        FOR_veg$family == "Lichenes" ,
      "lich",
      if_else(
        FOR_veg$code_layer_surface == "4" &
          FOR_veg$family != "Lichenes",
        "non_vasc",
        "vasc"
      )
    )
  )

#create plot ID_site without year
FOR_veg <- mutate(FOR_veg, ID_site = substr(FOR_veg$ID, 6, 12))
#add locations
FOR_veg <- left_join(FOR_veg, FOR_veg_locations)

#add fence, remove fenced plots (few present but not comparable)
fence <- dplyr::select(FOR_veg_site_info, code_fence, code_line)
testveg <- left_join(FOR_veg, fence, by="code_line")
testveg1 <- filter(testveg, code_fence != 1)
testveg1 <- dplyr::filter(testveg, !code_fence %in% 1)
testveg1$code_fence <- NULL
FOR_veg <- testveg1

#filter tree level data####
FOR_veg_trees <- filter(FOR_veg, code_layer_surface == "1")
FOR_veg_trees$species_name <-
  gsub("_", " ", FOR_veg_trees$species_name)
FOR_veg_trees <-
  filter(FOR_veg_trees, species_name != "Hedera helix") #interested in tree species, remove non
#tree species in canopy data (ivy).

decid_FOR_veg_trees <-
  c(
    "broadleaves",
    "Fagus sylvatica" ,
    "Carpinus betulus",
    "Quercus ilex" ,
    "Betula pubescens",
    "Betula pendula",
    "Quercus petraea" ,
    "Quercus robur",
    "Quercus cerris",
    "Castanea sativa",
    "Fraxinus excelsior",
    "Sorbus torminalis",
    "Populus tremula",
    "Tilia cordata",
    "Alnus glutinosa",
    "Prunus avium",
    "Quercus sp.",
    "Frangula alnus",
    "Acer campestre",
    "Acer pseudoplatanus",
    "Sorbus aucuparia",
    "Quercus pubescens",
    "Acer platanoides" ,
    "Populus nigra",
    "Sorbus aria",
    "Salix caprea",
    "Corylus avellana",
    "Fraxinus angustifolia",
    "Quercus rubra",
    "Quercus pyrenaica",
    "Quercus suber",
    "Quercus frainetto",
    "Salix cinerea",
    "Quercus pedunculiflora",
    "Tilia platyphyllos",
    "Sambucus nigra",
    "Populus x canescens",
    "Alnus incana",
    "Quercus x andegavensis",
    "Tilia sp.",
    "Acer sp.",
    "Fraxinus ornus",
    "Fagus sylvatica" ,
    "Carpinus betulus",
    "Quercus ilex" ,
    "Betula pubescens",
    "Quercus rotundifolia",
    "Betula pendula",
    "Quercus petraea" ,
    "Quercus faginea",
    "Quercus robur",
    "Quercus cerris",
    "Castanea sativa",
    "Fraxinus excelsior",
    "Quercus petraeaorrobur",
    "Sorbus torminalis",
    "Populus tremula",
    "Tilia cordata",
    "Alnus glutinosa",
    "Prunus avium",
    "Ulmus glabra",
    "Quercus sp.",
    "Frangula alnus",
    "Acer campestre",
    "Acer pseudoplatanus",
    "Sorbus aucuparia",
    "Acer platanoides" ,
    "Populus nigra",
    "Sorbus aria",
    "Salix caprea",
    "Corylus avellana",
    "Quercus rubra",
    "Eucalyptus sp",
    "Populus canescens" ,
    "Robinia pseudoacacia",
    "Tilia platyphyllos",
    "Sambucus nigra",
    "Populus x canescens"
  )
conif_FOR_veg_trees <-
  c(
    "conifers",
    "Abies alba" ,
    "Larix decidua" ,
    "Picea abies",
    "Pinus nigra",
    "Pinus sylvestris",
    "Pinus pinaster",
    "Pinus canariensis",
    "Pinus pinea",
    "Pinus halepensis",
    "Pinus cembra",
    "Pinus mugo",
    "Robinia pseudacacia",
    "Ulmus laevis",
    "Pinus halepensis",
    "Pinus radiata",
    "Picea sitchensis",
    "Pseudotsuga menziesii",
    "Abies grandis",
    "Pinus uncinata",
    "Pinus contorta",
    "Abies alba" ,
    "Abies borisiiregis" ,
    "Larix decidua" ,
    "Picea abies",
    "Pinus nigra",
    "Pinus sylvestris",
    "Pinus pinaster",
    "Juniperus oxycedrus",
    "Erica arborea",
    "Juniperus thurifera",
    "Pseudotsuga menziesii",
    "Abies grandis"
  )

FOR_veg_tree_groups <-
  FOR_veg_trees %>% dplyr::select(ID, species_cover, species_name) %>%
  group_by(ID) %>% top_n(1, species_cover) #take top cover tree species in each plot

#duplicated when removing eg line number
FOR_veg_tree_groups <-
  FOR_veg_tree_groups[!duplicated(FOR_veg_tree_groups$ID),]


FOR_veg_tree_groups2 <-
  FOR_veg_tree_groups %>% mutate(grp_tree_species = ifelse(
    species_name %in% decid_FOR_veg_trees,
    'broadleaves',
    ifelse(species_name %in% conif_FOR_veg_trees, 'conifers',
           'other')
  ))
FOR_veg_tree_groups <-
  filter(FOR_veg_tree_groups2, grp_tree_species != "other")
FOR_veg_tree_groups$species_cover <- NULL
#FOR_veg_tree_groups <- rename(FOR_veg_tree_groups, species_name2 = species_name)

#filter shrub level data####
FOR_veg_shrubs <-
  filter(FOR_veg, code_layer_surface  %in% c("2", "5", "6"))
FOR_veg_shrubs$species_name <-
  gsub("_", " ", FOR_veg_shrubs$species_name)
FOR_veg_shrubs <-
  filter(FOR_veg_shrubs, species_name != "Hedera helix") #interested in tree species, remove non
#tree species in canopy data (ivy).

decid_FOR_veg_shrubs <-
  c(
    "broadleaves",
    "Fagus sylvatica" ,
    "Carpinus betulus",
    "Quercus ilex" ,
    "Betula pubescens",
    "Betula pendula",
    "Quercus petraea" ,
    "Quercus robur",
    "Quercus cerris",
    "Castanea sativa",
    "Fraxinus excelsior",
    "Sorbus torminalis",
    "Populus tremula",
    "Tilia cordata",
    "Alnus glutinosa",
    "Prunus avium",
    "Quercus sp.",
    "Frangula alnus",
    "Acer campestre",
    "Acer pseudoplatanus",
    "Sorbus aucuparia",
    "Quercus pubescens",
    "Acer platanoides" ,
    "Populus nigra",
    "Sorbus aria",
    "Salix caprea",
    "Corylus avellana",
    "Fraxinus angustifolia",
    "Quercus rubra",
    "Quercus pyrenaica",
    "Quercus suber",
    "Quercus frainetto",
    "Salix cinerea",
    "Quercus pedunculiflora",
    "Tilia platyphyllos",
    "Sambucus nigra",
    "Populus x canescens",
    "Alnus incana",
    "Quercus x andegavensis",
    "Tilia sp.",
    "Acer sp.",
    "Fraxinus ornus",
    "Fagus sylvatica" ,
    "Carpinus betulus",
    "Quercus ilex" ,
    "Betula pubescens",
    "Quercus rotundifolia",
    "Betula pendula",
    "Quercus petraea" ,
    "Quercus faginea",
    "Quercus robur",
    "Quercus cerris",
    "Castanea sativa",
    "Fraxinus excelsior",
    "Quercus petraeaorrobur",
    "Sorbus torminalis",
    "Populus tremula",
    "Tilia cordata",
    "Alnus glutinosa",
    "Prunus avium",
    "Ulmus glabra",
    "Quercus sp.",
    "Frangula alnus",
    "Acer campestre",
    "Acer pseudoplatanus",
    "Sorbus aucuparia",
    "Acer platanoides" ,
    "Populus nigra",
    "Sorbus aria",
    "Salix caprea",
    "Corylus avellana",
    "Quercus rubra",
    "Eucalyptus sp",
    "Populus canescens" ,
    "Robinia pseudoacacia",
    "Tilia platyphyllos",
    "Sambucus nigra",
    "Populus x canescens"
  )
conif_FOR_veg_shrubs <-
  c(
    "conifers",
    "Abies alba" ,
    "Larix decidua" ,
    "Picea abies",
    "Pinus nigra",
    "Pinus sylvestris",
    "Pinus pinaster",
    "Pinus canariensis",
    "Pinus pinea",
    "Pinus halepensis",
    "Pinus cembra",
    "Pinus mugo",
    "Robinia pseudacacia",
    "Ulmus laevis",
    "Pinus halepensis",
    "Pinus radiata",
    "Picea sitchensis",
    "Pseudotsuga menziesii",
    "Abies grandis",
    "Pinus uncinata",
    "Pinus contorta",
    "Abies alba" ,
    "Abies borisiiregis" ,
    "Larix decidua" ,
    "Picea abies",
    "Pinus nigra",
    "Pinus sylvestris",
    "Pinus pinaster",
    "Juniperus oxycedrus",
    "Erica arborea",
    "Juniperus thurifera",
    "Pseudotsuga menziesii",
    "Abies grandis"
  )

FOR_veg_shrub_groups <-
  FOR_veg_shrubs %>% dplyr::select(ID, species_cover, species_name) %>%
  group_by(ID) %>% top_n(1, species_cover) #take top cover species in each plot


FOR_veg_shrub_groups2 <-
  FOR_veg_shrub_groups %>% mutate(grp_tree_species = ifelse(
    species_name %in% decid_FOR_veg_trees,
    'broadleaves',
    ifelse(species_name %in% conif_FOR_veg_trees, 'conifers',
           'other')
  ))
FOR_veg_shrub_groups <-
  filter(FOR_veg_shrub_groups2, grp_tree_species != "other")
FOR_veg_shrub_groups$species_cover <- NULL

#ctotal over per plot
tem_sh <- FOR_veg_shrubs
tem_sh <-
  transform(
    tem_sh,
    ID_site = paste0(
      tem_sh$code_country,
      sep = "_",
      tem_sh$code_plot
    )
  )

tem_sh <-
  tem_sh %>% group_by(survey_year, country, ID, ID_site,species_name) %>%
  summarise(mean_species = mean(species_cover))
#add species covers to get per plot canopy total (can sum to >100% due to overlap)
tem_sh <- tem_sh %>% group_by(survey_year, country, ID, ID_site) %>%
  summarise(sum_canopy = sum(mean_species))
FOR_shrub_site <- tem_sh
#per ID mean of plot cover
FOR_shrub_means <- FOR_shrub_site %>% group_by(survey_year, country, ID,ID_site) %>%
  summarise(mean_canopy = mean(sum_canopy))

#keep only moss and herb layers to allow combination with IM_veg
FOR_veg <- filter(FOR_veg, code_layer_surface %in% c("3", "4"))

FOR_veg <- FOR_veg %>%  dplyr::select(-c(survey_number,sample_id, line_nr, code_line)) %>% 
  group_by(ID, species_name) %>% 
  mutate(cover = mean(species_cover)) %>% dplyr::select(-species_cover) %>% 
  distinct()
# Import/clean IM VEG VS --------------------------------------------------

#(this file includes the extra data from Thomas (much bigger i.e. more recent)) but only up to 2012
VS_Assessment_2012 <-
  read_csv(
    "raw_data/Veg_IM/data_from_Thomas/VS_Assessment_2012_AT_DE_ES_SE.csv",
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
IM_veg <- VS_Assessment_2012
as.data.frame(table(IM_veg$SPECIES))
length(complete.cases(IM_veg$SPECIES) == TRUE)

#Don't need information on month just year. Add column for year only (and keep month just in case...
IM_veg$survey_year <- as.integer(substr(IM_veg$YYYYMM, 0, 4))
IM_veg$LISTMED <- NULL

#add country
IM_veg <-
  mutate(IM_veg, country = str_sub(IM_veg$AREA, 1, 2)) #keep first two letters of country code and use to match
IM_veg$country <-
  IM_veg$country %>% str_replace_all(
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
#total covers have NA for species, separate out these
IM_veg_totals <- filter(IM_veg,!is.na(SPECIES))
IM_veg <- filter(IM_veg,!is.na(SPECIES))

#harmonise labelling with Forests data

#species list provided (AASpecies_List_VS_VG_FINAL.xlsx) does NOT include bryophytes...
temp <- filter(IM_veg, PARAM == "COVE_B")
IM_bryo_spp_list <- unique(temp$SPECIES)
IM_bryo_spp_df <- as.data.frame(cbind(IM_bryo_spp_list, c(1:209)))
names(IM_bryo_spp_df) <- c("species_name", "number")

#fuzzy matching was not a success. Export to csv and do manually
#write.csv(IM_bryo_spp_df, file = "bryophyte_names_IM.csv")
#and import names file
bryophyte_names_IM <- read_csv("raw_data/bryophyte_names_IM.csv",
                               col_types = cols(X1 = col_skip()))

#write.csv(bryophyte_names_IM, file = "corrected_bryophyte_names_IM.csv")

#added bryophyte names to species list to make combined csv- import this
VG_VS_spp_list <-
  read_csv("raw_data/Veg_IM/data_from_Thomas/VG_spp_list_bryo_added.csv",
           trim_ws = TRUE)
VG_VS_spp_list$Genus <- NULL
VG_VS_spp_list$sp <- NULL

#match on full names in species list
new <- IM_veg
new[] <-
  VG_VS_spp_list$species_name[match(unlist(IM_veg$SPECIES), VG_VS_spp_list$SPECIES)]

#copy to veg dataframe
IM_veg <- mutate(IM_veg, species_name = new$SPECIES)
#cut to genus and species
IM_veg$species_name  <- word(IM_veg$species_name, 1, 2)
#replace spaces with underscore to match Forests
IM_veg$species_name <-  gsub("\\ ", "_", IM_veg$species_name)
#make new column ID (to year/site level)
IM_veg <-
  transform(IM_veg, ID = paste0(IM_veg$survey_year, sep = "_",  IM_veg$AREA))

#filter table for only parameters of interest to reduce size and keep canopy cover seperate
IM_veg_cover <-
  IM_veg %>% filter(IM_veg$PARAM %in% c("COVE_T", "COVE_T1", "COVE_T2", "COVE_T3"))
IM_veg <- IM_veg %>% filter(IM_veg$PARAM %in% c("COVE_F", "COVE_B"))

IM_veg_bak <- IM_veg
#create family column based on Forest veg species list
new <- IM_veg
new[] <-
  FOR_veg_species$family[match(unlist(IM_veg), FOR_veg_species$species_name)]
#copy to veg dataframe
IM_veg <- mutate(IM_veg, family = new$species_name)
#create new column for vascular/non vascular plants (also added to FOR data)
IM_veg <-
  mutate(IM_veg, vascular = if_else(IM_veg$PARAM == "COVE_F", "vasc", "non_vasc"))

#create new column for vascular/non vascular plants/Lichens (also added to forests data)
IM_veg <-
  mutate(
    IM_veg,
    vascular = if_else(
      IM_veg$PARAM == "COVE_B" &
        IM_veg$family == "Lichenes" ,
      "lich",
      if_else(
        IM_veg$PARAM == "COVE_B" &
          IM_veg$family != "Lichenes",
        "non_vasc",
        "vasc"
      )
    )
  )


#find rows with missing species names
nas <- IM_veg[is.na(IM_veg$species_name), ]
table(nas$species_name)
table(nas$country)

#add location
IM_locations <- IM_locations %>% mutate(AREA = ID_site)
IM_veg <- left_join(IM_veg, IM_locations, by = "AREA")

#tree canopy data####
#non tree species in canopy data. Mostly ivy.
IM_veg_cover$species_name <-
  gsub("_", " ", IM_veg_cover$species_name)

IM_veg_tree_groups <-
  IM_veg_cover %>% dplyr::select(ID, VALUE, species_name) %>%
  group_by(ID) %>% top_n(1, VALUE) #take top cover species in each plot

#some ties resulting in double entries
IM_veg_tree_groups <-
  IM_veg_tree_groups[!duplicated(IM_veg_tree_groups$ID),]

IM_veg_tree_groups <- IM_veg_tree_groups %>%
  mutate(grp_tree_species = ifelse(
    species_name %in% decid_FOR_veg_trees,
    'broadleaves',
    ifelse(species_name %in% conif_FOR_veg_trees, 'conifers', 'other')
  ))

IM_veg_tree_groups$VALUE <- NULL

IM_veg <- IM_veg %>% dplyr::select(-SCODE) %>% group_by(ID, species_name) %>%
  mutate(VALUE = mean(VALUE)) %>%  distinct()

# Combine IM and FOR VEG --------------------------------------------------
veg <- full_join(FOR_veg, IM_veg)
#species cover now in two columns, make one
veg <- mutate(veg, cover = coalesce(cover, VALUE))
veg$species_cover <- NULL
veg$VALUE <- NULL
veg$other_obs <- NULL
veg$q_flag <- NULL
veg$code_certainty <- NULL
#one entry is 500% cover for single species, error, remove
veg <- filter(veg, cover <= 100)
#check for NAs and deal with if neccessary
veg <- filter(veg,!is.na(cover))
veg <- filter(veg,!is.na(species_name))

#Problem- some plots are dated 1961. Obviously wrong (but what is real year?)
#Remove for now...
#unique(veg$survey_year)
veg_bak <- veg
veg <- filter(veg, survey_year != "1961")

#species which are ID'd only to genus cannot get Ellenbergs, drop them
remove.list <- paste(c("SP.", "sp.", "agg."), collapse = '|')
veg <- veg %>% filter(!grepl(remove.list, species_name))

#remove underscores until after ellenberg lookups with TR8/taxize
veg$species_name <-  gsub("_", " ", veg$species_name)

#modify names that have other possibilities to those versions found in ellenberg lists
#(from worth_finding)
veg$species_name <-
  gsub("Hieracium juranum", "Hieracium jurassicum", veg$species_name)
veg$species_name <-
  gsub("CAREX ALBA", "Carex alba", veg$species_name)
veg$species_name <-
  gsub("CARDAMINE TRIFOLIA", "Cardamine trifolia", veg$species_name)
veg$species_name <-
  gsub("Adenostyles alpina", "Adenostyles glabra", veg$species_name)
veg$species_name <-
  gsub("Acer obtusatum", "Acer opalus", veg$species_name)
veg$species_name <-
  gsub("Majanthemum bifolium",
       "Maianthemum bifolium",
       veg$species_name)
veg$species_name <-
  gsub("Stachys officinalis",
       "Betonica officinalis",
       veg$species_name)
veg$species_name <-
  gsub("Senecio ovatus", "Senecio fuchsii", veg$species_name)
veg$species_name <-
  gsub("Phegopteris connectilis",
       "Thelypteris phegopteris",
       veg$species_name)
veg$species_name <-
  gsub("Viola sylvestris", "Viola reichenbachiana", veg$species_name)
veg$species_name <-
  gsub("Ceratocapnos claviculata",
       "Corydalis claviculata",
       veg$species_name)
veg$species_name <-
  gsub("Orthodicranum montanum",
       "Dicranum montanum",
       veg$species_name)
veg$species_name <-
  gsub("Populus x canescens", "Populus alba", veg$species_name)
veg$species_name <-
  gsub("Carex curta", "Carex canescens", veg$species_name)
veg$species_name <-
  gsub("Calamagrostis epigeios",
       "Calamagrostis epigejos",
       veg$species_name)
veg$species_name <-
  gsub("ASTER BELLIDIASTRUM", "Aster bellidiastrum", veg$species_name)
#veg$species_name <- gsub("Brachythecium oedipodium", "Brachythecium starkei", veg$species_name)
veg$species_name <-
  gsub("Eurhynchium angustirete",
       "Eurhynchium striatulum",
       veg$species_name)
veg$species_name <-
  gsub("BUPHTHALMUM SALICIFOLIUM",
       "Buphthalmum salicifolium",
       veg$species_name)
veg$species_name <-
  gsub("Sesleria albicans", "Sesleria caerulea", veg$species_name)


#There are entries with NA NA for species, listed with explanation
#"provisional code during BioSoil BioDiv Evaluation".Discard...
veg <- filter(veg, species_name != "NA NA")

#few duplicate rows, remove
veg <- distinct(veg)

#remove entries with tiny cover values
veg <- filter(veg, cover>= 0.01)

#subset veg data to vascular/nonvascular
veg_bryo <- dplyr::filter(veg, vascular == "non_vasc")
veg_vasc <- dplyr::filter(veg, vascular == "vasc")
veg_lich <- dplyr::filter(veg, vascular == "lich")
#clean up, some vascular are in bryo and vice versa e.g tree seedlings recorded amongst moss layer
intersect(veg_bryo$species_name, veg_vasc$species_name)
veg_bryo <-
  filter(veg_bryo,!species_name %in% FOR_tree_species$species_name)
intersect(veg_bryo$species_name, veg_vasc$species_name)
vasc_remove_from_bryo <-
  c(
    "Frangula alnus",
    "Oxalis acetosella",
    "Crocus vernus",
    "Maianthemum bifolium",
    "Hedera helix",
    "Calluna vulgaris",
    "Rubus caesius",
    "Acacia dealbata",
    "Digitalis purpurea",
    "Vaccinium vitis-idaea",
    "Vaccinium oxycoccos",
    "Viburnum opulus",
    "Sambucus racemosa",
    "Sambucus ebulus",
    "Cirsium palustre",
    "Polygonum longipes",
    "Polystichum lonchitis",
    "Pteridium aquilinum",
    "Ptilium crista-castrensis",
    "Silene italica",
    "Poa trivialis",
    "Hieracium prenanthoides",
    "Ruscus aculeatus",
    "Corydalis claviculata",
    "Symphytum tuberosum",
    "Cytisus scoparius",
    "Lycopodium annotinum",
    "Cardamine chelidonia",
    "Lonicera periclymenum"
  )

move_to_vasc <- filter(veg_bryo,species_name %in% vasc_remove_from_bryo)
veg_bryo <- filter(veg_bryo,!species_name %in% vasc_remove_from_bryo)

bryo_remove_from_vasc <- intersect(veg_bryo$species_name, veg_vasc$species_name)
move_to_bryo <- filter(veg_vasc,species_name %in% bryo_remove_from_vasc)
veg_vasc <- filter(veg_vasc,!species_name %in% bryo_remove_from_vasc)

veg_vasc <- full_join(veg_vasc,move_to_vasc)
veg_bryo <- full_join(veg_bryo, move_to_bryo)

veg_vasc$vascular <- "vasc"
veg_bryo$vascular <- "non_vasc"

#possible duplicate rows, remove if any (none present)
veg_bryo <- distinct(veg_bryo)
veg_vasc <- distinct(veg_vasc)
veg_lich <- distinct(veg_lich)

#add species richness per site
veg_bryo %<>% group_by(ID) %>% mutate(richness = length(unique(species_name)))
veg_vasc %<>% group_by(ID) %>% mutate(richness = length(unique(species_name)))
veg_lich %<>% group_by(ID) %>% mutate(richness = length(unique(species_name)))

#add species richness per plot for combined dataframe (must come after subset richness above)
veg %<>% group_by(ID) %>% mutate(richness = length(unique(species_name)))
veg$survey_number <- NULL

#keep full versions, then..
veg_bryo_full <- veg_bryo
veg_vasc_full <- veg_vasc
veg_lich_full <- veg_lich

#join to make veg dataframe
veg <- full_join(veg_bryo_full, veg_vasc_full)
veg$survey_number <- NULL

veg_inc_lich <- full_join(veg, veg_lich)

# Ellenberg data ----------------------------------------------------------
#Compile a list of Ellenberg values for relevant species to use in later CWM analyses
#import BRYOATT csv of ellenberg values
# bryoatt_simple <-
#   read.csv("raw_data/bryoatt_simple.csv")
# bryophyte_ells <-  bryoatt_simple
# bryophyte_ells$Taxon.name <-
#   word(bryophyte_ells$Taxon.name, 1, 2, sep = " ") #cut out extra bits of taxon names (e.g sub-spp.)
# bryophyte_ells$Name_new <-
#   word(bryophyte_ells$Name_new, 1, 2, sep = " ")
# bryophyte_ells$Taxon.name <-
#   gsub(" ", "_", bryophyte_ells$Taxon.name)#change space to underscore to match FOR_veg
# bryophyte_ells$Name_new <-  gsub(" ", "_", bryophyte_ells$Name_new)
# bryophyte_ells <-
#   bryophyte_ells[!duplicated(bryophyte_ells$Taxon.name), ] #remove duplicates (sub-species entries)
# 
# #check whether names used in Forests/IM bryophyte data are "old" or "new" on BRYOATT
# #which new names are actually in the species list?
# temp.a <-
#   setdiff(bryophyte_ells$Name_new, bryophyte_ells$Taxon.name)
# intersect(temp.a, veg_bryo$species_name) #"Polytrichum_alpinum"
# #replace entry in Taxon.name with new name. Now bryophye_ells$Taxon.names matches taxonomy in veg_bryo
# bryophyte_ells$Taxon.name <-
#   gsub("Polytrichum_alpinum",
#        "Polytrichastrum_alpinum",
#        bryophyte_ells$Taxon.name)
# bryophyte_ells_bak <- bryophyte_ells
# 
# #Better to use BryForTrait and add extras from bryoatt
# #takes N values from BRYOATT but adds other functional traits
# BryForTrait <- read_csv(
#   "raw_data/BryForTrait.csv",
#   col_types = cols(
#     F = col_double(),
#     K = col_double(),
#     L = col_double(),
#     R = col_double(),
#     T = col_double(),
#     seta_length = col_double()
#   )
# )
# bryophyte_ells <- BryForTrait
# bryophyte_ells$species_name <-
#   word(bryophyte_ells$species_name, 1, 2, sep = " ") #cut out extra bits of taxon names (e.g sub-spp.)
# #bryophyte_ells$Taxon.name <-  gsub(" ", "_", bryophyte_ells$Taxon.name)#change space to underscore to match FOR_veg
# bryophyte_ells <-
#   bryophyte_ells[!duplicated(bryophyte_ells$species_name), ] #remove duplicates (sub-species entries)
# 
# #compare species lists for BRYOATT and BryForTrait
# bryophyte_ells_bak$Taxon.name <-
#   gsub("_", " ", bryophyte_ells_bak$Taxon.name)
# 
# BFT_names <- unique(bryophyte_ells$species_name)
# spp_list <- unique(veg$species_name)
# BRYOATT_names <- unique(bryophyte_ells_bak$Taxon.name)
# intersect(BFT_names, BRYOATT_names)
# extra_BRYOATT <- filter(bryophyte_ells_bak, Taxon.name %in%
#                           (setdiff(BRYOATT_names, BFT_names) %>% intersect(., spp_list))) #21
# #22 species in BRYOATT but not BryForTrait that occur in Veg spp list. Take these and add to BFT
# extra_BRYOATT$Name_new <- NULL
# extra_BRYOATT$S <- NULL
# extra_BRYOATT$HM <- NULL
# extra_BRYOATT <- rename(extra_BRYOATT, species_name = Taxon.name)
# 
# bryophyte_func <- bryophyte_ells
# bryophyte_ells <-
#   full_join(bryophyte_ells, extra_BRYOATT) %>% dplyr::select(species_name, L, T, K, F, R, N)
# 
# 
# #Take data from Ellenberg (Zeigerwerte)- No N data for bryophytes...
# #Use for vasc and lichens but combine with Bryoatt...
# Zeigerwerte_vasc_lich <-
#   read_csv(
#     "raw_data/Zeigerwerte_vasc_lich.csv",
#     trim_ws = TRUE,
#     col_types = cols(
#       F = col_double(),
#       K = col_double(),
#       N = col_double(),
#       R = col_double(),
#       T = col_double()
#     )
#   )
# ellenbergs_vasc_lich <- Zeigerwerte_vasc_lich
# ellenbergs_vasc_lich <-
#   rename(ellenbergs_vasc_lich, species_name = Name)
# ellenbergs_vasc_lich$species_name <-
#   word(ellenbergs_vasc_lich$species_name, 1, 2, sep = " ") #cut out extra bits of taxon names (e.g authority.)
# #ellenbergs_vasc_lich$Taxon.name <-  gsub(" ", "_", ellenbergs_vasc_lich$Taxon.name)#change space to underscore to match FOR_veg
# ellenbergs_vasc_lich <-
#   ellenbergs_vasc_lich[!duplicated(ellenbergs_vasc_lich$species_name), ] #remove duplicates
# ellenbergs_vasc_lich$F <-  gsub("=", "", ellenbergs_vasc_lich$F)
# ellenbergs_vasc_lich$F <-  gsub("~", "", ellenbergs_vasc_lich$F)
# 
# ellenbergs_vasc_lich$L <- as.integer(ellenbergs_vasc_lich$L)
# ellenbergs_vasc_lich$F <- as.integer(ellenbergs_vasc_lich$F)
# ellenbergs_vasc_lich$R <- as.integer(ellenbergs_vasc_lich$R)
# ellenbergs_vasc_lich$N <- as.integer(ellenbergs_vasc_lich$N)
# ellenbergs_vasc_lich$T <- as.integer(ellenbergs_vasc_lich$T)
# ellenbergs_vasc_lich$K <- as.integer(ellenbergs_vasc_lich$K)

#join to make one dataframe for all Ellenberg values
#ellenbergs <- full_join(bryophyte_ells, ellenbergs_vasc_lich)

#remove underscores until after ellenberg lookups with TR8/taxize
#ellenbergs$species_name <-  gsub("_", " ", ellenbergs$species_name)


# taxonomy check ----------------------------------------------------------
# 
# library(taxize)
# spp_list <- unique(veg$species_name)
# #tax_veg <- gnr_resolve(spp_list, canonical = TRUE, best_match_only = TRUE)
# ellenbergs_spp_list <- unique(ellenbergs$species_name)
# #tax_ellenbergs <- gnr_resolve(ellenbergs_spp_list, canonical = TRUE, best_match_only = TRUE)
# 
# #copy back check lists?
# #match on full names in species list
# # veg$species_name <-  gsub("_", " ", veg$species_name)
# # new <- veg
# # new[] <- tax_veg$matched_name2[match(unlist(veg), tax_veg$matched_name2)]
# 
# 
# #some missing values, check with Meditterranean databases for non-central europe spp.
# missing <-
#   setdiff(unique(veg$species_name), unique(ellenbergs$species_name))
# ## see the first lines of available_tr8 database
# library(TR8)
# available_tr8
# ##make vector of traits needed. Need more than one source to cover spp, will combine later
# #Database1
# my_ittraits_to_dload <-
#   c("ell_L_it", "ell_U_it", "ell_R_it", "ell_N_it", "ell_T_it")
# #Database2
# my_frtraits_to_dload <-
#   c("ell_L_fr", "ell_U_fr", "ell_R_fr", "ell_N_fr", "elle_T_fr")
# # #non ellenberg traits
# # non_ell_to_dload <- c("li_form_B", "canopy_height", "dispersal", "dispersal_morphology",
# #                       "seed_mass", "reprod_B", "strategy")
# 
# 
# ## now run tr8 and store the results in the my_traits object
# my_ittraits <-
#   tr8(species_list = missing, download_list = my_ittraits_to_dload)
# my_frtraits <-
#   tr8(species_list = missing, download_list = my_frtraits_to_dload)
# #extract to dataframe to allow combining
# it_ell = extract_traits(my_ittraits)
# fr_ell = extract_traits(my_frtraits)
# 
# #combine results to improve coverage and reduce NAs
# 
# 
# #first change column names to match
# it_ell <-
#   dplyr::rename(
#     it_ell,
#     L = ell_L_it,
#     F = ell_U_it,
#     R = ell_R_it,
#     N = ell_N_it,
#     T = ell_T_it
#   )
# fr_ell <-
#   dplyr::rename(
#     fr_ell,
#     L = ell_L_fr,
#     F = ell_U_fr,
#     R = ell_R_fr,
#     N = ell_N_fr,
#     T = elle_T_fr
#   )
# 
# # add the rownames as a proper column
# it_ell <- cbind(species_name = rownames(it_ell), it_ell)
# fr_ell <- cbind(species_name = rownames(fr_ell), fr_ell)
# 
# #ellenberg values only
# missing_ell <- left_join(it_ell, fr_ell, by = "species_name")
# missing_ell[1:11] <- lapply(missing_ell[1:11], as.character)
# missing_ell[2:11] <- lapply(missing_ell[2:11], as.double)
# missing_ell <-
#   missing_ell %>% mutate(L = coalesce(L.x, L.y)) %>% dplyr::select(-L.x,-L.y)
# missing_ell <-
#   missing_ell %>% mutate(N = coalesce(N.x, N.y)) %>% dplyr::select(-N.x,-N.y)
# missing_ell <-
#   missing_ell %>% mutate(R = coalesce(R.x, R.y)) %>% dplyr::select(-R.x,-R.y)
# missing_ell <-
#   missing_ell %>% mutate(F = coalesce(F.x, F.y)) %>% dplyr::select(-F.x,-F.y)
# missing_ell <-
#   missing_ell %>% mutate(T = coalesce(T.x, T.y)) %>% dplyr::select(-T.x,-T.y)
# 
# #keep only rows that aren't all NA, as in subset to all those that have at least two non-NA values
# #(two since the species name is always one non-NA value)
# missing_ell <-
#   missing_ell[rowSums(is.na(missing_ell)) < (length(missing_ell) - 1), ]
# 
# #add the "missing" data###
# ellenbergs <- full_join(ellenbergs, missing_ell)
# 
# #still missing
# missing <-
#   setdiff(unique(veg$species_name), unique(ellenbergs$species_name))
# #which ones occur often? can values be found for these from other sources?
# temp <- filter(veg, species_name %in% missing)
# worth_finding <-
#   temp %>% group_by(species_name) %>% tally() %>% filter(., n >= 50)

#saveRDS(ellenbergs,"ellenbergs.RDS")

ellenbergs <-  readRDS("output_data/ellenbergs.RDS")

# FOR CANOPY -------------------------------------------------
# total canopy calc -------------------------------------------------
FOR_veg_canopy <- FOR_veg_site_info %>% 
  dplyr::select(ID, tree_layer_cover) 

FOR_veg_canopy$tree_layer_cover <- as.numeric(as.character(FOR_veg_canopy$tree_layer_cover))
#3 x-99s recorded
FOR_veg_canopy$tree_layer_cover <- ifelse(FOR_veg_canopy$tree_layer_cover < 0, 99, FOR_veg_canopy$tree_layer_cover)
FOR_veg_canopy$tree_layer_cover <- ifelse(FOR_veg_canopy$tree_layer_cover == 0, NA, FOR_veg_canopy$tree_layer_cover)

saveRDS(FOR_veg_canopy,"FOR_veg_canopy.RDS")

tem <-
  FOR_veg_trees %>% group_by(survey_year, country, ID, species_name) %>%
  summarise(mean_species = mean(species_cover)) %>% ungroup()
#add species covers to get per plot canopy total
tem <- tem %>% group_by(survey_year, country, ID) %>%
  summarise(sum_canopy = sum(mean_species)) %>% ungroup()
FOR_CAN_site <- tem
#per ID mean of plot cover
FOR_CAN2_means <- FOR_CAN_site %>% group_by(survey_year, country, ID) %>%
  summarise(mean_canopy = mean(sum_canopy)) %>% ungroup()
#take out measurements below 1%
FOR_CAN2_means$mean_canopy <- replace(FOR_CAN2_means$mean_canopy,
                                      FOR_CAN2_means$mean_canopy < 1, NA)
t <- full_join(FOR_CAN2_means, FOR_veg_canopy, by = "ID")
#covert shading above 100% to 100% (apparent diff in countries? some use 100 as max, some not)
t$mean_canopy <- ifelse(t$mean_canopy > 100, 100, t$mean_canopy)
t$tree_layer_cover <- ifelse(t$tree_layer_cover > 100, 100, t$tree_layer_cover)

#prefer canopy cover data from level 2 veg plots summary data where available, but sometimes not there
#and never for IM plots, so use mean_canopy in those cases. OR reverse?
temp <- t %>%  mutate(canopy = coalesce(mean_canopy, tree_layer_cover))
FOR_canopy_means <- temp %>% group_by(ID) %>% summarise_at("canopy", mean)

#with shrubs
#In many cases shrub layer is available, this provides additional shading
#for bryophytes
#FOR_shrub_means <- ungroup(FOR_shrub_means)
tc <- full_join(FOR_canopy_means, dplyr::select(FOR_shrub_means, ID, ID_site, mean_canopy), by = "ID")
tc$mean_canopy <- ifelse(tc$mean_canopy > 100, 100, tc$mean_canopy)
tc$canopy <- round(tc$canopy,digits = 0)
tc$mean_canopy <- round(tc$mean_canopy,digits = 0)
#tc <- tc %>% group_by(ID_site) %>% arrange(ID) %>%
#  tidyr::fill(canopy, .direction="downup") %>%
#  tidyr::fill(mean_canopy, .direction="downup") %>%
#  ungroup()
#temporarliy change NA to zero to allow sum function to work
tc$canopy <-  tc$canopy %>% replace(is.na(.), 0)
tc$mean_canopy <-  tc$mean_canopy %>% replace(is.na(.), 0)
temp <- tc %>% mutate(canopy2 = canopy + mean_canopy) %>%
   ungroup()
temp <- temp %>% 
   dplyr::select(ID,ID_site,canopy2) %>% rename(canopy=canopy2)
#zeros back to NA
is.na(temp$canopy) <- !temp$canopy

FOR_canopy_means <- temp 


#without shrubs, mean canopy priority
# t1 <- t %>% group_by(ID_site.x) %>% arrange(ID) %>%
#   tidyr::fill(mean_canopy, .direction="downup") %>% 
#   ungroup()
# t1 <- t %>% rename(ID_site = ID_site.x) %>% 
#   mutate(canopy = coalesce(mean_canopy, tree_layer_cover)) 
# t1 <- t1 %>% group_by(ID_site) %>% arrange(ID) %>%
#   tidyr::fill(canopy, .direction="downup") %>% 
#   ungroup()

#without shrubs, tree_layer_cover priority
# t2 <- t %>% group_by(ID_site.x) %>% arrange(ID) %>%
#   tidyr::fill(tree_layer_cover, .direction="downup") %>% 
#   ungroup()
# t2<- t %>% rename(ID_site = ID_site.x) %>% 
#   mutate(canopy = coalesce(tree_layer_cover, mean_canopy)) 
# t2 <- t2 %>% group_by(ID_site) %>% arrange(ID) %>%
#   tidyr::fill(canopy, .direction="downup") %>% 
#   ungroup()
# 
# FOR_canopy_means <-t2 %>% group_by(ID) %>% summarise_at("canopy", mean)

saveRDS(FOR_canopy_means,"FOR_canopy_means.RDS")

#and IM data
tem <-
  IM_veg_cover %>% group_by(survey_year, country, ID, SPECIES) %>%
  summarise(mean_species = mean(VALUE)) %>% filter(survey_year >= 1994)
IM_veg_canopy <-
  tem %>% group_by(survey_year, country, ID) %>%
  summarise(sum_canopy = sum(mean_species))
IM_veg_canopy_means <-
  IM_veg_canopy %>% ungroup() %>% group_by(survey_year, country, ID) %>%
  summarise(mean_canopy = mean(sum_canopy)) %>%
  ungroup %>% dplyr::select(ID,mean_canopy) %>% rename(canopy = mean_canopy)

saveRDS(IM_veg_canopy_means,"IM_veg_canopy_means.RDS")

#combine
canopy_cover <-
  full_join(FOR_canopy_means, IM_veg_canopy_means) 
canopy_cover$canopy <- round(canopy_cover$canopy,digits = 0)

saveRDS(canopy_cover,"canopy_cover.RDS")

# VEG CWM Ellenbergs------------------------------------------------------------
# Non vascular plants site level-----------------------------------------------------
bryo_means <- veg_bryo %>% ungroup() %>% 
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover)) %>%
  spread(species_name, avg)

bryo_means_long <- veg_bryo %>%
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover))

#check for plots that have at least 2 species with data on Ellenberg values
Ell_check <- bryo_means_long %>% filter(species_name %in% ellenbergs$species_name) %>% 
  group_by(ID) %>% summarise(count = n_distinct(species_name)) %>%  
  filter(count < 2) 

remove_list <- Ell_check$ID

#take average cover per species per year/site combo (most plots have only one observation
#per year/site but some have multiple)
veg_bryo_e <- filter(veg_bryo, ID %!in% remove_list)
bryo_means_e <- veg_bryo_e %>% ungroup() %>% 
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover)) %>%
  spread(species_name, avg)

bryo_means_e <- as_tibble(bryo_means_e)
bryo_means_e[is.na(bryo_means_e)] <- 0 #replace NA with zero
bryo_means_e <-
  bryo_means_e[rowSums(bryo_means_e[, -c(1:3)]) != 0, ] #remove rows that sum to zero

# colnames(bryo_means[,404]) #last col is NA, remove
# bryo_means[, 404] <- NULL
bryo_means_matrix_e2 <- column_to_rownames(bryo_means, var = "ID") %>%
  select(-c(survey_year,country))
bryo_matrix_groups_e <- select(bryo_means_e, ID, survey_year, country)
bryo_means_matrix_e <- bryo_means_e[, -c(1:3)] #remove factor columns

#and without removing plots wth few ellenbergs but keep only those with at least 2 spp
veg_bryo2 <- filter(veg_bryo, richness >= 2)

bryo_means <- veg_bryo2 %>% ungroup() %>% 
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover))  %>% 
  spread(species_name, avg)

rich_mean_bryo <- veg_bryo %>%
  group_by(ID) %>%
  summarise(rich = mean(richness))
rich_mean_bryo$ID <- as.factor(rich_mean_bryo$ID)

bryo_means <- as_tibble(bryo_means)
bryo_means[is.na(bryo_means)] <- 0 #replace NA with zero
bryo_means <-
  bryo_means[rowSums(bryo_means[, -c(1:3)]) != 0, ] #remove rows that sum to zero

# colnames(bryo_means[,404]) #last col is NA, remove
# bryo_means[, 404] <- NULL
bryo_means_matrix2 <- column_to_rownames(bryo_means, var = "ID") %>%
  select(-c(survey_year,country))
bryo_matrix_groups <- select(bryo_means, ID, survey_year, country)
bryo_means_matrix <- bryo_means[, -c(1:3)] #remove factor columns

#make matrix of mean ellenberg values####
emptycols = c(which(colSums(bryo_means_matrix_e) == 0))
emptynames = names(emptycols)
nozeros = bryo_means_matrix_e %>% dplyr::select(-one_of(emptynames))
nozerosnames = colnames(nozeros)
# nozeros <- dplyr::select(nozeros, one_of(ellenbergs$species_name))
# nozeros <- filter(nozeros, rowSums(nozeros) > 0)
traits_nozeros <-
  ellenbergs[ellenbergs$species_name %in% nozerosnames,]
rownames(traits_nozeros) <- traits_nozeros$species_name

N = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "N",
  keyname = "species_name",
  method = 'mean'
)

T = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "T",
  keyname = "species_name",
  method = 'mean'
)

L = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "L",
  keyname = "species_name",
  method = 'mean'
)

F = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "F",
  keyname = "species_name",
  method = 'mean'
)

R = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "R",
  keyname = "species_name",
  method = 'mean'
)

temp <- cbind(N,T,L,F,R)
colnames(temp) <-  c("N","T","L","F","R")
temp <- as_tibble(temp)
temp <- dplyr::select(temp, L,F,R,N,T)
#CWM ellenbergs
CWM_bryo <- temp

#add back info columns
CWM_bryo$survey_year <-  as.numeric(as.character(bryo_means_e$survey_year))
CWM_bryo$country <- as.factor(bryo_means_e$country)
CWM_bryo$ID <- as.character(bryo_means_e$ID)

#some rows sum to zero, remove
CWM_bryo <- CWM_bryo %>% rownames_to_column()
zero_rows <- as.character(which(rowSums(CWM_bryo[, 2:6]) == 0))
CWM_bryo_nozero <- filter(CWM_bryo, rowname %!in% zero_rows)
CWM_bryo_nozero$rowname <- NULL
CWM_bryo <- CWM_bryo_nozero
CWM_bryo <-
  filter(CWM_bryo, survey_year >= "1994") #very few values from earlier, and those are low, potential skew

#Shannon diversity
library(vegan)

div_bryo_vec <- bryo_means_matrix %>% diversity(index = "shannon")

div_bryo_site <- cbind(div_bryo_vec,as.character(bryo_means$ID)) %>% as_tibble() %>% 
  dplyr::select(div_bryo_vec, V2) %>% rename(ID = V2)
div_bryo_site$div_bryo_vec <- as.numeric(div_bryo_site$div_bryo_vec)
div_bryo_site <- div_bryo_site %>% rename(div = div_bryo_vec)

#Simpson diversity
simp_div_bryo_site <- bryo_means_matrix %>% diversity(index = "simpson")
simp_div_bryo_site <- cbind(simp_div_bryo_site,as.character(bryo_means$ID)) %>% as_tibble() %>% rename(ID = V2) %>%
  rename( div_simp = simp_div_bryo_site)
simp_div_bryo_site$div_simp <- as.numeric(simp_div_bryo_site$div_simp)

div_bryo <- left_join(div_bryo_site, rich_mean_bryo, by = "ID") %>% 
  left_join(simp_div_bryo_site, by = "ID")
# Pielou's evenness (J):
div_bryo$J <- div_bryo$div/log(div_bryo$rich)

#add diversity to CWM ellenbergs
CWM_bryo <- full_join(CWM_bryo, div_bryo, by = "ID")
#CWM_bryo <- left_join(CWM_bryo, div_bryo, by = "ID") 

#colnames(CWM_bryo)[9] <- "div"
str(CWM_bryo)
CWM_bryo$ID <- as.factor(CWM_bryo$ID)
#CWM_bryo$div <- as.numeric(CWM_bryo$div)
str(CWM_bryo)

## Vascular plants ---------------------------------------------------------
#take average abundances per year/site combo
vasc_means <- veg_vasc %>%
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover)) %>%
  spread(species_name, avg)

rich_mean_vasc <- veg_vasc %>%
  group_by(ID) %>%
  summarise(rich = mean(richness))
rich_mean_vasc$ID <- as.factor(rich_mean_vasc$ID)

vasc_means_long <- veg_vasc %>%
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover))

vasc_means <- as_tibble(vasc_means)
vasc_means[is.na(vasc_means)] <- 0 #replace NA with zero
vasc_means <-
  vasc_means[rowSums(vasc_means[, -c(1:3)]) != 0, ] #remove rows that sum to zero

vasc_means_matrix <- vasc_means[, -c(1:3)] #remove factor columns

vasc_means_matrix[is.na(vasc_means_matrix)] <-
  0 #replace NA with zero
plot(colSums(vasc_means_matrix))
#plot(rowSums(vasc_means_matrix))
#a few outliers on colsums, investigate
dat <- vasc_means_matrix
dat[order(colSums(dat), decreasing = T)]
#it's a few very common species that are responsible for a large proportion of coverage
#eg vacciuniums, des flex...

#CWM:make matrix of mean ellenberg values per country/year
emptycols = c(which(colSums(vasc_means_matrix) == 0))
emptynames = names(emptycols)
nozeros = vasc_means_matrix %>% dplyr::select(-one_of(emptynames))
nozerosnames = colnames(nozeros)
traits_nozeros <-
  ellenbergs[ellenbergs$species_name %in% nozerosnames,]
rownames(traits_nozeros) <- traits_nozeros$species_name

y_vasc = isc(
  veg = as.matrix(vasc_means_matrix),
  trait.db = traits_nozeros,
  ivname = c("L", "F", "R", "N", "T"),
  keyname = "species_name",
  method = 'mean'
)

N = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "N",
  keyname = "species_name",
  method = 'mean'
)

T = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "T",
  keyname = "species_name",
  method = 'mean'
)

L = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "L",
  keyname = "species_name",
  method = 'mean'
)

F = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "F",
  keyname = "species_name",
  method = 'mean'
)

R = isc(
  veg = as.matrix(nozeros),
  trait.db = ellenbergs,
  ivname = "R",
  keyname = "species_name",
  method = 'mean'
)

temp <- cbind(N,T,L,F,R)
colnames(temp) <-  c("N","T","L","F","R")
temp <- as_tibble(temp)
temp <- dplyr::select(temp, L,F,R,N,T)
#CWM ellenbergs
CWM_vasc <- temp

CWM_vasc$survey_year <-
  as.numeric(as.character(vasc_means$survey_year))
CWM_vasc$country <- as.factor(vasc_means$country)
CWM_vasc$ID <- as.character(vasc_means$ID)

#some rows sum to zero, remove
CWM_vasc <- CWM_vasc %>% rownames_to_column()
zero_rows <- as.character(which(rowSums(CWM_vasc[, 2:6]) == 0))
CWM_vasc_nozero <- filter(CWM_vasc, rowname %!in% zero_rows)
CWM_vasc_nozero$rowname <- NULL
CWM_vasc <- CWM_vasc_nozero
CWM_vasc <-
  filter(CWM_vasc, survey_year >= "1994") #very few values from earlier, and those are low, potential skew

#CWM_vasc <- filter(CWM_vasc, survey_year >= "1995") #very few values from earlier, and those are low, potential skew

#Shannon diversity
div_vasc_site <- vasc_means_matrix %>% diversity(index = "shannon")
div_vasc_site <- cbind(div_vasc_site, vasc_means$ID) %>% as_tibble()
div_vasc_site <-rename(div_vasc_site,ID = V2)
div_vasc_site <-rename(div_vasc_site,div = div_vasc_site)
div_vasc_site$div <- as.numeric(div_vasc_site$div)
#add diversity as column
str(div_vasc_site)
CWM_vasc <- left_join(CWM_vasc, div_vasc_site) 
#colnames(CWM_vasc)[9] <- "div"
#str(CWM_vasc)
CWM_vasc$ID <- as.factor(CWM_vasc$ID)
CWM_vasc$div <- as.numeric(CWM_vasc$div)
str(CWM_vasc)


#Combine data for analysis####
#site level
VP <- CWM_vasc
VB <- CWM_bryo

VP2 <- left_join(CWM_vasc, rich_mean_vasc, by = "ID")
VB2 <- CWM_bryo

#merge datasets# add dep data
X <- full_join(VP2, VB2, by = "ID") %>%
  left_join (dep, by = "ID")

X$ID_site <- NULL

X <- mutate(X, ID_site = substr(X$ID, 6, 12))

X <- X %>%  left_join(location, by = "ID_site")
X <- dplyr::rename(X, country = country.x)
#x is vascular, y is bryophytes!
X <- as_tibble(X)

survey_year <- coalesce(X$survey_year.x, X$survey_year.y)
X$survey_year <- survey_year
X$survey_year.x <- NULL
X$survey_year.y <- NULL

#can't take mean of country, add country code number
X <- left_join(X, FOR_country, by = "country")
#X$country <- NULL

#add shading data
tem <-
  ungroup(canopy_cover) %>% dplyr::select(., ID, canopy)
X <- left_join(X, tem, by = "ID")
X <- rename(X, sum_canopy = canopy)
X <-
  X %>% mutate_if(is.character, as.factor) #change character variables to factors

#add temp/precip data
X2 <- left_join(X, meteo3, by = "ID")
X2$lon <- NULL
X2$lat <- NULL
#X2$mean_summer_temp <- X2$mean_summer_temp / 10
#X2$mean_winter_temp <- X2$mean_winter_temp / 10
#X2$mean_temp <- X2$mean_temp / 10
X <- X2
rm(X2)
temp2 <- X
temp2$ID_site.y <- NULL
temp2$country.y <- NULL
temp2 <- temp2 %>% rename(ID_site = ID_site.x)

#add conifer/decid group data from veg_tree_groups
veg_tree_groups <-
  full_join(FOR_veg_tree_groups, IM_veg_tree_groups)

temp2 <-
  left_join(temp2,
            dplyr::select(veg_tree_groups, species_name, grp_tree_species, ID),
            by = "ID")

#still have some missing data for grp_tree_species. Try imputing missing values from
#nearest neighbour (using lat long)####
#filter out the (very few) rows with no lat long data
temp2 <- temp2 %>% filter(!is.na(latitude))

data <-
  data.frame(
    latitude = temp2$latitude,
    longitude = temp2$longitude,
    grp_tree_species = temp2$grp_tree_species,
    ID = temp2$ID
  )
#data$grp_tree_species <-  as.character(data$grp_tree_species)
#Split the data frame up 
query = data[is.na(data$grp_tree_species), ]
tree_groups = data[!is.na(data$grp_tree_species), ]
library(FNN)
neighs = get.knnx(tree_groups[, c("latitude", "longitude")], query[, c("latitude", "longitude")], k =
                    1)
#Now insert the replacement classification directly into the data dataframe:

data[is.na(data$grp_tree_species), "grp_tree_species"] = tree_groups$grp_tree_species[neighs$nn.index]
# plot(data$longitude,
#      data$latitude,
#      col = data$grp_tree_species,
#      pch = 19)

extra_grp <- dplyr::select(data, grp_tree_species, ID)
temp2$grp_tree_species <- extra_grp$grp_tree_species
temp2 <- distinct(temp2)
temp2$grp_tree_species <- as.factor(temp2$grp_tree_species)
temp2$species_name <- as.factor(temp2$species_name)
temp2$country <- as.factor(temp2$country)
temp2$ID <- as.factor(temp2$ID)
temp2$ID_site <- as.factor(temp2$ID_site)

#add plain site ID
temp2$ID_siteonly <- temp2$ID_site
  temp2$ID_siteonly <- str_replace(temp2$ID_siteonly, "AT", "AT_") 
  temp2$ID_siteonly <- str_replace(temp2$ID_siteonly, "DE", "DE_") 
  temp2$ID_siteonly <- str_replace(temp2$ID_siteonly, "SE", "SE_") 

temp2$ID_site_only <- word(temp2$ID_siteonly, 2, sep = "_")
temp2$ID_siteonly <- temp2$ID_site_only
temp2$ID_site_only <- NULL
temp2$ID_siteonly <- as.factor(temp2$ID_siteonly)

#remove Canary islands - very little data, acts as outlier in many ways
temp2 <- filter(temp2, longitude>=-10)
#DONE!
  #We now have temp2, which is everything combined####

#Scale numeric variables####
non_nums <- select_if(temp2, negate(is.numeric))
nums <- as.data.frame(select_if(temp2, is.numeric) %>% scale())
temp2_scaled <- bind_cols(non_nums, nums)

#Functional diversity####
#By plot values to add alongside div.y etc in temp2####
bryophyte_func <-  readRDS("output_data/bryophyte_func.RDS")

#check for plots that have at least 2 species with data on Efunctional values
fun_check <- bryo_means_long %>% filter(species_name %in% bryophyte_func$species_name) %>% 
  group_by(ID) %>% summarise(count = n_distinct(species_name)) %>%  
  filter(count < 2) 

remove_list <- fun_check$ID

#take average cover per species per year/site combo (most plots have only one observation
#per year/site but some have multiple)
veg_bryo_f <- filter(veg_bryo, ID %!in% remove_list)
bryo_means_f <- veg_bryo_e %>% ungroup() %>% 
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover)) %>%
  spread(species_name, avg)

bryo_means_f <- as_tibble(bryo_means_f)
bryo_means_f[is.na(bryo_means_f)] <- 0 #replace NA with zero
bryo_means_f <-
  bryo_means_f[rowSums(bryo_means_f[, -c(1:3)]) != 0, ] #remove rows that sum to zero

#species data
# tem <- bryo_means
# tem <- filter(tem, ID %in% temp2$ID)#temp2?
# tem$survey_year <- as.numeric(as.character(tem$survey_year))
tem_matrix <- rename(bryo_means_f, rowname = ID) %>% column_to_rownames() %>%
  dplyr::select(.,-c(1:2))

x <- tem_matrix 
emptycols <-  c(which(colSums(x) == 0))
emptynames <- names(emptycols)
nozeros <- x %>% dplyr::select(-one_of(emptynames))
nozerosnames <-  colnames(nozeros)
bryo_func <- dplyr::select(bryophyte_func, -c(division, Family, L, T, K, F, R, N)) %>% 
  mutate_if(is.character, as.factor) %>% 
  rename(rowname = species_name) %>% arrange(rowname) %>% 
  column_to_rownames()
#bryo_func %>%  
traits_nozeros <- bryo_func[rownames(bryo_func) %in% nozerosnames, ]
traits_nozeros <- traits_nozeros[order(row.names(traits_nozeros)),]

#species with no available trait data
notrait <- setdiff(nozerosnames,rownames(traits_nozeros))
#dplyr::filter(bryo_means_long, species_name %in% notrait) %>% group_by(species_name) %>% 
#count() %>% arrange(desc(n)) %>% View()
#Plagiothecium curvifolium	
#remove species with no data from matrix
nozeros2 <- dplyr::select(nozeros, -one_of(notrait))
#now have some empty rows, remove
emptyrows <- c(which(rowSums(nozeros2) ==0))
emptynames <- names(emptyrows)
nozeros2 <- nozeros2[rowSums(nozeros2)!=0,]
#redo traits list to match nozeros2
nozerosnames <- colnames(nozeros2)
traits_nozeros <- bryophyte_func[bryophyte_func$species_name %in% nozerosnames, ]
traits_nozeros <- traits_nozeros %>% dplyr::select(-c(division, Family, L, T, K, F, R, N))
traits_nozeros <-  rename(traits_nozeros, rowname = species_name)
traits_nozeros <-  column_to_rownames(traits_nozeros) #traits_nozeros <- as.matrix(traits_nozeros)
traits_nozeros <- traits_nozeros %>% rownames_to_column() %>% as_tibble() %>% column_to_rownames()
traits_nozeros <- traits_nozeros %>%  rownames_to_column()  %>%
  mutate_if(is.character, as.factor) %>% 
  column_to_rownames()
traits_nozeros <- traits_nozeros[order(row.names(traits_nozeros)),]


traits_nozeros$life_strat2 <- NULL #mostly na's, not informative and creates problems with analysis
traits_nozeros$dist_gamet2 <- NULL
traits_nozeros <- dplyr::select(traits_nozeros,c(7:11))

#turn columns into factors unless binary which must remain numeric
col_names <- sapply(traits_nozeros, function(col) length(unique(col)) > 2)
traits_nozeros[ , col_names] <- lapply(traits_nozeros[ , col_names] , factor)

traits_nozeros$shoot_length <- as.numeric(traits_nozeros$shoot_length)

#some communities have <3 functionally singular species so FRic is problematic and gives warnings
#y <-  dbFD(traits_nozeros,nozeros2, corr = "lingoes", calc.FRic = FALSE)
#only morphological traits? as in not autoecological or sexual/veg regenerative traits...
traits_nozeros2 <- dplyr::select(traits_nozeros, growth_form, life_form, life_strat)
y2 <-  dbFD(traits_nozeros2,nozeros2, corr = "lingoes", calc.FRic = FALSE)

# FD_site <- as.data.frame(cbind(y$FEve, y$FDis, y$RaoQ)) 
# colnames(FD_site) <- c("FEve","FDis","RaoQ")
# FD_site <- rownames_to_column(FD_site) %>% rename(ID = rowname)

FD_site <- as.data.frame(cbind(y2$FEve, y2$FDis, y2$RaoQ)) 
colnames(FD_site) <- c("FEve","FDis","RaoQ")
FD_site <- rownames_to_column(FD_site) %>% rename(ID = rowname)

#try different distances: mahalanobis
#species by species trait distance
library(dummies)
dum_trait <- dummy.data.frame(traits_nozeros2)
x <- dum_trait
x.dist1 <- vegdist(x, method ="mahalanobis")
#x.dist <- FD::gowdis(dum_trait)
x.dist <- FD::mahaldis(as.matrix(dum_trait))

y2 <-  dbFD(x.dist,nozeros2, corr = "none", calc.FRic = FALSE)
FD_site2 <- as.data.frame(cbind(y2$FEve, y2$FDis, y2$RaoQ)) 
colnames(FD_site2) <- c("FEve_m","FDis_m","RaoQ_m")
FD_site2 <- rownames_to_column(FD_site2) %>% rename(ID = rowname)

FD_site <- left_join(FD_site, FD_site2, by ="ID")
rm(FD_site2)
saveRDS(FD_site, file = "FD_site.RDS")

#separate "environmental variables" dataframe for analyses####
env <- dplyr::select(temp2, country,ID_site, ID, n_nh4, n_no3,
       latitude, longitude, code_country, sum_canopy, mean_temp,
       mean_precip, survey_year, 
       grp_tree_species)

env %>% group_by(ID) %>% filter(n() > 1)#no duplictates

env_means <- env
env_means$plot.x <- NULL
env_means$ID_siteplot <- NULL
env_means$species_name <- NULL
env_means <- env_means %>% dplyr::group_by(ID) %>% summarise_if(is.numeric, mean) #put back country names and grp_tree?
env_means <- env_means[!duplicated(env_means), ]
env_means_scaled <- scale(dplyr::select(env_means, -c("ID")))  %>%
  as.data.frame() %>% mutate(ID = env_means$ID) 
env_means_scaled$n <- NULL
env_means_scaled$p <- NULL
env_means_scaled$N.P <- NULL

bryo_env_means <- env_means %>% filter(ID %in% bryo_means$ID)

#temp2 <- rename(temp2, div.y = div)
#Generate data3, per site values####
data3 <- temp2 %>% ungroup() %>% group_by(ID) %>%
  dplyr::select(
    L.y,
    R.y,
    N.y,
    F.y,
    T.y,
    L.x,
    rich.y,
    div.y,
    div_simp,
    J,
    n_nh4,
    n_no3,
    latitude,
    longitude,
    survey_year,
    mean_precip, 
    mean_temp,
    sum_canopy,
    ID_site,
    ID_siteonly,
    ID,
    country
  ) %>% 
  left_join(., FD_site, by = "ID") %>% 
  filter(., latitude >=45.0) %>% filter (longitude > -2.5) %>%
  filter(longitude < 30) %>% 
  #drop_na(N.y) %>%
  drop_na(n_nh4)
#data3 <- left_join(data3, canopy_cover)
data3 <- ungroup(data3)
#data3 <- rename(data3, sum_canopy = canopy)
#data3 <- left_join(data3, p_mean)
data3$ID <- as.factor(data3$ID)
data3 <- dplyr::select(data3, country, ID_site, ID_siteonly, ID, everything())

#data3 <- filter(data3, L.y > 1)
#data3 <- filter(data3, N.y > 1)
#data3 <- filter(data3, L.x > 1)
#data3 <- filter(data3, RaoQ > 0)
data3 <- filter(data3, div.y > 0.01)

#add group of dominant tree species, 1 = broadleaf 2 = conifer
Mode1 <- function(codes){
  which.max(tabulate(codes))
}
t <- dplyr::select(temp2, ID, grp_tree_species) %>% group_by(ID) %>% summarise(tree = Mode1(grp_tree_species))
t <- t %>% filter(!duplicated(.))

data3 <- left_join(data3, t, by = "ID")
data3$ID <- as_factor(data3$ID)
data3$tree <- as_factor(data3$tree)
data3$div.y <- as.numeric(as.character(data3$div.y))

#Add age etc data from site info SI file for Forests####
si_stand_decription <- read_delim("raw_data/Veg_Forests/si_stand_decription.csv", 
                                  ";", escape_double = FALSE, col_types = cols(change_date = col_skip(), 
                                                                               code_line = col_skip(), code_nfi_status = col_skip(), 
                                                                               code_plot_status = col_skip(), line_nr = col_skip(), 
                                                                               other_obs = col_skip(), partner_code = col_skip()), 
                                  trim_ws = TRUE)
si_stand_decription <-
  transform(
    si_stand_decription,
    ID = paste0(
      si_stand_decription$survey_year,
      sep = "_",
      si_stand_decription$code_country,
      sep = "_",
      si_stand_decription$code_plot
    )
  )

si_stand_decription <-
  transform(
    si_stand_decription,
    ID_site = paste0(
      si_stand_decription$code_country,
      sep = "_",
      si_stand_decription$code_plot
    )
  )
# library(DataExplorer)
# create_report(select(si_stand_decription, -c(ID))) 
# lots of missing data for site info
site_info <- dplyr::select(si_stand_decription, survey_year, ID, ID_site,
                           code_mean_age, code_forest_type, canopy_closure)
site_info <- rename(site_info, mean_age=code_mean_age)
#site_info$survey_year <- as.Date(site_info$survey_year)

#mixed age stands, code 8, how to deal with?
#%>% mutate(mean_age = replace(mean_age, which(mean_age == 8), NA)) 
tem <- site_info  %>%
  mutate(mean_age = replace(mean_age, which(survey_year == 2011 & mean_age == 1), NA)) %>% 
  group_by(variable=factor(ID_site)) %>% 
  filter(mean_age != 0) %>% 
  # Keep the first row of each variable, after sorting by Date
  # This gives us the first non-zero row
  arrange(survey_year) %>% 
  slice(1) %>% 
  ungroup() %>% 
  complete(ID_site)

s <- c("1994-09-08")
e <- c("2018-09-08")
tem$start <- s
tem$end <- e
tem$date1 <- as.Date(as.character(tem$survey_year), format = "%Y")
tem$start <- as.Date(tem$start)
tem$end <- as.Date(tem$end)
#tem <- tem %>% mutate(mean_age = replace(mean_age, which(mean_age == 8), NA)) #8 is mixed stands

  # sequence of monthly dates for each corresponding start, end elements
tem2 <-  tem %>%  mutate(ID_site, date1 = map2(start, end, seq, by = "1 year")) %>%
  # unnest the list column
  unnest %>% 
  # remove any duplicate rows
  distinct() %>% dplyr::select(.,-c(start, end))

tem2$survey_year2=(substr(tem2$date1,0,4))
tem2$survey_year2 <- as.numeric(tem2$survey_year2)
#construct age classes related to the earliest class we have in the data
tem2$dif <- tem2$survey_year2 -tem2$survey_year
#assume midpoint of age class range, then 11 years up will take it up a class, 11 years down, down.
#or conservative assumption and require 20 years diff? histogram looks much more like site_info with latter.

#8 is not a real age class so cannot be adjusted in the same way as 1-7
temna <- filter(tem2, mean_age == 8)
tem2 <- filter(tem2, mean_age != 8)
tem2 <- tem2 %>% mutate(
  mean_age2 = case_when(dif > 20  ~ mean_age+1,
          dif < -20 ~ mean_age-1))
tem2$mean_age3 <- coalesce(tem2$mean_age2, tem2$mean_age)
#mean age cannot be less than 1 or greater than 7 (class 8 is mixed age stands, 
#genuine 8's were removed earlier)
tem2$mean_age3[tem2$mean_age3 > 7] <- 7
tem2$mean_age3[tem2$mean_age3 < 1] <- 1
tem2$mean_age <- NULL
tem2 <- rename(tem2, mean_age=mean_age3)
#add back genuine 8's ? 
tem2 <- full_join(tem2, temna)
tem2 <- dplyr::select(tem2,c(ID_site, code_forest_type, canopy_closure,survey_year2, mean_age)) %>% 
  rename(survey_year=survey_year2) 
tem2$mean_age <- as.character(tem2$mean_age)
#convert to factor using mid-point of each range
tem2$mean_age <- recode(tem2$mean_age,"1"="10", "2"="30","3"="50","4"="70","5"="90","6"="110",
      "7"="130","8" = "mixed")

#regenerate ID (needed)
tem2 <- transform(
  tem2,
  ID = paste0(
    tem2$survey_year,
    sep = "_",
    tem2$ID_site
  ))

#now with age class
data4 <- left_join(data3, dplyr::select(tem2, mean_age, ID), by= "ID")
#manually add some missing ICP IM ages
data4 <- data4 %>% mutate(mean_age=replace(mean_age, ID_site=="SE14", "130")) 
data4 <- data4 %>% mutate(mean_age=replace(mean_age, ID_site=="SE15", "130")) 
data4 <- data4 %>% mutate(mean_age=replace(mean_age, ID_site=="SE16", "130")) 
data4 <- data4 %>% mutate(mean_age=replace(mean_age, ID_site=="SE04", "130")) 
data4 <- data4 %>% mutate(mean_age=replace(mean_age, ID_site=="AT01", "130")) 
data4 <- data4 %>% mutate(mean_age=replace(mean_age, ID_site=="DE01", "130")) 
#make factor?
data4$mean_age <- as.factor(data4$mean_age) 
#data4$mean_age <- ordered(data4$mean_age, levels = c("10","30","50","70","90","110","130"))

#check missing canopy
missing_can <- filter(data4, is.na(sum_canopy))

# data5bak <- data5
data5 <- data4
q <- quantile(data5$n_nh4, .99, na.rm = TRUE)
data5 <- data5 %>% filter(n_nh4 < q)

q <- quantile(data5$n_no3, .99, na.rm = TRUE)
data5 <- data5 %>% filter(n_no3 < q)

#keep only columns we need
#datc <- dplyr::select(data5,country,ID,ID_site,ID_siteonly,N.y,div.y,rich.y,RaoQ,RaoQ_m,n_nh4,n_no3,
#               sum_canopy,div_simp,J,L.x,mean_temp,mean_precip,latitude,longitude,survey_year,tree,mean_age)

#unique species by country- find and remove countries with poor coverage of bryophytes
# Dsum<-veg_bryo %>% group_by(country) %>% summarise(unique_spp = n_distinct(species_name)) %>% 
#   arrange(unique_spp)%>% mutate(row = row_number())
# ggplot(Dsum,aes(x=reorder(country,row),y=unique_spp))+geom_point()+geom_hline(yintercept=20)+xlab("country")+
#   theme(axis.text.x=element_text(angle=90, hjust=1))
# ggplot(Dsum,aes(x=reorder(country,row),y=log(unique_spp)))+geom_point()+geom_hline(yintercept=log(20))+xlab("country")+
#   theme(axis.text.x=element_text(angle=90, hjust=1))

dat.all.countries <- data5
#dat.all.countries$sum_canopy[dat.all.countries$sum_canopy >100] <- 100#
dat.all.countries <- dat.all.countries[!is.na(dat.all.countries$n_nh4),]

#ghost levels in factors to remove (e.g countries not in dataset)
dat.all.countries$country <-  dat.all.countries$country %>% as.character() %>% as_factor()
dat.all.countries$ID_site <-  dat.all.countries$ID_site %>% as.character() %>% as_factor()
dat.all.countries$ID_siteonly <-  dat.all.countries$ID_siteonly %>% as.character() %>% as_factor()
dat.all.countries$ID <-  dat.all.countries$ID %>% as.character() %>% as_factor()

str(dat.all.countries)
library(VIM)
matrixplot(dat.all.countries, interactive = F, sortby = "sum_canopy") #VIM

# q <- quantile(dep_test$n_nh4_s, .99, na.rm = TRUE)
# data5 <- dep_test %>% filter(n_nh4_s < q)
# 
# q <- quantile(dep_test$n_no3_s, .99, na.rm = TRUE)
# data5 <- dep_test %>% filter(n_no3_s < q)
# data5bak <- dep_test

#Export data####
#write_csv2(dat.all.countries, file="dat_all_countries2.csv")
# write_csv2(dat.select.countries, file="dat_select_countries.csv")
# saveRDS(temp2, file = "temp2.RDS")
# saveRDS(data3, file = "data3.RDS")
# saveRDS(data4, file = "data4.RDS")
# saveRDS(data5, file = "data5.RDS")
#saveRDS(dat.select.countries, file = "dat.select.countries.RDS")
saveRDS(dat.all.countries, file = "dat.all.countries.RDS")
#saveRDS(dat.all.countries, file = "dat.all.countries.noex.RDS")
# saveRDS(veg, file = "veg.RDS")
# saveRDS(veg_bryo, file = "veg_bryo.RDS")
# saveRDS(env_means, file = "env_means.RDS")
saveRDS(bryo_means, file = "bryo_means.RDS")
#saveRDS(bryo_env_means_scaled, file = "bryo_env_means_scaled.RDS")

# saveRDS(vasc_plots_wide, file = "vasc_plots_wide.RDS")
# saveRDS(vasc_plots_matrix, file = "vasc_plots_matrix.RDS")
# saveRDS(bryo_plots_wide, file = "bryo_plots_wide.RDS")
# saveRDS(bryo_plots_matrix, file = "bryo_plots_matrix.RDS")
#saveRDS(meteo2, file = "meteo2.RDS")

dat <- dat.all.countries#
dat <- filter(dat, mean_age != "mixed") %>%
  mutate(mean_age = as.numeric(as.character(mean_age)))#
dat <- dplyr::select(dat, -c(ID_siteonly,country, rich.y, J, L.x, div.y, FEve,
                             FDis, FEve_m, FDis_m, RaoQ_m, L.y, T.y, R.y, F.y))
dropdat <- filter(dat, !is.na(sum_canopy))
dat <- dropdat#
#few plots have rao of zero- if all species have the same functional values they occupy the
#same space in trait space so have the same dispersion. But we log transform so need to 
#add tiny value to allow this
#dat$RaoQ <- ifelse(dat$RaoQ < 0.000001, 0.000001, dat$RaoQ)
#dat$RaoQ_m <- ifelse(dat$RaoQ_m < 0.000001, 0.000001, dat$RaoQ_m)
#dat$Rao_logit <- log(dat$RaoQ/(1-dat$RaoQ))
#dat$div_simp_logit <- log(dat$div_simp/(1-dat$div_simp))

#change names to make plots nicer
dat_names <- dat %>% rename(NH4=n_nh4, NO3=n_no3, canopy=sum_canopy,
                            precipitation = mean_precip,
                            temperature = mean_temp, age = mean_age,
                            diversity = div_simp, Ellenberg_N = N.y )

#remove outlier 4_1406?
dat<-dat[dat$ID_site!="4_1404",] #
dat_names<-dat_names[dat_names$ID_site!="4_1404",] 

saveRDS(dat, file = "dat.RDS")
saveRDS(dat_names,file = "dat_names.RDS")


#bryo_means
# bryo_means <- veg_bryo %>% ungroup() %>% 
#   group_by(survey_year, country, ID, species_name) %>%
#   summarise(avg = mean(cover)) %>%
#   spread(species_name, avg)

saveRDS(bryo_means, file = "bryo_means.RDS")

# bryo_env_means_scaled$n <- NULL
# bryo_env_means_scaled$p <- NULL
# bryo_env_means_scaled$N.P <- NULL

#dropdat <- filter(dat, !is.na(sum_canopy))
#datcan <- dat
#datcan$sum_canopy <-ifelse(datcan$sum_canopy > 100, 100, datcan$sum_canopy)
#dat <- drop_na(dat)
#dat <- dat.select.countries

#Outlier plots
# create the subset
#Plot.1404 <- subset(dat, ID_site =="4_1404")
#outliers <- filter(dat, n_nh4 >25 & n_no3 >20)
#outlier <- subset(dat, ID %in% outliers$ID)
# plot the data
#ggplot(dat,aes(x=n_nh4,y=n_no3))+geom_point()+geom_rug() +  
#  geom_point(data=Plot.1404, colour="red") +  
#  geom_text(data=Plot.1404, label="Plot.1404", vjust=1) #+

#Plot.1404 <- subset(dat.all.countries, ID_site =="4_1404")
#outliers <- filter(dat.all.countries, n_nh4 >25 & n_no3 >20)
#outlier <- subset(dat.all.countries, ID %in% outliers$ID)
# plot the data
#ggplot(dat.all.countries,aes(x=n_nh4,y=n_no3))+geom_point()+geom_rug() +  
#  geom_point(data=Plot.1404, colour="red") +  
#  geom_text(data=Plot.1404, label="Plot.1404", vjust=1) #+
#geom_point(data=outlier, colour="yellow") +  
#  geom_text(data=outlier, label="outs", vjust=1) 

#canopy change plot
#ggplot(dat, aes(x=survey_year, y= sum_canopy))+
#  geom_point()+
#  geom_smooth(method = "loess", se =FALSE) + facet_wrap( ~ ID_site)

#ggplot(dat.all.countries.old, aes(x=survey_year, y= N.y))+
#  geom_point()+
#  geom_smooth(method = "loess", se =FALSE) + facet_wrap( ~ ID_site)



#STOP####