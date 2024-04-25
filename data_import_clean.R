library(tidyverse)
library(lubridate)
library(vegdata)
#library(FD)
#library(labdsv)

`%!in%` = Negate(`%in%`)

# Import location data ####
#Import IM locations (in decimal degrees)
IM_locations <- read_csv("raw_data/From_paper2/Veg_IM/IM_locations.csv")
IM_locations <-
  dplyr::select(IM_locations, Code, Latitude, Longitude) %>%
  rename(ID_site = Code,
         latitude = Latitude,
         longitude = Longitude)

#Import VG FD (Intensive plots and Forest damage) #1990 - 2018 ####
FD_VG <- read_csv("raw_data/ICPIM_FD_VG.csv", col_types = cols(Institute = col_skip()))
VG_orig <- filter(FD_VG, Subprog == "VG")

VG <- read_csv("raw_data/VG_2020.csv", col_types = cols(Institute = col_skip()))
VG <- VG %>% rename("Value" = "SampleValValue")
#break down to FD and VG
#FD <- filter(FD_VG, Subprog == "FD")

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
      "CH" = "Switzerland",
      "CZ" = "Czechia",
      "DE" = "Germany",
      "EE" = "Estonia",
      "ES" = "Spain",
      "FI" = "Finland",
      "IT" = "Italy",
      "LT" = "Lithuania",
      "LV" = "Latvia",
      "SE" = "Sweden",
      "IE" = "Ireland",
      "NO" = "Norway",
      "PL" = "Poland",
      "RU" = "Russia"
      
    )
  )
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

###Make exclude list e.g sites with only one year of data
exclude <- c("AT01", "CH02", "CZ02", "IT01","IT02","IT03","IT04","IT05","IT06",
             "IT07","IT08","IT09","IT10","IT11","IT12", "LT03", "ES02", "DE02", "RU04", "RU13")
VG <- VG %>% filter(!ID_site %in% exclude)

#there is a species named "NULL", remove
VG <- VG %>% filter(!Medium == "NULL")

#filter to keep only field and bottom layers
VG_ground <- VG %>% filter(VG$Parameter %in% c("COVE_F", "COVE_B")) %>%
  rename(Param = Parameter)
#VG_ground <- select(VG_ground, -c("AreaName", "StationName", "Descrption", "country", "SpeciesName", "Class"))
VG_ground <- select(VG_ground, -c("AreaName", "StationName", "country", "SpeciesName"))
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

#Make data wide format
VG_g_wide <- VG_ground %>% group_by(ID,ID_fine,ID_fine2) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Cover)
VG_g_wide <- VG_g_wide %>% replace(is.na(.), 0) %>% ungroup()

VG_g_wide$ID_plot <- as.factor(VG_g_wide$ID_plot)
VG_g_wide$ID_subplot <- as.factor(VG_g_wide$ID_subplot)

#Also group by intensive plot and use means rather than individual small squares####
VG_ground_pl <- VG_ground %>% group_by(survey_year,ID,ID_site,ID_plot,ID_fine, Species) %>% 
  summarise(Cover = mean(Cover)) %>% ungroup()

#Make wide format
VG_g_pl_wide <- VG_ground_pl %>% group_by(ID,ID_fine) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Cover)
VG_g_pl_wide <- VG_g_pl_wide %>% replace(is.na(.), 0) %>% ungroup()

#find sites that have no data for subplots (are they reporting overall mean??)
few_g <- VG_g_wide %>% group_by(ID_fine) %>% summarise(distinctID2 = n_distinct(ID_fine2)) %>% filter(distinctID2 < 2)

#make factors
VG_g_wide$ID <- as_factor(VG_g_wide$ID)
VG_g_wide$ID_site <- as_factor(VG_g_wide$ID_site)
VG_g_wide$ID_fine <- as_factor(VG_g_wide$ID_fine)
VG_g_wide$ID_fine2 <- as_factor(VG_g_wide$ID_fine2)

#later functions complain about zero sum rows, check and remove
VG_g_wide <- VG_g_wide[rowSums(VG_g_wide[,-c(1:7)]) != 0, ]
VG_g_wide$StationCode <- NULL
VG_g_wide$SpatialPool <- NULL
VG_g_wide$Param <- NULL
VG_gm <- column_to_rownames(VG_g_wide, var = "ID_fine2")
VG_gm <- select(VG_gm, -c(1:6))

#And the same for the per intensive plot version
#make factor
VG_g_pl_wide$ID <- as_factor(VG_g_pl_wide$ID)
VG_g_pl_wide$ID_site <- as_factor(VG_g_pl_wide$ID_site)
VG_g_pl_wide$ID_fine <- as_factor(VG_g_pl_wide$ID_fine)
VG_g_pl_wide$ID_plot <- as_factor(VG_g_pl_wide$ID_plot)

#later functions complain about zero sum rows, check and remove
VG_g_pl_wide <- VG_g_pl_wide[rowSums(VG_g_pl_wide[,-c(1:5)]) != 0, ]
VG_gm2 <- dplyr::select(VG_g_pl_wide, -c(1:4))
VG_gm2 <- column_to_rownames(VG_gm2, var = "ID_fine")

#rename for use in ordinations
#Plot level
ordispe_pl <- VG_gm2
ordispe_pl_extras <- VG_g_pl_wide
#Subplot level
ordispe_sub <- VG_gm
ordi_extras <-  select(ungroup(VG_g_wide), ID_fine2, ID_plot, ID_site, survey_year,ID, ID_fine) %>% column_to_rownames(var="ID_fine2")


# STOP HERE for exploratory veg data only #
# Continue for environmental variables as well #

#Environmental Variables#####
ICPIM_AC_AM_TF <- read_csv("raw_data/ICPIM_AC_AM_TF.csv", 
                           col_types = cols(Day = col_skip(), DataQualityFlag = col_skip(), 
                                            StatusFlag = col_skip(), NeedleAgeCode = col_skip(), 
                                            ParameterInfo = col_skip(), PretreInfo = col_skip(), 
                                            DeterInfo = col_skip()))
library(lubridate)
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
IM_TF<- IM_TF %>% group_by(ID,country,ID_site,survey_year, YearMonth,StationCode, ParameterCode) %>%
  summarise_at(vars(Value), mean, na.rm = TRUE) 

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
  summarise_at(vars(NH4M, NO3M, SO4SM, PREC), sum, na.rm = TRUE) 

IM_dep_means <- IM_dep_means %>% group_by(ID,ID_site,survey_year) %>%
  summarise_at(vars(NH4M, NO3M, SO4SM, PREC), mean, na.rm = TRUE) 

IM_dep2 <- drop_na(IM_dep_means) 
IM_dep2$survey_year <- NULL

#Temp
IM_AM <- select(IM_AM,ID, ID_site,survey_year, country, StationCode,
                YearMonth, ParameterCode, Value)
#take monthly mean for variables with more than one measurement in a month recorded
IM_AM <- IM_AM %>% group_by(ID,country,ID_site,survey_year, YearMonth,StationCode,
                            ParameterCode) %>%
  summarise_at(vars(Value), mean, na.rm = TRUE) 
IM_AM <- pivot_wider(IM_AM, names_from = ParameterCode, values_from = Value)
#IM_AM is monthly if you need eg summer vs winter avergae temps
#make annual mean
IM_meteo <- IM_AM %>% group_by(ID,ID_site,survey_year) %>%
  summarise_at(vars(TEMP), mean, na.rm = TRUE) %>% ungroup()

##Temp for missing data from copernicus database
# library(ncdf4)
# mycdf <- nc_open(file.choose(), verbose = TRUE, write = FALSE)
# timedata <- ncvar_get(mycdf,'time')

#library(RNCEP)
##Define limits for latitude and longitude
# min_lat <- min(all_data$latitude, na.rm = TRUE)
# max_lat <- max(all_data$latitude, na.rm = TRUE)
# 
# min_lon <- min(all_data$longitude, na.rm = TRUE)
# max_lon <- max(all_data$longitude, na.rm = TRUE)

# define arguments for latitude and longitude
# lat_range <- c(min_lat, max_lat)
# lon_range <-c(min_lon, max_lon)
# 
# # get monthly air temperature between 2001 and 2018
# weather <- NCEP.gather(variable = "air.sig995", level = "surface", months.minmax = c(1,12),
#                        years.minmax = c(1998,2018), lat.southnorth =lat_range,
#                        lon.westeast = lon_range)
# #annual means
# weather2 <- NCEP.aggregate(weather, YEARS = TRUE, MONTHS = FALSE,
#                            DAYS = FALSE, HOURS = FALSE, fxn = 'mean')
# t98 <- as.data.frame(weather[,,1])
# #etc
# 
# tempa <- select(all_data, ID_fine, latitude, longitude, survey_year) %>% distinct()
# #round to nearest 2.5
# mround <- function(x,base){ 
#   base*round(x/base) 
# } 
# 
# tempa$latituder <- mround(tempa$latitude, 2.5)
# tempa$longituder <- mround(tempa$longitude, 2.5)
# tempa$TEMP2 <- as.double(c(1:85))
# tempa[85,7] <- 267.5
# bak <- tempa
# #in kelvin, convert to centigrade
# tempa$TEMP_C <- tempa$TEMP2
# tempa$TEMP_C <- tempa$TEMP_C - 273.15


#interpolate

# lat <- all_data$latitude
# long <- all_data$longitude
# year <- all_data$survey_year
# ID_fine <- as.character(all_data$ID_fine)
# 
# temps <- cbind(ID_fine, lat,long, year) %>% as.data.frame()
# temps$date <- ymd_hms(paste0(temps$year, '-01-01 00:00:01'))
# lat <- as.numeric(temps$lat)
# lon <- as.numeric(temps$long)
# dt <- as.character(temps$date)
#or apply interpolate function
# interp <- NCEP.interp('air.sig995','surface', lat, lon, dt, reanalysis2 = FALSE,
#             interpolate.space = TRUE, interpolate.time = TRUE,
#             keep.unpacking.info = FALSE, return.units = TRUE,
#             status.bar=TRUE) 
#in kelvin, convert to centigrade
# interp_c <- interp - 273.15
# temps$temp2 <- interp_c
# temps2 <- distinct(temps)

# extract longitude & latitude based on created weather dataset
# lat <- dimnames(weather)[[1]] # in increments of 2.5
# lon <- dimnames(weather)[[2]] # in increments of 2.5

#ENV combo####
IM_env <- left_join(IM_dep2, select(IM_meteo, ID, TEMP), by = "ID") %>% ungroup() %>% 
  left_join(., IM_locations) %>% 
  select(., -c(ID_site))


#create data - species and environemntal matrixes####
#There is species data back to 90 but only env data from 98 on...
test <- left_join(VG_g_wide, IM_env, by = "ID") %>% 
  #filter(NH4M > 0) %>% #POOR MATCH, MANY HAVE VG BUT NO DEP DATA??
  ungroup()

#test <- filter(test, ID_plot %!in% c("CH02_0001"))
#drop unused facotr levels
test$ID <- as.factor(test$ID)
test$ID_fine2 <- droplevels(test$ID_fine2)
test$ID_fine <- droplevels(test$ID_fine)
test$ID_site <- droplevels(test$ID_site)
test$ID_plot <- droplevels(test$ID_plot)

#date range for plots
last <- test %>% group_by(ID_plot) %>% summarise(last=max(survey_year))
first <- test %>% group_by(ID_plot) %>% summarise(first=min(survey_year))
range <- left_join(first, last)
range <- range %>% mutate(diff= last-first)
View(range)

#make matrix for ordinations
ordi_extras <-  select(ungroup(test), ID_fine2, ID_plot, ID_site, NH4M, NO3M, SO4SM, TEMP, latitude, longitude, 
                   PREC, survey_year) %>% column_to_rownames(var="ID_fine2")

ordi_extras.df <-  select(ungroup(test), ID_fine2, ID_plot, ID_site, NH4M, NO3M, SO4SM, TEMP, latitude, longitude, 
                       PREC, survey_year)
extras <- test %>% 
  select(ID, ID_site,ID_fine, ID_fine2, ID_plot, ID_subplot,survey_year)

ordispe.df <- select(test, -c(survey_year,ID,ID_site,ID_plot,ID_subplot,ID_fine, 
                              NH4M,NO3M,SO4SM,PREC,TEMP,latitude,longitude)) %>%
  column_to_rownames(var="ID_fine2")
                     
#ordispe <- select(test, -c(2:7, 331:337)) %>% column_to_rownames(var="ID_fine2") 
ordispe.df <- ordispe.df[,colSums(ordispe.df !=0)>0] #drop zero sum columns (spp not present anywhere)

# generate species sums these will be NA is any are missing
csum <- colSums(ordispe.df)
# check if any are missing
any(is.na(csum))
# yes, some missing, so which ones?
which(is.na(csum))

ordienv <-  select(test, ID_fine2, NH4M, NO3M, SO4SM, latitude, longitude, 
                   PREC, survey_year) %>% column_to_rownames(var="ID_fine2")
ordispe <- as.matrix(ordispe.df)
ordienv <- as.matrix(ordienv)

#Scale numeric variables####
ordienv_scaled <- ordienv %>% scale()

