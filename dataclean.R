`%!in%` = Negate(`%in%`)

library(tidyverse)
library(lubridate)
library(vegdata)
#library(runjags)
#library(mcmcplots)
library(FD)
library(labdsv)
library(vegan)
library(ggfortify)
library(dave)
library(TR8)
library(magrittr)
#library(car)#vif
#install.packages("DataExplorer")
#library(DataExplorer)
#library(mice)#dealing with missing data

##RUN FUNCTIONS FILE FIRST!####

# FC - Foliar chemistry data ----------------------------------------------

# Import/clean IM FC ------------------------------------------------------

ICPIM_FC_data <- read.csv("raw_data/IM_FC_data.csv")

#make new dataframe with only relevant columns
IM_FC <-
  dplyr::select(
    ICPIM_FC_data,
    AreaCode,
    AreaName,
    Medium,
    YearMonth,
    ParameterCode,
    Value,
    Unit,
    NeedleAgeCode,
    Sample_ID
  )
#ParameterCode has whitespace after entries, strip out
IM_FC$ParameterCode <-  trimws(IM_FC$ParameterCode, "right")
#Don't need information on month just year? Add column for year only and keep both just in case...
IM_FC$Year <- as.integer(substr(IM_FC$YearMonth, 0, 4))

#once subset into species and nutrient types there are unequal p and n observations
#i.e 7 more p: can remove "extras" to allow more analyses
n <- filter(IM_FC, ParameterCode == "NTOT")
p <- filter(IM_FC, ParameterCode == "PTOT")
cutlist <- setdiff(p$Sample_ID, n$Sample_ID)
p <- filter(p, Sample_ID %!in% cutlist)
#put back together
IM_FC_npequal <- full_join(p, n)

#make data wide and create NP ration column
n.p.ratio <-
  spread(IM_FC_npequal, ParameterCode, Value) %>% mutate(N.P = NTOT / PTOT)

#set main dataframe to be the new version with ratio column
IM_FC <- n.p.ratio

#harmonise labelling with Forests data.
IM_FC <-
  IM_FC %>% rename(
    p = PTOT,
    n = NTOT,
    code_leaves_type = NeedleAgeCode,
    survey_year = Year,
    species_name = Medium
  )

#remove factor
IM_FC$species_name <- as.character(IM_FC$species_name)
IM_FC$code_leaves_type <- as.character(IM_FC$code_leaves_type)

IM_FC$species_name <-
  IM_FC$species_name %>% str_replace_all(
    c(
      "PICE ABI" = "Picea abies",
      "FAGU SYL" = "Fagus sylvatica",
      "PINU SYL" = "Pinus sylvestris",
      "QUER ILE" = "Quercus ilex",
      "BETU PUB" = "Betula pubescens",
      "QUER PET" = "Quercus petraea",
      "QUER CER" = "Quercus cerris"
    )
  )

decid_im <-
  c(
    "Fagus sylvatica" ,
    "Quercus ilex" ,
    "Betula pubescens",
    "Quercus petraea" ,
    "Quercus cerris"
  )
IM_FC <-
  mutate(IM_FC,
         grp_tree_species = ifelse(species_name %in% decid_im, "broadleaves", "conifers")) #add grouping column

#add country
IM_FC <-
  mutate(IM_FC, country = str_sub(IM_FC$AreaCode, 1, 2)) #keep first two letters of country code and use to match
IM_FC$country <-
  IM_FC$country %>% str_replace_all(
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

#NULL for IM means not a needle i.e. broadleaves. Forests data records age regardless, so to combine
#we count NULL as same year (which they are by definition in deciduous trees)
IM_FC$code_leaves_type <- IM_FC$code_leaves_type %>%
  str_replace_all(c(
    "C" = "0",
    "P" = "1",
    "E" = "2",
    "NULL" = "0"
  )) %>%
  as.integer(.)

#make line ID to match with deposition data
IM_FC <- transform(IM_FC,ID = paste0(IM_FC$survey_year, sep = "_",  IM_FC$AreaCode))


# Combine IM and FOR FC ---------------------------------------------------

FC <- full_join(IM_FC, FOR_FC)
#FC <- FOR_FC
#plot(jitter(FC$survey_year,2), jitter(FC$N.P, 2), cex = 0.1) #quick graphical check
#there is a small cluster of minus values for n and p,from germany (plot 308, 611, 613) must be in error, remove them
FC <- FC[FC$p >= 0,]
FC <- FC[FC$n >= 0,]
FC <- as_tibble(FC)

FC <-
  FC %>% filter(rowSums(is.na(.)) != ncol(.)) #remove all all NA rows

#sort out dates
FC$survey_date <-
  as.Date(paste0("0101", FC$survey_year), format = "%d%m%Y")

FC$survey_year <-  year(FC$survey_date)#lubridate
#FC$survey_year <- as.factor(FC$survey_year)#messes up some ggplots when factor
fc.factor.year <- as.factor(FC$survey_year)

#some samples (up to 2001) are conifers but have NA for code year. Assume current year?
#some broadleafs are recorded NA, some 0. Change all <1 to 0?
filter(FC, grp_tree_species == "conifers" & is.na(code_leaves_type))
#replace na with 0
FC <-
  FC %>% mutate(code_leaves_type = replace(code_leaves_type, which(is.na(code_leaves_type) &
                                                                     survey_year > 0) , 0))

#& grp_tree_species == "conifers") 

#make some columns into factors?
# factors <- colnames(FC)
# factors <- factors[-c(4,5,6,8,10,11,14:17,19)] #choose which to leave out of factorising
# FC[factors] <- lapply(FC[factors], factor)
# sapply(FC, class) #check results

#limit values to plausible ranges https://www.icp-forests.org/pdf/manual/2016/ICP_Manual_2017_01_part12.pdf
#for N
a <-
  FC %>% filter(species_name == "Picea abies") %>% filter(., between(n, 9.47, 16.68))
b <-
  FC %>% filter(species_name == "Fagus sylvatica") %>%  filter(., between(n, 20.41, 29.22))
c <-
  FC %>% filter(species_name == "Pinus sylvestris") %>%  filter(., between(n, 10.94, 20.41))
d <-
  FC %>% filter(species_name == "Quercus robur") %>% filter(., between(n, 20.31, 30.69))
e <-
  FC %>% filter(species_name == "Quercus petraea") %>%  filter(., between(n, 19.75, 29.84))
f <-
  FC %>% filter(species_name == "Quercus ilex") %>%  filter(., between(n, 11.95, 17.24))
g <-
  FC %>% filter(species_name == "Abies alba") %>%  filter(., between(n, 11.55, 16.46))
h <-
  FC %>% filter(species_name == "Pinus pinaster") %>%  filter(., between(n, 6.25, 13.71))
i <-
  FC %>% filter(
    species_name %!in% c(
      "Picea abies",
      "Fagus sylvatica",
      "Pinus sylvestris",
      "Quercus robur",
      "Quercus petraea",
      "Quercus ilex",
      "Abies alba",
      "Pinus pinaster"
    )
  )

temp <-
  full_join(a, b) %>% full_join(., c) %>% full_join(., d) %>% full_join(., e) %>%
  full_join(., f) %>% full_join(., g) %>% full_join(., h) %>% full_join(., i)
FC_bak <- FC
FC <- temp
#and for P
a <-
  FC %>% filter(species_name == "Picea abies") %>% filter(., between(p, 0.81, 2.1))
b <-
  FC %>% filter(species_name == "Fagus sylvatica") %>%  filter(., between(p, 0.89, 1.86))
c <-
  FC %>% filter(species_name == "Pinus sylvestris") %>%  filter(., between(p, 1, 2.06))
d <-
  FC %>% filter(species_name == "Quercus robur") %>% filter(., between(p, 0.97, 2.55))
e <-
  FC %>% filter(species_name == "Quercus petraea") %>%  filter(., between(p, 0.9, 1.85))
f <-
  FC %>% filter(species_name == "Quercus ilex") %>%  filter(., between(p, 0.69, 1.22))
g <-
  FC %>% filter(species_name == "Abies alba") %>%  filter(., between(p, 0.86, 2.23))
h <-
  FC %>% filter(species_name == "Pinus pinaster") %>%  filter(., between(p, 0.4, 1.38))
i <-
  FC %>% filter(
    species_name %!in% c(
      "Picea abies",
      "Fagus sylvatica",
      "Pinus sylvestris",
      "Quercus robur",
      "Quercus petraea",
      "Quercus ilex",
      "Abies alba",
      "Pinus pinaster"
    )
  )

temp <-
  full_join(a, b) %>% full_join(., c) %>% full_join(., d) %>% full_join(., e) %>%
  full_join(., f) %>% full_join(., g) %>% full_join(., h) %>% full_join(., i)

FC <- temp

#temp <- FC %>% group_by(species_name, survey_year) %>% summarise(meannp = mean(N.P))
#temp2 <- FC_bak %>% group_by(species_name, survey_year) %>% summarise(meannp = mean(N.P))

FC <- rename(FC, species_name_FC = species_name)

#average per plot
FC_means <-
  FC %>% group_by(ID,
                  survey_year,
                  country,
                  grp_tree_species,
                  species_name_FC,
                  code_leaves_type) %>%
  summarise_at(vars(n, p, s, ca, mg, k, N.P), funs(mean), na.rm = TRUE)
#include code_leaves_type?

dups_FC <-
  FC %>% group_by(ID) %>% filter(n_distinct(species_name_FC) > 1)
dups_FC_means <-
  FC_means %>% group_by(ID, survey_year) %>% filter(n_distinct(species_name_FC) >
                                                      1) %>%
  filter(code_leaves_type == 0)

#filter(n() > 1)

# FOR_FC_means <- FOR_FC %>% group_by(ID, survey_year, country, grp_tree_species, species_name, code_leaves_type) %>%
#   summarise_at(vars(n, p, s, ca, mg, k, N.P),funs(mean), na.rm = TRUE)


#plot(jitter(FC$survey_year,2), FC$N.P , cex = 0.1) 



# Import/clean IM and FOR DEP ---------------------------------------------

# DEP-Deposition dat ------------------------------------------------------

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
parameters <- c("NH4N", "NO3N", "PH")
IM_dep <- filter(IM_dep_long, ParameterCode %in% parameters) %>%
  spread(., ParameterCode, Value)

# rename columns to match FOR_dep
IM_dep <- IM_dep %>% rename(n_nh4 = NH4N, n_no3 = NO3N, ph = PH)

# average per site/year combo
IM_dep_means <- IM_dep %>% group_by(ID) %>%
  summarise_at(vars(ph, n_nh4, n_no3), funs(mean), na.rm = TRUE) #mean and sd?

# sum per site/year combo
IM_dep_sums <- IM_dep %>% group_by(ID) %>%
  summarise_at(vars(n_nh4, n_no3), funs(sum), na.rm = TRUE) #mean and sd?

# Combine IM and FOR dep --------------------------------------------------

str(FOR_dep)
str(IM_dep)
dep <- full_join(IM_dep, FOR_dep)
N_total <- dplyr::select(FOR_dep, ID, n_total) %>% drop_na %>% group_by(ID) %>% summarise_at(vars(n_total), funs(sum))
dep <- filter(dep, !is.na(survey_year))
dep <- filter(dep, !is.na(n_nh4))
dep <- filter(dep, !is.na(n_no3))

depbak <- dep

#remove some outliers- there are a very few plots with crazy high levels, probably
#near a point source (or error in recording). Remove those observations outside the
#upper limit with the value of 99th %ile

q <- quantile(dep$n_nh4, .99, na.rm = TRUE)
dep <- dep %>% filter(n_nh4 < q)

q <- quantile(dep$n_no3, .99, na.rm = TRUE)
dep <- dep %>% filter(n_no3 < q)

#depmeans_bak <- dep_means
dep_means <-
  dep %>% group_by(ID) %>% summarise_at(vars(ph, n_nh4, n_no3, n_total), funs(mean), na.rm = TRUE)
dep_means$n_total_sd <- NULL
dep_means$n_total_mean <- NULL

dep_sums<-
  dep %>% group_by(ID) %>% summarise_at(vars(n_nh4, n_no3, n_total), funs(sum), na.rm = TRUE)
dep_sums <- dep_sums %>% rename(nh4_sum = n_nh4, no3_sum = n_no3, n_total_sum = n_total)


# Import location data ----------------------------------------------------

#Import IM locations (in decimal degrees)
IM_locations <- read_csv("raw_data/Veg_IM/IM_locations.csv")
IM_locations <-
  dplyr::select(IM_locations, Code, Latitude, Longitude) %>%
  rename(ID_site = Code,
         latitude = Latitude,
         longitude = Longitude)

location <- full_join(FOR_veg_locations, IM_locations)


# VEG- Vegetation Data ----------------------------------------------------


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

#some multiple observations of same species at one site in spain. Take mean value.
IM_veg <-
  IM_veg %>% group_by(
    AREA,
    SCODE,
    SIZE,
    YYYYMM,
    SPOOL,
    SPECIES,
    PARAM,
    UNIT,
    survey_year,
    country,
    ID,
    vascular,
    latitude,
    longitude,
    species_name
  ) %>% summarise(VALUE = mean(VALUE))

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

#which(rowSums(veg, na.rm = TRUE) == 0) #which rows usm to zero

#add plot ID (only)
# veg <-
#   veg %>% mutate(plot = coalesce(as.character(survey_number), as.character(SCODE)))
# veg$plot <- as.factor(veg$plot)

#few duplicate rows, remove
veg <- distinct(veg)

#subset veg data to vascular/nonvascular
veg_bryo <- dplyr::filter(veg, vascular == "non_vasc")
veg_vasc <- dplyr::filter(veg, vascular == "vasc")
veg_lich <- dplyr::filter(veg, vascular == "lich")
#clean up, some vascular are in bryo and vice versa e.g tree seedlings recorded amongst ground layer
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
    "Pteridium aquilinum",
    "Silene italica",
    "Poa trivialis",
    "Hieracium prenanthoides",
    "Ruscus aculeatus",
    "Corydalis claviculata",
    "Symphytum tuberosum",
    "Cytisus scoparius",
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
# #Better to use BryForTrait and add extras from bryoatt?
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


#Canopy
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

saveRDS(canopy_cover,"canopy_cover.RDS")



# export data -------------------------------------------------------------


# veg_xport <- veg %>% dplyr::select(survey_year, species_name, vascular, country, ID, ID_fine, cover, latitude, longitude)
# 
# veg_mean_xport <- veg %>%  group_by(survey_year, vascular, country, ID, species_name, latitude, longitude) %>%
#   dplyr::summarise(avg = mean(cover))
# # #replace spaces with underscore
# # veg_xport$species_name <-  gsub("\\ ", "_", veg_xport$species_name)
# dep_xport <- dep %>% dplyr::select(survey_year, ID, country, n_nh4, n_no3)
# fc_xport <- FC %>% dplyr::select(survey_year, ID, country, n, p, N.P, species_name, code_leaves_type, grp_tree_species)
# #these need veg CWM calculations run first
# vasc_CWM_ell <- test_vasc
# bryo_CWM_ell <- test
# 
# write_csv(veg, "vegetation.csv")
# write_csv(dep_xport, "deposition.csv")
# write_csv(fc_xport, "foliar_chemistry.csv")
# write_csv(CWM_vasc_plot, "vasc_CWM_ell_shan_plot.csv")
# write_csv(CWM_bryo_plot, "nonvasc_CWM_ell_shan_plot.csv")
# write_csv(canopy_cover_plot, "shading.csv")
# write_csv(location, "location.csv")
# 
#write_csv(temp2,"combined_data_7_sept.csv")

# VEG CWM Ellenbergs------------------------------------------------------------
# Non vascular plants site level-----------------------------------------------------

#take average cover per species per year/site combo
bryo_means <- veg_bryo %>% ungroup() %>% 
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover)) %>%
  spread(species_name, avg)

rich_mean_bryo <- veg_bryo %>%
  group_by(ID) %>%
  summarise(rich = mean(richness))
rich_mean_bryo$ID <- as.factor(rich_mean_bryo$ID)

bryo_means_long <- veg_bryo %>%
  group_by(survey_year, country, ID, species_name) %>%
  summarise(avg = mean(cover))

bryo_means <- as_tibble(bryo_means)
bryo_means[is.na(bryo_means)] <- 0 #replace NA with zero
bryo_means <-
  bryo_means[rowSums(bryo_means[, -c(1:3)]) != 0, ] #remove rows that sum to zero
# colnames(bryo_means[,404]) #last col is NA, remove
# bryo_means[, 404] <- NULL
bryo_means_matrix <- bryo_means[, -c(1:3)] #remove factor columns

#plot(colSums(bryo_plots_matrix))
#plot(rowSums(bryo_means_matrix))

#make matrix of mean ellenberg values####
emptycols = c(which(colSums(bryo_means_matrix) == 0))
emptynames = names(emptycols)
nozeros = bryo_means_matrix %>% dplyr::select(-one_of(emptynames))
nozerosnames = colnames(nozeros)
# nozeros <- dplyr::select(nozeros, one_of(ellenbergs$species_name))
# nozeros <- filter(nozeros, rowSums(nozeros) > 0)
traits_nozeros <-
  ellenbergs[ellenbergs$species_name %in% nozerosnames,]
rownames(traits_nozeros) <- traits_nozeros$species_name
y = isc(
  veg = as.matrix(nozeros),
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
CWM_bryo <- temp

CWM_bryo$survey_year <-
  as.numeric(as.character(bryo_means$survey_year))
CWM_bryo$country <- as.factor(bryo_means$country)
CWM_bryo$ID <- as.character(bryo_means$ID)

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

div_bryo_site <- bryo_means_matrix %>% diversity(index = "shannon")
div_bryo_site <- cbind(div_bryo_site,as.character(bryo_means$ID)) %>% as_tibble() %>% 
  dplyr::select(div_bryo_site, V2) %>% rename(ID = V2)

#Simpson diversity
simp_div_bryo_site <- bryo_means_matrix %>% diversity(index = "simpson")
simp_div_bryo_site <- cbind(simp_div_bryo_site,as.character(bryo_means$ID)) %>% as_tibble() %>% rename(ID = V2) %>%
  rename( div_simp = simp_div_bryo_site)

#add diversity as column
CWM_bryo <- left_join(CWM_bryo, div_bryo_site) 
colnames(CWM_bryo)[9] <- "div"
str(CWM_bryo)
CWM_bryo$ID <- as.factor(CWM_bryo$ID)
CWM_bryo$div <- as.numeric(CWM_bryo$div)
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
div_vasc_site <- cbind(div_vasc_site,vasc_means$ID) %>% as_tibble() %>%
  select(div_vasc_site, 'vasc_means$ID') %>%  rename(ID = 'vasc_means$ID')


#add diversity as column
str(div_vasc_site)
CWM_vasc <- left_join(CWM_vasc, div_vasc_site) 
colnames(CWM_vasc)[9] <- "div"
str(CWM_vasc)
CWM_vasc$ID <- as.factor(CWM_vasc$ID)
CWM_vasc$div <- as.numeric(CWM_vasc$div)
str(CWM_vasc)

#meterological data####

library(raster)
library(sp)

#r <- getData("worldclim",var="bio",res=10)

# r <- r[[c(1,12)]]
# names(r) <- c("Temp","Prec")
# points <- cbind.data.frame(veg$ID, veg$latitude, veg$longitude)
# points <- points[!duplicated(points), ]
# values <- extract(r,points[,2:3])
# meteo <- cbind.data.frame(points,values)
# colnames(meteo) <-  c("ID", "latitude", "longitude", "temp", "prec")


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

#higher resolution?
#r <- getData("worldclim",var="bio",res=2.5) #uncomment for first run!!!XXXXXX
# r <- r[[c(1,10,11,12)]]
# names(r) <- c("mean_temp", "mean_summer_temp", "mean_winter_temp", "mean_precip")
# points <- cbind.data.frame(veg$ID, veg$latitude, veg$longitude)
# points <- points[!duplicated(points), ]
#  values <- extract(r,points[,2:3])
#  meteo <- cbind.data.frame(points,values)
#  colnames(meteo) <-  c("ID", "latitude", "longitude", "temp", "prec")

# Create a data.frame with sample site coordinates
# site <- veg$ID
# lon <- veg$longitude
# lat <- veg$latitude
# samples <- data.frame(site, lon, lat)
# points2 <- samples[!duplicated(samples), ]
# rownames(points2) <- points2$site
# points2$site <- NULL
# values2 <- extract(r,points2)
# meteo2 <- cbind.data.frame(points2,values2)
# meteo2$ID <- rownames(meteo2)
# rownames(meteo2) <- NULL
#
#saveRDS(meteo2, file = "meteo2.RDS")
meteo2 <-  readRDS("output_data/meteo2.RDS")


#
# #higher resolution
# # temp1 <- raster("wc2.0_30s_tavg_01.tif")
# # temp2 <- raster("wc2.0_30s_tavg_02.tif")
# # temp3 <- raster("wc2.0_30s_tavg_03.tif")
# # temp4 <- raster("wc2.0_30s_tavg_04.tif")
# # temp5 <- raster("wc2.0_30s_tavg_05.tif")
# # temp6 <- raster("wc2.0_30s_tavg_06.tif")
# # temp7 <- raster("wc2.0_30s_tavg_07.tif")
# # temp8 <- raster("wc2.0_30s_tavg_08.tif")
# # temp9 <- raster("wc2.0_30s_tavg_09.tif")
# # temp10 <- raster("wc2.0_30s_tavg_10.tif")
# # temp11 <- raster("wc2.0_30s_tavg_11.tif")
# # temp12 <- raster("wc2.0_30s_tavg_12.tif")
#
#


#Combine data for analysis####
#site level
VP <- CWM_vasc
VB <- CWM_bryo

VP2 <- left_join(CWM_vasc, rich_mean_vasc, by = "ID")
VB2 <- left_join(CWM_bryo, rich_mean_bryo, by = "ID")

#merge datasets# add dep data
X <- full_join(VP2, VB2, by = "ID") %>%
  left_join (dep_means, by = "ID")

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
X2 <- left_join(X, meteo2, by = "ID")
X2$lon <- NULL
X2$lat <- NULL
X2$mean_summer_temp <- X2$mean_summer_temp / 10
X2$mean_winter_temp <- X2$mean_winter_temp / 10
X2$mean_temp <- X2$mean_temp / 10
X <- X2
rm(X2)

FC_means2 <- FC_means
FC_means2$grp_tree_species <- NULL
#FC means has data for current and older needles and in some cases multiple species
#per plot
#To avoid many duplicate rows in a combined dataframe, make version with only
#current year and mean values across all species for a by plot mean value (without
#species column) and use full version where e.g. species specific reponses are required
FC_means_for_df <-
  filter(ungroup(FC_means), code_leaves_type == 0) %>% dplyr::select(-c(species_name_FC, s, ca , mg, k)) %>%
  group_by(ID) %>% summarise_at(vars(n:N.P), mean)

#add FC data so we can relate canopy N.P ratio to things#
temp <- X
temp2 <- left_join(temp, FC_means_for_df, by = "ID")

temp2$country.y <- NULL
temp2$ph <- NULL

#combined richness for vasc and non
#temp2 <- mutate(temp2, rich.combi = rich.x + rich.y)

#add conifer/decid group data from veg_tree_groups
veg_tree_groups <-
  full_join(FOR_veg_tree_groups, IM_veg_tree_groups)

temp2 <-
  left_join(temp2,
            dplyr::select(veg_tree_groups, species_name, grp_tree_species, ID),
            by = "ID")

#some gaps that can be filled with data from FC
# temp2$grp_tree_species2 <- ifelse(is.na(temp2$grp_tree_species), temp2$species_name_FC, temp2$grp_tree_species)
# temp2$grp_tree_species3 <- ifelse(temp2$grp_tree_species2 %in% decid_FOR_veg_trees,
#                                   'broadleaves', ifelse(temp2$grp_tree_species2 %in% conif_FOR_veg_trees,
#                                   'conifers', temp2$grp_tree_species2))
# temp2$grp_tree_species2 <- NULL
# temp2 <- rename(temp2, grp_tree_species_bak = grp_tree_species)
# temp2 <- rename(temp2, grp_tree_species = grp_tree_species3)

# joint_grp <- coalesce(temp2$grp_tree_species.x, temp2$grp_tree_species.y)
# temp2$joint_grp <- joint_grp
# temp2$grp_tree_species.x <- NULL
# temp2$grp_tree_species.y <- NULL
# temp2 <- rename(temp2, grp_tree_species = joint_grp)

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
#Split the data frame up - you could probably do this in one go but this makes things a bit more obvious:
query = data[is.na(data$grp_tree_species), ]
tree_groups = data[!is.na(data$grp_tree_species), ]
library(FNN)
neighs = get.knnx(tree_groups[, c("latitude", "longitude")], query[, c("latitude", "longitude")], k =
                    1)
#Now insert the replacement colours directly into the data dataframe:

data[is.na(data$grp_tree_species), "grp_tree_species"] = tree_groups$grp_tree_species[neighs$nn.index]
# plot(data$longitude,
#      data$latitude,
#      col = data$grp_tree_species,
#      pch = 19)

extra_grp <- dplyr::select(data, grp_tree_species, ID)
temp2$grp_tree_species <- extra_grp$grp_tree_species


#
# joint_spp <- coalesce(temp2$species_name_FC, temp2$species_name)
# joint_spp2 <- coalesce(temp2$species_name, temp2$species_name_FC)
# temp2$joint_spp <- joint_spp
# temp2$species_name_FC <- NULL
# temp2$species_name <- NULL
# temp2 <- rename(temp2, species_name = joint_spp)
#
# #make combo dep variable
# temp2 <- temp2 %>% mutate(N_dep = n_nh4 + n_no3)

#make ID_siteplot
#temp2 <- mutate(temp2, ID_siteplot = substr(temp2$ID_fine, 6, 20))

#fill missing code_leaves with zeros
#temp2$code_leaves_type <- replace_na(temp2$code_leaves_type, 0)

temp2 <- distinct(temp2)

#ICP IM generally only has N and P, drop_na() will lose them. Instead drop these
#columns (can come back if we decide to use these variables)
#temp2 <- dplyr::select(temp2, -c(mg, s, ca, k))

#temp2 <- temp2 %>% rename(country = country.x)


#temp2$code_leaves_type <- as.factor(temp2$code_leaves_type)
temp2$grp_tree_species <- as.factor(temp2$grp_tree_species)
temp2$species_name <- as.factor(temp2$species_name)
#temp2$species_name_FC <- as.factor(temp2$species_name_FC)
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
#species data
tem <- bryo_means
tem <- filter(tem, ID %in% temp2$ID)#temp2?
tem$survey_year <- as.numeric(as.character(tem$survey_year))
tem_matrix <- rename(tem, rowname = ID) %>% column_to_rownames() %>% dplyr::select(.,-c(1:2))

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
#filter(bryo_plots_long, species_name %in% notrait) %>% group_by(ID,species_name) %>% 
#count() %>% arrange(desc(n)) %>% View
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
y2 <-  dbFD(traits_nozeros2,nozeros2, corr = "none", calc.FRic = FALSE)

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

# simple.FD <- left_join(simple, FD_site, by = "ID")
# simple.FD$ID <- as.factor(simple.FD$ID)
# 
# saveRDS(simple.FD, file = "simple.FD.RDS")
# 

#Subsets of temp2####
#Drop NAs
temp2_complete <- drop_na(temp2)

#temp2_range####
#specify date ranges and find which (and how many) plots are in both data ranges####
t1 <- temp2 %>% filter(survey_year <= 2000)
t2 <- temp2 %>% filter(survey_year >= 2005)
length(unique(t1$ID_site))#678
length(unique(t2$ID_site))#646
length(unique(intersect(t1$ID_site, t2$ID_site)))#334
ID_list <- unique(intersect(t1$ID_site, t2$ID_site))#334

#filter full data to only plots found in both date ranges
temp2_range <- filter(temp2, ID_site %in% ID_list)#7814
#temp2_range$survey_year <- as.factor(temp2_range$survey_year)
#Scale numeric variables####
non_nums <- select_if(temp2_range, negate(is.numeric))
nums <-
  as.data.frame(select_if(temp2_range, is.numeric) %>% scale())
temp2_range_scaled <- bind_cols(non_nums, nums)

#subset both date ranges for comparison
temp2_range.t1 <-
  temp2_range %>% filter(survey_year <= 2000)#2001 observations (ID_fine)
temp2_range.t2 <-
  temp2_range %>% filter(survey_year >= 2005)#1161 observations (ID_fine)

#only sites that have N dep data
temp2N <- filter(temp2, n_no3 > 0)
#specify date ranges and find which (and how many) plots are in both data ranges####
t1N <- temp2N %>% filter(survey_year <= 2000)
t2N <- temp2N %>% filter(survey_year >= 2005)
length(unique(t1N$ID_site))#292
length(unique(t2N$ID_site))#350
length(unique(intersect(t1N$ID_site, t2N$ID_site)))#201
ID_list <- unique(intersect(t1N$ID_site, t2N$ID_site))#2119

#filter full data to only plots found in both date ranges
temp2N_range <- filter(temp2N, ID_site %in% ID_list)#1114
t1N <- temp2N_range %>% filter(survey_year <= 2000)
t2N <- temp2N_range %>% filter(survey_year >= 2005)

#filter to both bryo and dep data available
temp2N_range_comp <- filter(temp2N_range, N.y > 0)


#separate "environmental variables" dataframe for analyses####
env <- dplyr::select(temp2, country,ID_site, ID, n_nh4, n_no3,
                     latitude, longitude, code_country, sum_canopy, mean_temp,
                     mean_precip, n, p, N.P, survey_year, 
                     grp_tree_species)

env %>% group_by(ID) %>% filter(n() > 1)#no duplictates

env_means <- env
env_means$plot.x <- NULL
env_means$ID_siteplot <- NULL
env_means$species_name <- NULL
env_means <- env_means %>% dplyr::group_by(ID) %>% summarise_if(is.numeric, mean) #put back country names and grp_tree?
env_means <- env_means[!duplicated(env_means), ]
env_means_scaled <- scale(dplyr::select(env_means, -c("ID")))  %>% as.data.frame() %>% mutate(ID = env_means$ID)

bryo_env_means <- env_means %>% filter(ID %in% bryo_means$ID)

#Export data####
saveRDS(temp2, file = "temp2.RDS")
saveRDS(temp2_range, file = "temp2_range.RDS")
saveRDS(temp2N_range, file = "temp2N_range.RDS")
saveRDS(temp2N_range_comp, file = "temp2_range_comp.RDS")
saveRDS(temp2_range.t1, file = "temp2_range.t1.RDS")
saveRDS(temp2_range.t2, file = "temp2_range.t2.RDS")
saveRDS(t1N, file = "t1N.RDS")
saveRDS(t2N, file = "t2N.RDS")
saveRDS(simple, file = "simple.RDS")
saveRDS(simple_scaled, file = "simple_scaled.RDS")

saveRDS(veg, file = "veg.RDS")
saveRDS(veg_vasc, file = "veg_vasc.RDS")
saveRDS(veg_bryo, file = "veg_bryo.RDS")
saveRDS(ellenbergs, file = "ellenbergs.RDS")
#saveRDS(BryForTrait, file = "BryForTrait.RDS")
saveRDS(bryophyte_func, file = "bryophyte_func.RDS")
saveRDS(dep, file = "dep.RDS")
saveRDS(dep_means, file = "dep_means.RDS")
saveRDS(FC, file = "FC.RDS")
saveRDS(FC_means, file = "FC_means.RDS")
saveRDS(env, file = "env.RDS")
saveRDS(env_means, file = "env_means.RDS")
saveRDS(env_means_scaled, file = "env_means_scaled.RDS")
# saveRDS(vasc_plots_wide, file = "vasc_plots_wide.RDS")
# saveRDS(vasc_plots_matrix, file = "vasc_plots_matrix.RDS")
# saveRDS(bryo_plots_wide, file = "bryo_plots_wide.RDS")
# saveRDS(bryo_plots_matrix, file = "bryo_plots_matrix.RDS")
#saveRDS(meteo2, file = "meteo2.RDS")
temp2 <- rename(temp2, div.y = div)
#Generate data3, per site means####
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
    n_nh4,
    n_no3,
    N.P,
    latitude,
    longitude,
    survey_year,
    mean_precip, 
    mean_temp,
    ID_site,
    ID_siteonly,
    ID,
    country
  ) %>% 
  left_join(., FD_site, by = "ID") %>% 
  filter(., latitude >=45.0) %>% filter (longitude > -2.5) %>% filter(longitude < 30) %>% 
  drop_na(N.y) %>% drop_na(n_nh4)
data3 <- left_join(data3, canopy_cover)
data3 <- ungroup(data3)
data3 <- rename(data3, sum_canopy = canopy)
#data3 <- left_join(data3, p_mean)
data3$ID <- as.factor(data3$ID)
data3 <- dplyr::select(data3, country, ID_site, ID_siteonly, ID, everything())

data3 <- filter(data3, L.y > 0)
data3 <- filter(data3, N.y > 0)
data3 <- filter(data3, L.x > 0)
data3 <- filter(data3,n_no3<6) #one site with outlier high no3 left after 99th percentile, remove
#keep only sites with at least 2 species of bryophyte present
data3 <- filter(data3, rich.y >= 2)

#create index of year
data3$year.i <- I(data3$survey_year - 1994)

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

#add simpson
data3 <- left_join(data3, simp_div_bryo_site, by = "ID")
data3$div_simp <- as.numeric(data3$div_simp)

#add n_total
data3_nt <- left_join(data3, N_total, by ="ID")

#library(VIM)
#matrixplot(data3, interactive = F, sortby = "sum_canopy") #VIM

saveRDS(data3, file = "data3.RDS")

#drop NAs
dat <- data3
#keep only columns we need
datc <- select(dat,country,ID,ID_site,ID_siteonly,N.y,div.y,rich.y,RaoQ,RaoQ_m,n_nh4,n_no3,
               sum_canopy,L.x,mean_temp,mean_precip,latitude,longitude,survey_year,tree)
#drop na 
datc <- drop_na(datc)

#drop data above 99th percentile in each explanatory variable
not_used <- select(datc, -c(L.x, mean_precip, mean_temp, n_nh4, n_no3, sum_canopy))
p99 <- datc %>% select(L.x, mean_precip, mean_temp, n_nh4, n_no3, sum_canopy) %>% 
  gather(key, value) %>%
  group_by(key) %>%
  mutate(q99 = quantile(value, 0.99), row = row_number()) %>% 
  filter(value <= q99) %>%
  select(-q99) %>%
  spread(key, value) %>%
  select(-row)
new <- cbind(not_used, p99)  
newd <- drop_na(new) 

#check
a <- select(datc, ID, N.y, mean_precip, mean_temp, L.x)
b <- select(newd, ID, N.y, mean_precip, mean_temp, L.x)
c <- left_join(a,b, by= "ID")

datd <- newd

#select/subset based on data quality####
#data quality
qual.s <- veg_bryo %>% group_by(country,ID) %>% summarise(Unique_spp = n_distinct(species_name))
qual.p <- veg_bryo %>% group_by(ID_fine) %>% summarise(Unique_spp = n_distinct(species_name))
hist(qual.s$Unique_spp)
hist(qual.p$Unique_spp)
qual.c <-veg_bryo %>% group_by(country) %>% summarise(Unique_spp = n_distinct(species_name))
hist(qual.c$Unique_spp)
qual2 <- veg_bryo %>% group_by(ID) %>% mutate(rich_s = mean(richness)) %>% ungroup() %>% 
  group_by(country) %>% mutate(rich_c = mean(richness))
hist(qual2$richness)
hist(qual2$rich_s)
table(qual2$country)


df.sum <- qual2 %>% ungroup() %>% group_by(country) %>% 
  select(country, rich_s, rich_c, richness) %>% 
  summarise_each(funs(min = min, 
                      max = max,
                      mean = mean, 
                      sd = sd)) %>% 
  select(country, rich_c_mean,
         rich_s_min, rich_s_max, rich_s_mean,
         richness_min, richness_max, richness_mean)

rich.country <-dat1 %>% group_by(country) %>% summarise(rich = mean(rich.y))

#unique species by country
Dsum<-veg_bryo %>% group_by(country) %>% summarise(unique_spp = n_distinct(species_name)) %>% 
  arrange(unique_spp)%>% mutate(row = row_number())
ggplot(Dsum,aes(x=reorder(country,row),y=unique_spp))+geom_point()+geom_hline(yintercept=20)+xlab("country")+
  theme(axis.text.x=element_text(angle=90, hjust=1))
#log
ggplot(Dsum,aes(x=reorder(country,row),y=log(unique_spp)))+geom_point()+geom_hline(yintercept=log(20))+xlab("country")+
  theme(axis.text.x=element_text(angle=90, hjust=1))

#try removing countries that record less than 20 species
dat1.qual.country <- filter(datc, country %!in% c("Estonia", "Luxembourg", "Romania", "Czech Republic"))



#STOP####

#rm(list = ls())


#
# test <- left_join(bryo_means, bryo_env_means, by = "ID") #join and divide to get order of rows right
# b.sp.ord <- dplyr::select(test, -c(ID, n_nh4, n_no3, latitude, longitude,country, code_country, sum_canopy, mean_temp,
#                                     mean_precip, n, p, N.P, survey_year.x))
# b.env.ord <- dplyr::select(test, ID, n_nh4, n_no3, latitude, longitude,country, code_country, sum_canopy, mean_temp,
#                    mean_precip, n, p, N.P, survey_year.x)
#check all NA rows at start??
# test1 <- cca(b.sp.ord,na.action = na.exclude)
# test2 <- cca(b.sp.ord, b.env.ord, na.action = na.exclude)
# test3 <- rda(b.sp.ord, b.env.ord, na.action = na.exclude)

#ordiplot(test1, type = "text")
##sites 147 (1996_9_1) and 355 (1999_50_13) are crazy outliers
# test <- bryo_means[-c(147,355),] #remove
# test.matrix <- test[,-c(1:3)]
# test1 <- cca(test.matrix)
# ordiplot(test1, type = "text")

#outliers test####

#outliers <- outlier(bryo_means_matrix, thresh = 0.2, y = 0.5)
#outliers.v <- outlier(veg_means_matrix, thresh = 0.2, y = 0.5)
# test <- cca(outliers$new.data)
# ordiplot(test, type = "text")
# 
# de_bryo_means <- filter(bryo_means, country == "Germany")
# de_bryo_means_matrix <- de_bryo_means[,-c(1:3)]
# outliers_de <- outlier(de_bryo_means_matrix, thresh = 0.2, y = 0.5)
# test1 <- metaMDSdist(outliers_de$new.data)
# ordiplot(test1, type = "text")
# 
# testdat <- dat1[,c(1:5)]
# testdat <- drop_na(testdat)
# 
# test <- decorana(testdat)
# ordiplot(test, type = "text")
# 
# test <- decorana(dat1[,c(12:16)])
# ordiplot(test, type = "text")

# #Missing data#####
# #missing data check
# summary(env)
# pMiss <- function(x){sum(is.na(x))/length(x)*100}
# apply(env,2,pMiss)
# apply(env,1,pMiss)
# 
# library(mice)
# md.pattern(env)
# 
# library(VIM)
# aggr_plot <- aggr(env, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
#                   labels=names(env), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
# 
# #compare distributions of complete and missing to see if missing is randomly distributed (MCAR)
# marginplot(env[c(6,5)])

#impute missing data
# tempData <- mice(env,m=5,maxit=20,meth='pmm',seed=500)
# summary(tempData)
# completedEnv <- complete(tempData,1)
# 
# densityplot(tempData)
# xyplot(tempData,n_nh4 ~ n_no3+mean_temp+sum_canopy,pch=18,cex=1)
# stripplot(tempData, pch = 20, cex = 1.2)
# 
# impute.env <- completedData
# impute.env.scaled <- scale(completedData)