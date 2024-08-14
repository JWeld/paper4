# Skript för att läsa VG data från internationella IM-databasen

setwd("C:/Users/grandin/Documents/Ulf IVM/IM/.Internationella ordförandeskapet/Strategi/Ellenberg/HelaIM")

library(vegan)
library(readxl)
library(tidyverse)
library(ggplot2)
library(rkt)
library(vegdata)


# Read data ----
df<-as.data.frame(read_xlsx("VG.xlsx", sheet = "Blad1",col_names = TRUE))
#names(df)

# skapa ID för varje provyta
df$Area_Stn<-paste(df$AreaCode,df$StationCode, sep = "_")

# ta bort alla med färre än 4 år
df<-df %>% 
filter(Area_Stn %in% c("LT01_100", "LT01_102", "LT03_100", "NO01_1", "NO02_1", "SE04_2", "SE14_1", "SE15_1", "SE15_2", "SE16_1", "SE16_2"))

# ta bort allt som inte är COVE_B eller COVE_F
df<-df %>% 
  filter(Parameter %in% c("COVE_B", "COVE_F"))

#Lägg till år
df$YEAR<-as.factor(substr(df$YearMonth,1,4))
df$YEAR<-as.character(df$YEAR)


# Ta bort rader utan angiven art
df<-df[!is.na(df$Medium),]

# ta bort år 2015 från LT01_100 och LT01_102
df<-df[-which(df$Area_Stn == "LT01_100" & df$YEAR == 2015),]
df<-df[-which(df$Area_Stn == "LT03_100" & df$YEAR == 2015),]

# ta bort år 2022 från SE14_1
df<-df[-which(df$Area_Stn == "SE14_1" & df$YEAR == 2022),]


## Dela upp på respektive område
df1<-split(df, df$Area_Stn, drop = TRUE)


# Transform long to wide
lapply(seq_along(df1), function(i) assign(paste0(names(df1)[i],"w"), df1[[i]] %>% 
         pivot_wider(names_from = Medium, values_from = Value, values_fill = 0), envir = .GlobalEnv))



# Ellenberg-värden
Ell<-read_excel("C:/Users/grandin/Documents/Ulf IVM/IM/Rapportering/Årsrapporter/R för I2/Ellenberg.xlsx", sheet = 1)



# LT01_100w----
#names(LT01_100w)
område<-"LT01_100"
year<- levels(as.factor(LT01_100w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(LT01_100w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  LT01_100w_y<-as.data.frame(LT01_100w[LT01_100w$YEAR == year[i], 26:(ncol(LT01_100w))]) 
  
  Ell_R<-cwm(veg = LT01_100w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-paste("EllR",område, sep = "_")

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-paste("mk",område, sep = "_")

# diversitet
H<-as.data.frame(diversity(LT01_100w[, 26:(ncol(LT01_100w))], index = "shannon", groups = LT01_100w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-paste("H",område, sep = "_")

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-paste("mkH",område, sep = "_")



# LT01_102w----
names(LT01_102w)
område<-"LT01_102"
year<- levels(as.factor(LT01_102w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(LT01_102w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  LT01_102w_y<-as.data.frame(LT01_102w[LT01_102w$YEAR == year[i], 26:(ncol(LT01_102w))]) 
  
  Ell_R<-cwm(veg = LT01_102w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(LT01_102w[, 26:(ncol(LT01_102w))], index = "shannon", groups = LT01_102w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# LT03_100w----
names(LT03_100w)
område<-"LT03_100"
year<- levels(as.factor(LT03_100w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(LT03_100w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  LT03_100w_y<-as.data.frame(LT03_100w[LT03_100w$YEAR == year[i], 26:(ncol(LT03_100w))]) 
  
  Ell_R<-cwm(veg = LT03_100w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(LT03_100w[, 26:(ncol(LT03_100w))], index = "shannon", groups = LT03_100w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# NO01_1w----
names(NO01_1w)
område<-"NO01_1"
year<- levels(as.factor(NO01_1w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(NO01_1w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  NO01_1w_y<-as.data.frame(NO01_1w[NO01_1w$YEAR == year[i], 26:(ncol(NO01_1w))]) 
  
  Ell_R<-cwm(veg = NO01_1w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(NO01_1w[, 26:(ncol(NO01_1w))], index = "shannon", groups = NO01_1w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# NO02_1w----
names(NO02_1w)
område<-"NO02_1"
year<- levels(as.factor(NO02_1w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(NO02_1w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  NO02_1w_y<-as.data.frame(NO02_1w[NO02_1w$YEAR == year[i], 26:(ncol(NO02_1w))]) 
  
  Ell_R<-cwm(veg = NO02_1w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(NO02_1w[, 26:(ncol(NO02_1w))], index = "shannon", groups = NO02_1w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# SE04_2w----
names(SE04_2w)
område<-"SE04_2"
year<- levels(as.factor(SE04_2w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(SE04_2w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")


for (i in 1:length(year)) {
  SE04_2w_y<-as.data.frame(SE04_2w[SE04_2w$YEAR == year[i], 26:(ncol(SE04_2w))]) 
  
  Ell_R<-cwm(veg = SE04_2w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(SE04_2w[, 26:(ncol(SE04_2w))], index = "shannon", groups = SE04_2w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# SE14_1w----
names(SE14_1w)
område<-"SE14_1"
year<- levels(as.factor(SE14_1w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(SE14_1w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  SE14_1w_y<-as.data.frame(SE14_1w[SE14_1w$YEAR == year[i], 26:(ncol(SE14_1w))]) 
  
  Ell_R<-cwm(veg = SE14_1w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(SE14_1w[, 26:(ncol(SE14_1w))], index = "shannon", groups = SE14_1w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# SE15_1w----
names(SE15_1w)
område<-"SE15_1"
year<- levels(as.factor(SE15_1w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(SE15_1w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  SE15_1w_y<-as.data.frame(SE15_1w[SE15_1w$YEAR == year[i], 26:(ncol(SE15_1w))]) 
  
  Ell_R<-cwm(veg = SE15_1w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(SE15_1w[, 26:(ncol(SE15_1w))], index = "shannon", groups = SE15_1w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# SE15_2w----
names(SE15_2w)
område<-"SE15_2"
year<- levels(as.factor(SE15_2w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(SE15_2w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  SE15_2w_y<-as.data.frame(SE15_2w[SE15_2w$YEAR == year[i], 26:(ncol(SE15_2w))]) 
  
  Ell_R<-cwm(veg = SE15_2w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(SE15_2w[, 26:(ncol(SE15_2w))], index = "shannon", groups = SE15_2w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# SE16_1w----
names(SE16_1w)
område<-"SE16_1"
year<- levels(as.factor(SE16_1w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(SE16_1w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  SE16_1w_y<-as.data.frame(SE16_1w[SE16_1w$YEAR == year[i], 26:(ncol(SE16_1w))]) 
  
  Ell_R<-cwm(veg = SE16_1w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(SE16_1w[, 26:(ncol(SE16_1w))], index = "shannon", groups = SE16_1w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))


# SE16_2w----
names(SE16_2w)
område<-"SE16_2"
year<- levels(as.factor(SE16_2w$YEAR))
EllR<- matrix(NA, nrow=nlevels(as.factor(SE16_2w$YEAR)), ncol=4)
colnames(EllR)<-c("mean", "sd", "Year", "Site")

for (i in 1:length(year)) {
  SE16_2w_y<-as.data.frame(SE16_2w[SE16_2w$YEAR == year[i], 26:(ncol(SE16_2w))]) 
  
  Ell_R<-cwm(veg = SE16_2w_y, trait.db = Ell, ivname = "R", method = "mean", keyname = "Fältskikt")
  Ell_R[!(apply(Ell_R, 1, function(y) any(y == 0))),] # Ta bort nollor
  
  EllR[i, 1] <- mean(Ell_R)
  EllR[i, 2] <- sd(Ell_R)
  EllR[i, 3] <- year[i]
  EllR[i, 4] <- område
}
assign(paste("EllR",område, sep = "_"), as.data.frame(EllR))
EllR_namn<-append(EllR_namn, paste("EllR",område, sep = "_"))

# Mann-Kendall
EllR<-as.data.frame(EllR)
EllR[,1:3]<- EllR[,1:3] %>% mutate_if(is.character, as.numeric)
mk<-rkt(EllR$Year, EllR$mean) 
mk
mk<- unlist(mk) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mk",område, sep = "_"), as.data.frame(mk)))
mk_namn<-append(mk_namn, paste("mk",område, sep = "_"))

# diversitet
H<-as.data.frame(diversity(SE16_2w[, 26:(ncol(SE16_2w))], index = "shannon", groups = SE16_2w$Year))
colnames(H)[1]<-"Shannon"
H$Year<-as.numeric(rownames(H))
H$Site<-område
H
(assign(paste("H",område, sep = "_"), as.data.frame(H)))
H_namn<-append(H_namn, paste("H",område, sep = "_"))

mkH<-rkt(H$Year, H$Shannon) 
mkH
mkH<- unlist(mkH) %>% t() %>%  as.data.frame() %>% 
  select(sl, B) %>%  mutate(Site = område)
(assign(paste("mkH",område, sep = "_"), as.data.frame(mkH)))
mkH_namn<-append(mkH_namn, paste("mkH",område, sep = "_"))



############################################################
# Slå samman Mann-Kendall på EllR till en tabell ---------
mkR_alla<-do.call(rbind, mget(mk_namn))
#colnames(mkR_alla)<-c("p(sl)","B","Site")
mkR_alla %>% 
  filter(sl<0.05)

# Slå samman diversitet på EllR till en tabell ---------
H_alla<-do.call(rbind, mget(H_namn))
H_alla
##colnames(mkR_alla)<-c("p(sl)","B","Site")

# Slå samman Mann-Kendall på Shannon till en tabell ---------
mkH_alla<-do.call(rbind, mget(mkH_namn))
#colnames(mkR_alla)<-c("p(sl)","B","Site")
mkH_alla %>% 
  filter(sl<0.05)

# Alla Ell_R i samma plot -----------
EllR_alla<-do.call(rbind, mget(EllR_namn))
EllR_alla[,1:3]<- EllR_alla[,1:3] %>% mutate_if(is.character, as.numeric)
EllR_alla$Year<-as.integer(EllR_alla$Year)
#str(EllR_alla)

# Plot Shannon
ggplot(H_alla, aes(x=Year, y=Shannon, group = Site)) +
  scale_x_continuous(breaks=seq(min(H_alla$Year),max(H_alla$Year),2)) +
  scale_shape_manual(values=1:nlevels(as.factor(H_alla$Site))) +
  geom_point(aes(shape=Site, color=Site),size = 1) +
  geom_line(aes(color=Site),linewidth = 0.3) +
  ylab("Shanon diversity index") +
  stat_smooth(aes(colour = Site),method=lm, se= FALSE, linewidth = 0.5,
              data = H_alla[which(H_alla$Site %in% c("LT01_100","SE04_2","SE14_1","SE15_1")),]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Shannon index u FI.png", width = 15,  height = 10,  units = "cm")



#Plot Ellenberg ----

ggplot(EllR_alla, aes(x=Year, y=mean, group = Site)) +
  scale_x_continuous(breaks=seq(min(EllR_alla$Year),max(EllR_alla$Year),2)) +
  scale_shape_manual(values=1:nlevels(as.factor(EllR_alla$Site))) +
  geom_point(aes(shape=Site, color=Site),size = 1) +
  geom_line(aes(color=Site),linewidth = 0.3) +
  ylab("Ellenberg R index") +
  stat_smooth(aes(colour = Site),method=lm, se= FALSE, linewidth = 0.5,data = EllR_alla[which(EllR_alla$Site %in% c("SE14_1","SE15_1","SE15_2","SE16_1","SE16_2")),]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Ellenberg R index u FI.png", width = 15,  height = 10,  units = "cm")

