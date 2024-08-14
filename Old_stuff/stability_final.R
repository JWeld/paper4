source("file:///C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Scripts/SpUseDD.R")

library(reshape)
library(vegan)
library(tidyverse)
library(magrittr)
library(zoo)

lakedetails <- read.csv("file:///C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Data Files/LakeInfo109.csv", header=TRUE)
#setwd("C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/LN_fromEmma/To Hannah")
#map_map <- readOGR("limnorrl.shp")
#map_wgs84 <- spTransform(map_map, CRS("+proj=longlat +datum=WGS84"))
#centroid <- gCentroid(map_wgs84)
#mapImageData1 <- get_stamenmap(bbox = c(left=8,bottom=69,right=25,top=55), color = "color", maptype = "terrain", zoom=4)
#polys = fortify(map_wgs84)


Full <- read.csv("file:///C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Data Files/BiologicalData/MZBOct2018.csv", sep=",", header=TRUE)
Full$Year <- as.factor(Full$ProvtagningsÂr)
Full$Month <- as.factor(Full$ProvtagningsmÂnad)
Full$Lake <- Full$Stationsnamn
Full2 <- Full %>%
  filter(Month %in% c(8, 9, 10, 11), Vattenzon.P.L.SP. %in% "Litoral", Delprogram %in% c("Trendsjˆar", "NM÷"), Projekt %in% c("Bottenfauna i trendsjˆar","Sjˆar trendstationer"))


TM <- cast(Full2, Lake + Year + Stationskoordinat.N.X + Stationskoordinat.E.Y ~ Taxonnamn, value="Medelantal.per.prov", fun.aggregate=sum)
TM2 <- as.data.frame(TM)
SY <- table(TM2$Stationskoordinat.N.X, factor(TM2$Year))
SY2 <- as.data.frame.matrix(SY)
NYS <- SY2[rowSums(SY2)>9,]
AYsmp <- as.factor(rownames(NYS)) 
AYSmp2 <- TM2[TM2$Stationskoordinat.N.X %in% AYsmp,]
AYSmp2$Lake <- factor(AYSmp2$Lake)
SpUse <- append(SpUse, c("Lake", "Year", "Stationskoordinat.N.X", "Stationskoordinat.E.Y"))
AYSmp3 <- AYSmp2[,which(names(AYSmp2) %in% SpUse)]
AYSmp3[is.na(AYSmp3)] <- 0
rownames(AYSmp3) <- paste(AYSmp3$Lake, AYSmp3$Year)
AYSmp3$LakeYear <- paste(AYSmp3$Lake, AYSmp3$Year)
AYSmp3$Richness <- specnumber(AYSmp3[,5:233])
AYSmp3$id  <- 1:nrow(AYSmp3)

AYSmp4 <- AYSmp3 %>% 
  group_by(Lake, Year, Stationskoordinat.E.Y, Stationskoordinat.E.Y) %>% 
  select_if(function(x)(sum(x!=0, na.rm = T) > 106))

#AYSmp4 <- AYSmp3[,c(TRUE, TRUE, TRUE, TRUE, colSums(AYSmp3[,5:233]!=0) > 106)] 
AYSmp4$LakeYear <- paste(AYSmp4$Lake, AYSmp4$Year)
AYSmp4$Richness <- specnumber(AYSmp4[,5:123])
AYSmp4$id  <- 1:nrow(AYSmp4)

AYSmp4 <- AYSmp4[AYSmp4$Lake!="Storvindeln",]

#names(AYSmp4)[c(22,74,59,24,21,39,89)] <- c("C. luctuosa", "L. vespertina", "H. stagnalis", "C. luteolum","C. horaria", "D. vulneratus", "O. testacea")

morelakes <- merge(AYSmp4, lakedetails, by="Lake")

write.csv2(morelakes, "file:///C:/Users/hhen0001/Documents/QuantumPaper/Data/InvertData.csv")


MLsort <- morelakes[order(morelakes$id), ]
MLsort$GLlakes <- c(rep("FALSE", 123), rep("GL1", 23), rep("FALSE", 696), rep("GL2", 11), rep("FALSE",1279))

##CA for inertia
MZBca <- cca(downweight(AYSmp4[,5:123]))

#Hellinger transformed data
Hellinger <- decostand(AYSmp4[,5:123], "hellinger") 

#down <-downweight(MLsort[,5:233], fraction = 5)
#names(AYSmp4)[c(18,70,55,20,17,23,31,35,85,92)] <- c("C. luctuosa", "L. vespertina", "H. stagnalis", "C. luteolum","C. horaria", "Cladotanytarsus", "Cryptochironomus", "D. vulneratus", "O. testacea","Parakiefferiella")

DCAall <- decorana(AYSmp4[,5:123], iweigh=1)

PCAall <- prcomp(Hellinger)

DCAall$evals/5.603096
SiteSc <- as.data.frame(scores(DCAall, choices=c(1,2)))


#strings <- strsplit(names(AYSmp4[,5:123]), " ")
#genus <- lapply(strings, function(x) x[1])
#families <- tax_name(genus,get="family", db = 'itis')
#order <- tax_name(genus,get="order", db = 'itis')
#taxainfo <- cbind(names(AYSmp4[,5:123]), families[3], order[3])

## Matching environmental data to species data from here
CP3 <- read.csv("file:///C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Data Files/PhysChemData/ChemDataOct2018.csv", sep=",", header=TRUE)
CP3$Year <- as.factor(CP3$ProvtagningsÂr)
CP3$Month <- as.factor(CP3$ProvtagningsmÂnad)
CP3$Lake <- CP3$Stationsnamn
CP4 <- CP3 %>%
  filter(Year %in% c("1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017"),
         Month %in% c(8,9,10,11),
         Stationskoordinat.N.X %in% AYSmp4$Stationskoordinat.N.X)

co.var <- function(x) (100*sd(x)/mean(x))
HHenv <- CP3 %>%
  filter(Year %in% c("1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017"),
         Stationskoordinat.N.X %in% AYSmp4$Stationskoordinat.N.X)
HHenv$Measurements <- as.numeric(levels(HHenv$V‰rde.Koncentration))[HHenv$V‰rde.Koncentration]
HHenvfix <- cast(HHenv, Lake + Year + Month + Stationskoordinat.N.X + Stationskoordinat.E.Y ~ Parameter, value="Measurements", fun.aggregate=mean)
HHenvfix2<- HHenvfix[complete.cases(HHenvfix[c(6:9)]),]

CVs <- cbind(aggregate(`Alk/Acid` ~ Lake + Year, data = HHenvfix2, co.var),
             aggregate(TOC ~ Lake + Year, data = HHenvfix2, co.var)[3],
             aggregate(`Tot-P` ~ Lake + Year, data = HHenvfix2, co.var)[3],
             aggregate(Vattentemperatur ~ Lake + Year, data = HHenvfix2, co.var)[3])

CVs2 <- cbind(aggregate(`Alk/Acid` ~ Lake , data = HHenvfix2, co.var),
              aggregate(TOC ~ Lake, data = HHenvfix2, co.var)[2],
              aggregate(`Tot-P` ~ Lake, data = HHenvfix2, co.var)[2],
              aggregate(Vattentemperatur ~ Lake, data = HHenvfix2, co.var)[2])

write.csv2(HHenvfix2, "file:///C:/Users/hhen0001/Documents/Stability Course/Paper/Data/WC_allMO.csv")
write.csv2(CVs, "file:///C:/Users/hhen0001/Documents/Stability Course/Paper/Data/WC_CVs.csv")
write.csv2(CVs2, "file:///C:/Users/hhen0001/Documents/Stability Course/Paper/Data/WC_CVs2.csv")
write.table(AYSmp4[,5:123], "file:///C:/Users/hhen0001/Documents/Stability Course/Paper/Data/OTUs.csv", sep=",")

CP4$Measurements <- as.numeric(levels(CP4$V‰rde.Koncentration))[CP4$V‰rde.Koncentration]
CP5 <- cast(CP4, Lake + Year + Stationskoordinat.N.X + Stationskoordinat.E.Y ~ Parameter, value="Measurements", fun.aggregate=mean)
CP5$LakeYear <- paste(CP5$Lake, CP5$Year)
CP6 <- subset(CP5, CP5$LakeYear %in% AYSmp4$LakeYear)

CP6$BelowLN <- morelakes$BelowLN
CP6$Size <- morelakes$Area
CP6$Richness <- morelakes$Richness
CP6$Latfix <- morelakes$POINT_Y

CP6$Vattentemperatur<- na.approx(CP6$Vattentemperatur)
CP6$TOC<- na.approx(CP6$TOC)
CP6$"Tot-P"<- na.approx(CP6$"Tot-P")

CP6$`Alk/Acid`[CP6$`Alk/Acid`<=0] <- 0.0001

CP6$Lat2 <- scale(CP6$Stationskoordinat.N.X)
CP6$TOC2 <- scale(CP6$TOC)
CP6$TP <- scale(CP6$"Tot-P")
CP6$WT <- scale(CP6$Vattentemperatur)
CP6$Alk <- scale(CP6$`Alk/Acid`)
CP6$Size2 <- scale(CP6$Size)
CP6$Rich2 <- scale(CP6$Richness)

rownames(CP6) <- paste(CP6$Lake, CP6$Year)
Biochem <- cbind(CP6, SiteSc)
CP6$MeanRichness <- rep(tapply(CP6$Richness, factor(CP6$Lake), mean), times=table(factor(CP6$Lake)))
Biochem <- Biochem[Biochem$Lake!="Storvindeln",]
Biochem$Trecord <- rep(data.frame(table(MLsort$Lake))[c(1:86, 88:106),2], times=data.frame(table(MLsort$Lake))[c(1:86, 88:106),2])
InvBel <- Biochem[Biochem$BelowLN=="YES",]
BInoGud <- Biochem[Biochem$BelowLN=="NO",]

SpScDCA1 <- scores(DCAall, choices=1, display="species")
SpDC1 <- SpScDCA1[order(-abs(SpScDCA1))]
orderDC1 <- tax_name(names(SpDC1[1:15]),get="order", db = 'itis')
SpScDCA2 <- scores(DCAall, choices=2, display="species")
SpDC2 <- SpScDCA2[order(-abs(SpScDCA2))]
orderDC2 <- tax_name(names(SpDC2[1:15]),get="order", db = 'itis')

write.csv2(cbind(SpDC1[1:15], orderDC1, colSums(AYSmp3[names(SpDC1[1:15])]),names(SpDC2[1:15]),orderDC2, SpDC2[1:15],colSums(AYSmp3[names(SpDC2[1:15])])), "file:///C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Data Files/15SpDC1.csv" )
write.csv2(InvBel, "file:///C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Data Files/InvBel.csv" )


var <- tapply(Biochem$DCA1, factor(Biochem$Lake), var)
rich <- tapply(Biochem$Richness, factor(Biochem$Lake), mean)
lat<- tapply(Biochem$Latfix, factor(Biochem$Lake), mean)
yaaas <- cbind(var, lat, rich)


