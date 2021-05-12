#BETAPART STUFF###############################################################
library(betapart)
library(tidyverse)
library(vegan)
#Use sign function to convert abundance data 
#into prescence abscence matrix that can be fed to Betapart. Can't include factors (non numeric)
dat <- ordispe
dat.pa=sign(dat) #exclude them. works but lose first factor columns
dat.pa <- as.data.frame(dat.pa)
dat.pa <- rownames_to_column(dat.pa, var = "ID_fine2")
dat.pa$ID_fine2 <- as.factor(dat.pa$ID_fine2)
dat.pa=cbind(ID_site=0, dat.pa)
dat.pa=cbind(ID_plot=0, dat.pa)
dat.pa=cbind(ID_subplot=0, dat.pa) #add back some empty columns in that case...
dat.pa=cbind(survey_year=0, dat.pa)
dat.pa=cbind(ID_fine=0,dat.pa)
#dat.pa=cbind(ID=0, dat.pa)
dat.pa$survey_year <- test$survey_year #copy the values
dat.pa$ID_site <- test$ID_site
dat.pa$ID_plot <- test$ID_plot
dat.pa$ID_subplot <- test$ID_subplot
dat.pa$ID_fine <- test$ID_fine
#to put them back in.
str(dat.pa)

select(dat.pa, -c(2:6)) %>% group_by(ID_fine) %>% group_modify(beta.pair(.,index.family = "jaccard"))

#Stack overflow
repex <- select(dat.pa, 1, 6:16) %>% slice(.,1:10)
repex <- repex %>% rename(Site_ID = ID_fine, Plot_ID = ID_fine2)
repex$Site_ID <- c(1,1,1,1,1,2,2,2,2,2)
repex$Plot_ID <- c(1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5)
repex[1,3] <- 1
repex[6,6] <- 1
repex[8,7] <- 1

beta.pair(repex[,3:12])
repex %>% group_by(Site_ID) %>% summarise(beta = rowsum(.))
str(repex)
#Betapart doesn't like year and plot columns. Subset by year...
#fieldpresabs_1 created above
period1 <- c(1996,1997,1998,2000)
period2 <- c(2001,2002,2004,2005)
period3 <- c(2006,2007,2008,2010)
period4 <- c(2011,2012,2013,2015)

dat.pa.1 <- dat.pa %>% filter(Year %in% period1)  
dat.pa.2 <- dat.pa %>% filter(Year %in% period2)  
dat.pa.3 <- dat.pa %>% filter(Year %in% period3)   
dat.pa.4 <- dat.pa %>% filter(Year %in% period4)  

#and then subset by site
Gd1 <- filter(dat.pa.1, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Gd2 <- filter(dat.pa.2, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID) 
Gd3 <- filter(dat.pa.3, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Gd4 <- filter(dat.pa.4, Site == "Gd") %>% select(-Year, -Plot, -Site, -Tree, -ID)

An1 <- filter(dat.pa.1, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)
An2 <- filter(dat.pa.2, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)
An3 <- filter(dat.pa.3, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)
An4 <- filter(dat.pa.4, Site == "An") %>% select(-Year, -Plot, -Site, -Tree, -ID)

Ki1 <- filter(dat.pa.1, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ki2 <- filter(dat.pa.2, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ki3 <- filter(dat.pa.3, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ki4 <- filter(dat.pa.4, Site == "Ki") %>% select(-Year, -Plot, -Site, -Tree, -ID)

Ga1 <- filter(dat.pa.1, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ga2 <- filter(dat.pa.2, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ga3 <- filter(dat.pa.3, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)
Ga4 <- filter(dat.pa.4, Site == "Ga") %>% select(-Year, -Plot, -Site, -Tree, -ID)

#check things look right 
str(Gd1)
#remove zero sum columns #
Ga1 <- Ga1[, colSums(Ga1 != 0) > 0]
Ga2 <- Ga2[, colSums(Ga2 != 0) > 0]
Ga3 <- Ga3[, colSums(Ga3 != 0) > 0]
Ga4 <- Ga4[, colSums(Ga4 != 0) > 0]

Ki1 <- Ki1[, colSums(Ki1 != 0) > 0]
Ki2 <- Ki2[, colSums(Ki2 != 0) > 0]
Ki3 <- Ki3[, colSums(Ki3 != 0) > 0]
Ki4 <- Ki4[, colSums(Ki4 != 0) > 0]

An1 <- An1[, colSums(An1 != 0) > 0]
An2 <- An2[, colSums(An2 != 0) > 0]
An3 <- An3[, colSums(An3 != 0) > 0]
An4 <- An4[, colSums(An4 != 0) > 0]

Gd1 <- Gd1[, colSums(Gd1 != 0) > 0]
Gd2 <- Gd2[, colSums(Gd2 != 0) > 0]
Gd3 <- Gd3[, colSums(Gd3 != 0) > 0]
Gd4 <- Gd4[, colSums(Gd4 != 0) > 0]


# make betapart core objects to use in later analyses
Gd1c <- betapart.core(Gd1) 
Gd2c <- betapart.core(Gd2) 
Gd3c <- betapart.core(Gd3) 
Gd4c <- betapart.core(Gd4) 

An1c <- betapart.core(An1) 
An2c <- betapart.core(An2) 
An3c <- betapart.core(An3) 
An4c <- betapart.core(An4) 

Ki1c <- betapart.core(Ki1) 
Ki2c <- betapart.core(Ki2) 
Ki3c <- betapart.core(Ki3) 
Ki4c <- betapart.core(Ki4) 

Ga1c <- betapart.core(Ga1) 
Ga2c <- betapart.core(Ga2) 
Ga3c <- betapart.core(Ga3) 
Ga4c <- betapart.core(Ga4) 

#Returns three values-turnover and nestedness components, + total multi-plot dissimilarity across the site
(Ga1m <- beta.multi(Ga1c))#turnover=0.7966574, nestedness=0.06651038, betadiv=0.8631678
(Ga2m <- beta.multi(Ga2c))#turnover=0.8189655 nestedness=0.05272491, betadiv=0.8716904
(Ga3m <- beta.multi(Ga3c))#turnover=0.7884615 nestedness=0.07204571, betadiv=0.8605072
(Ga4m <- beta.multi(Ga4c))#turnover=0.808399 nestedness=0.06228307, betadiv=0.870682

(Gd1m <- beta.multi(Gd1c))#turnover=0.7558528, nestedness=0.0876627, betadiv=0.8435155
(Gd2m <- beta.multi(Gd2c))#turnover=0.7903226 nestedness=0.05601312, betadiv=0.8463357
(Gd3m <- beta.multi(Gd3c))#turnover=0.7916667 nestedness=0.05190875, betadiv=0.8435754
(Gd4m <- beta.multi(Gd4c))#turnover=0.6836158 nestedness=0.1351544, betadiv=0.8187702

(Ki1m <- beta.multi(Ki1c))#turnover=0.6054422, nestedness=0.2163704, betadiv=0.8218126
(Ki2m <- beta.multi(Ki2c))#turnover=0.5909091 nestedness=0.2266585, betadiv=0.8175676
(Ki3m <- beta.multi(Ki3c))#turnover=0.6462264 nestedness=0.1548346, betadiv=0.801061
(Ki4m <- beta.multi(Ki4c))#turnover=0.625 nestedness=0.1504011, betadiv=0.7754011

(An1m <- beta.multi(An1c))#turnover=0.5541401, nestedness=0.2324452, betadiv=0.7865854
(An2m <- beta.multi(An2c))#turnover=0.5510204 nestedness=0.2271309, betadiv=0.7781513
(An3m <- beta.multi(An3c))#turnover=0.6551724 nestedness=0.1876567, betadiv=0.8428291
(An4m <- beta.multi(An4c))#turnover=0.7286822 nestedness=0.1092808, betadiv=0.837963


#Make a beta.pair distance matrix for each beta.core object
Ga1p = beta.pair(Ga1c)
Ga2p = beta.pair(Ga2c)
Ga3p = beta.pair(Ga3c)
Ga4p = beta.pair(Ga4c)

Gd1p = beta.pair(Gd1c)
Gd2p = beta.pair(Gd2c)
Gd3p = beta.pair(Gd3c)
Gd4p = beta.pair(Gd4c)

An1p = beta.pair(An1c)
An2p = beta.pair(An2c)
An3p = beta.pair(An3c)
An4p = beta.pair(An4c)

Ki1p = beta.pair(Ki1c)
Ki2p = beta.pair(Ki2c)
Ki3p = beta.pair(Ki3c)
Ki4p = beta.pair(Ki4c)

#Compare first year with last
#Gammtratten
#turnover
Ga1p.sim <- Ga1p$beta.sim
Ga4p.sim <- Ga4p$beta.sim
t.test(Ga1p$beta.sim,Ga4p$beta.sim)
mantel(Ga1p.sim, Ga4p.sim, method = "pearson", permutations = 9999, na.rm = TRUE) # ** i.e. similar!
#nestedness
Ga1p.sne <- Ga1p$beta.sne
Ga4p.sne <- Ga4p$beta.sne
t.test(Ga1p$beta.sne,Ga4p$beta.sne)
mantel(Ga1p.sne, Ga4p.sne, method = "pearson", permutations = 9999, na.rm = TRUE) # **
#overall
Ga1p.sor <- Ga1p$beta.sor
Ga4p.sor <- Ga4p$beta.sor
t.test( Ga1p$beta.sor, Ga4p$beta.sor)
mantel(Ga1p.sor, Ga4p.sor, method = "pearson", permutations = 9999, na.rm = TRUE) # **

#Gårdsjön
#turnover
Gd1p.sim <- Gd1p$beta.sim
Gd4p.sim <- Gd4p$beta.sim
t.test(Gd1p$beta.sim,Gd4p$beta.sim)
mantel(Gd1p.sim, Gd4p.sim, method = "pearson", permutations = 9999, na.rm = TRUE) # NS i.e. dissimilar!
#nestedness
Gd1p.sne <- Gd1p$beta.sne
Gd4p.sne <- Gd4p$beta.sne
t.test(Gd1p$beta.sne,Gd4p$beta.sne)
mantel(Gd1p.sne, Gd4p.sne, method = "pearson", permutations = 9999, na.rm = TRUE) #  NS i.e. dissimilar!
#overall
Gd1p.sor <- Gd1p$beta.sor
Gd4p.sor <- Gd4p$beta.sor
t.test(Gd1p$beta.sor,Gd4p$beta.sor)
mantel(Gd1p.sor, Gd4p.sor, method = "pearson", permutations = 9999, na.rm = TRUE) # *

#Aneboda
#turnover
An1p.sim <- An1p$beta.sim
An4p.sim <- An4p$beta.sim
t.test(An1p$beta.sim,An4p$beta.sim)
mantel(An1p.sim, An4p.sim, method = "pearson", permutations = 9999, na.rm = TRUE) # NS i.e. dissimilar!
#nestedness
An1p.sne <- An1p$beta.sne
An4p.sne <- An4p$beta.sne
t.test(An1p$beta.sne,An4p$beta.sne)
mantel(An1p.sne, An4p.sne, method = "pearson", permutations = 9999, na.rm = TRUE) #  NS i.e. dissimilar!
#overall
An1p.sor <- An1p$beta.sor
An4p.sor <- An4p$beta.sor
t.test(An1p$beta.sor,An4p$beta.sor)
mantel(An1p.sor, An4p.sor, method = "pearson", permutations = 9999, na.rm = TRUE) # NS i.e dissimilar!

#Kindla
#turnover
Ki1p.sim <- Ki1p$beta.sim
Ki2p.sim <- Ki2p$beta.sim
Ki3p.sim <- Ki3p$beta.sim
Ki4p.sim <- Ki4p$beta.sim
t.test( Ki1p$beta.sim, Ki4p$beta.sim)
mantel(Ki1p.sim, Ki4p.sim, method = "pearson", permutations = 9999, na.rm = TRUE) # NS i.e. dissimilar!
#nestedness
Ki1p.sne <- Ki1p$beta.sne
Ki2p.sne <- Ki2p$beta.sne
Ki3p.sne <- Ki3p$beta.sne
Ki4p.sne <- Ki4p$beta.sne
t.test( Ki1p$beta.sne, Ki4p$beta.sne)
mantel(Ki1p.sne, Ki4p.sne, method = "pearson", permutations = 9999, na.rm = TRUE) #  *
#overall
Ki1p.sor <- Ki1p$beta.sor
Ki2p.sor <- Ki2p$beta.sor
Ki3p.sor <- Ki3p$beta.sor
Ki4p.sor <- Ki4p$beta.sor
t.test(Ki1p$beta.sor,Ki4p$beta.sor)
mantel(Ki1p.sor, Ki4p.sor, method = "pearson", permutations = 9999, na.rm = TRUE) # *


#add beta.pair results as vectors to a dataframe to plot results
#Gammtratten
`2000` <- as.vector(Ga1p$beta.sim)
`2005` <- as.vector(Ga2p$beta.sim)
`2010` <- as.vector(Ga3p$beta.sim)
`2015` <- as.vector(Ga4p$beta.sim)
Ga.sim <- tibble(`2000`,`2005`,`2010`,`2015`) %>%
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sim")
`2000` <- as.vector(Ga1p$beta.sne)
`2005` <- as.vector(Ga2p$beta.sne)
`2010` <- as.vector(Ga3p$beta.sne)
`2015` <- as.vector(Ga4p$beta.sne)
Ga.sne <- tibble(`2000`,`2005`,`2010`,`2015`) %>% 
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sne")
`2000` <- as.vector(Ga1p$beta.sor)
`2005` <- as.vector(Ga2p$beta.sor)
`2010` <- as.vector(Ga3p$beta.sor)
`2015` <- as.vector(Ga4p$beta.sor)
Ga.sor <- tibble(`2000`,`2005`,`2010`,`2015`) %>% 
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sor") 
Ga_plot <- bind_cols(Ga.sim, beta.sne=Ga.sne$beta.sne,beta.sor=Ga.sor$beta.sor)
Ga_plot$Year <- as_factor(Ga_plot$Year)

#Gårdsjön
`2000` <- as.vector(Gd1p$beta.sim)
`2005` <- as.vector(Gd2p$beta.sim)
`2010` <- as.vector(Gd3p$beta.sim)
`2015` <- as.vector(Gd4p$beta.sim)
Gd.sim <- tibble(`2000`,`2005`,`2010`,`2015`) %>%
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sim")
`2000` <- as.vector(Gd1p$beta.sne)
`2005` <- as.vector(Gd2p$beta.sne)
`2010` <- as.vector(Gd3p$beta.sne)
`2015` <- as.vector(Gd4p$beta.sne)
Gd.sne <- tibble(`2000`,`2005`,`2010`,`2015`) %>% 
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sne")
`2000` <- as.vector(Gd1p$beta.sor)
`2005` <- as.vector(Gd2p$beta.sor)
`2010` <- as.vector(Gd3p$beta.sor)
`2015` <- as.vector(Gd4p$beta.sor)
Gd.sor <- tibble(`2000`,`2005`,`2010`,`2015`) %>% 
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sor") 
Gd_plot <- bind_cols(Gd.sim, beta.sne=Gd.sne$beta.sne,beta.sor=Gd.sor$beta.sor)
Gd_plot$Year <- as_factor(Gd_plot$Year)


#Kindla
`2000` <- as.vector(Ki1p$beta.sim)
`2005` <- as.vector(Ki2p$beta.sim)
`2010` <- as.vector(Ki3p$beta.sim)
`2015` <- as.vector(Ki4p$beta.sim)
Ki.sim <- tibble(`2000`,`2005`,`2010`,`2015`) %>%
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sim")
`2000` <- as.vector(Ki1p$beta.sne)
`2005` <- as.vector(Ki2p$beta.sne)
`2010` <- as.vector(Ki3p$beta.sne)
`2015` <- as.vector(Ki4p$beta.sne)
Ki.sne <- tibble(`2000`,`2005`,`2010`,`2015`) %>% 
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sne")
`2000` <- as.vector(Ki1p$beta.sor)
`2005` <- as.vector(Ki2p$beta.sor)
`2010` <- as.vector(Ki3p$beta.sor)
`2015` <- as.vector(Ki4p$beta.sor)
Ki.sor <- tibble(`2000`,`2005`,`2010`,`2015`) %>% 
  pivot_longer(everything(),names_to = "Year", values_to = "beta.sor") 
Ki_plot <- bind_cols(Ki.sim, beta.sne=Ki.sne$beta.sne,beta.sor=Ki.sor$beta.sor)
Ki_plot$Year <- as_factor(Ki_plot$Year)

# Anebodaxc 
#violin plots
gg <- ggplot(Gd_plot, aes(x = Year,
                          y = beta.sne))
gg +        
  geom_violin(fill='grey', color="darkgrey")+
  stat_summary(fun=median, geom="point", shape = 21, colour = "black", fill = "white", size = 3, stroke = 1)+
  ggtitle("") +
  scale_x_discrete() + 
  scale_y_continuous(limits = c(0, 1))+
  theme_bw() +
  theme(legend.position="none")+
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        plot.title = element_text(size=14, face="bold"))




#boxplot turnover (pass 1 to function), nestedness (2), or total (3)
#Gammtratten
fbeta.pair= function(x){
  y=c(Ga1p[x],Ga2p[x], Ga3p[x], Ga4p[x])
  boxplot(y)
  title(main="GA,Beta.pair,4 yrs")
}

par(mfrow=c(1,3))
fbeta.pair(1) 

#Gårdsjön
fbeta.pair= function(x){
  y=c(Gd1p[x],Gd2p[x], Gd3p[x], Gd4p[x])
  boxplot(y)
  title(main="Gd,Beta.pair,4 yrs")
}

par(mfrow=c(1,3))
fbeta.pair(1) 

#Kindla
fbeta.pair= function(x){
  y=c(Ki1p[x],Ki2p[x], Ki3p[x], Ki4p[x])
  boxplot(y)
  title(main="Ki")
}

par(mfrow=c(1,3))
fbeta.pair(1) 

#Aneboda
fbeta.pair= function(x){
  y=c(An1p[x],An2p[x], An3p[x], An4p[x])
  boxplot(y)
  title(main="An,Beta.pair,4 yrs")
}

par(mfrow=c(1,3))
fbeta.pair(1) 

#boxplot turnover (pass 1 to function), nestedness (2), or total (3)
fbeta.multi= function(x){
  y=c(Ga1m[x], Ga2m[x], Ga3m[x], Ga4m[x])
  boxplot(y)
  title(main="Beta.multi, all 4 sample years")
}

fbeta.multi(1)

t.test(An1p[1], An4p[1])
v1 <- vegdist(Gd1p$beta.sor)
v2 <- vegdist(Gd4p$beta.sor)
mantel(Gd1p$beta.sor, Gd4p$beta.sor, method = "pearson", permutations = 9999, na.rm = TRUE) 
mantel(xdis = v1, ydis = v2, method = "pearson", permutations = 9999, na.rm = TRUE) 

#temporal change between years
#apply betapart temp (must be same sized matrices)
#ID means different site names for each year. Remove rownames and run using row numbers only.
Gd1t <- Gd1
rownames(Gd1t) <- c(1:20)
Gd4t <- Gd4
rownames(Gd4t) <- c(1:20)

Ga1t <- Ga1
rownames(Ga1t) <- c(1:20)
Ga4t <- Ga4
rownames(Ga4t) <- c(1:20)

Ki1t <- Ki1
rownames(Ki1t) <- c(1:20)
Ki4t <- Ki4
rownames(Ki4t) <- c(1:20)

An1t <- An1
rownames(An1t) <- c(1:20)
An4t <- An4
rownames(An4t) <- c(1:20)

betapart_t_Gd <-  beta.temp(Gd1t, Gd4t)
betapart_t_An <-  beta.temp(An1t, An4t)
betapart_t_Ki <-  beta.temp(Ki1t, Ki4t)
betapart_t_Ga <-  beta.temp(Ga1t, Ga4t)

par(mfrow=c(1,4))
boxplot(betapart_t_Gd)
title("betapart_t_Gd")
boxplot(betapart_t_An)
title("betapart_t_An")
boxplot(betapart_t_Ki)
title("betapart_t_Ki")
boxplot(betapart_t_Ga)
title("betapart_t_Ga")

#plot probability distributions####
# sampling across equal sites (change to the beta.core object for the site of interest)
samp_1 <- beta.sample(Ki1c, sites=15, samples=100, index.family="sor")
samp_2 <- beta.sample(Ki2c, sites=15, samples=100, index.family="sor") 
samp_3 <- beta.sample(Ki3c, sites=15, samples=100, index.family="sor") 
samp_4 <- beta.sample(Ki4c, sites=15, samples=100, index.family="sor") 

# plotting the distributions of components 
dist_1 <- samp_1$sampled.values
dist_2 <- samp_2$sampled.values 
dist_3 <- samp_3$sampled.values 
dist_4 <- samp_4$sampled.values 

#Aneboda
par(mfrow=c(1,1))
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0,40), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Aneboda" )
# Add a legend
legend(0.6, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.6, 12, legend=c("1997", "2002", "2007", "2012"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

# Compute the p-value of difference in beta.SIM
p.value.beta.SIM<-length(which(samp_1$sampled.values$beta.SIM>
                                 samp_4$sampled.values$beta.SIM))/100

p.value.beta.SNE<-length(which(samp_1$sampled.values$beta.SNE<
                                 samp_4$sampled.values$beta.SNE))/100

p.value.beta.SOR<-length(which(samp_1$sampled.values$beta.SOR>
                                 samp_4$sampled.values$beta.SOR))/100
p.value.beta.SIM # 10 0.04 * # 15 0.001 **
p.value.beta.SNE # 10 0.05 # 15 0.001 **
p.value.beta.SOR # 10 0.13 # 15 0.013 *

#Gårdsjön
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0, 17), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Gårdsjön" )
# Add a legend
legend(0.7, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.7, 12, legend=c("1996", "2001", "2006", "2011"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

# Compute the p-value of difference in beta.SIM
p.value.beta.SIM<-length(which(samp_1$sampled.values$beta.SIM<
                                 samp_4$sampled.values$beta.SIM))/100

p.value.beta.SNE<-length(which(samp_1$sampled.values$beta.SNE>
                                 samp_4$sampled.values$beta.SNE))/100

p.value.beta.SOR<-length(which(samp_1$sampled.values$beta.SOR<
                                 samp_4$sampled.values$beta.SOR))/100
p.value.beta.SIM #10 0.13 # 15 0.031*
p.value.beta.SNE #10 0.16 # 15 0.027*
p.value.beta.SOR #10 0.19 # 15 0.088 

#Gammtratten
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0, 18), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Gammtratten" )
# Add a legend
legend(0.5, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.5, 12, legend=c("2000", "2005", "2010", "2015"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

# Compute the p-value of difference in beta.SIM
p.value.beta.SIM<-length(which(samp_1$sampled.values$beta.SIM>
                                 samp_4$sampled.values$beta.SIM))/100

p.value.beta.SNE<-length(which(samp_1$sampled.values$beta.SNE<
                                 samp_4$sampled.values$beta.SNE))/100

p.value.beta.SOR<-length(which(samp_1$sampled.values$beta.SOR>
                                 samp_4$sampled.values$beta.SOR))/100

p.value.beta.SIM # 10 0.33 # 15 0.336
p.value.beta.SNE # 10 0.38 # 15 0.4
p.value.beta.SOR # 10 0.35 # 15 0.304

#Kindla
plot(density(
  dist_1$beta.SNE), xlim=c(0,0.9), ylim=c(0, 35), xlab='Beta diversity', main='', col='grey60',lwd=2)
lines(density(dist_1$beta.SOR),col='grey60', lwd=2) 
lines(density(dist_1$beta.SNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.SIM),col='grey60', lty=2, lwd=2)
lines(density(dist_2$beta.SOR),col='red', lwd=2) 
lines(density(dist_2$beta.SNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.SIM),col='red', lty=2,lwd=2) 
lines(density(dist_3$beta.SOR),col='green', lwd=2) 
lines(density(dist_3$beta.SNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.SIM),col='green', lty=2, lwd=2)
lines(density(dist_4$beta.SOR),col='blue', lwd=2) 
lines(density(dist_4$beta.SNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.SIM),col='blue', lty=2, lwd=2)
title("Kindla" )
# Add a legend
legend(0.5, 15, legend=c("Turnover", "Nestedness"), lty= c(1,2), cex=0.8)
legend(0.5, 12, legend=c("1998", "2004", "2008", "2013"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)

# Compute the p-value of difference in beta.SIM
p.value.beta.SIM<-length(which(samp_1$sampled.values$beta.SIM>
                                 samp_4$sampled.values$beta.SIM))/100

p.value.beta.SNE<-length(which(samp_1$sampled.values$beta.SNE<
                                 samp_4$sampled.values$beta.SNE))/100

p.value.beta.SOR<-length(which(samp_1$sampled.values$beta.SOR<
                                 samp_4$sampled.values$beta.SOR))/100

p.value.beta.SIM # 10 0.47 # 15 0.405
p.value.beta.SNE # 10 0.13 # 15 0.072
p.value.beta.SOR # 10 0.11 # 15 0.033

#alternative Jaccard
par(mfrow=c(1,1))
plot(density(
  dist_1$beta.JNE), xlim=c(0,0.9), ylim=c(0, 22), xlab='Beta diversity', main='', col='grey60',lwd=2)
#lines(density(dist_1$beta.JAC),col='grey60', lwd=2) 
lines(density(dist_1$beta.JNE),col='grey60', lty=1, lwd=2) 
lines(density(dist_1$beta.JTU),col='grey60', lty=2, lwd=2)
#lines(density(dist_2$beta.JAC),col='red', lwd=2) 
lines(density(dist_2$beta.JNE),col='red', lty=1,lwd=2) 
lines(density(dist_2$beta.JTU),col='red', lty=2,lwd=2) 
#lines(density(dist_3$beta.JAC),col='green', lwd=2) 
lines(density(dist_3$beta.JNE),col='green', lty=1, lwd=2) 
lines(density(dist_3$beta.JTU),col='green', lty=2, lwd=2)
#lines(density(dist_4$beta.JAC),col='blue', lwd=2) 
lines(density(dist_4$beta.JNE),col='blue', lty=1, lwd=2) 
lines(density(dist_4$beta.JTU),col='blue', lty=2, lwd=2)
#title("gray=1,r=2, g=3, b=4, right=total, left=nest, dash=turn" )
title("Gårdsjön")
# Add a legend
legend(0.7, 18, legend=c("1997", "2002", "2007", "2012"),
       col=c("grey60","red","green", "blue"), lty=1, cex=0.8)
text(0.15, 18, labels = "Turnover")
text(0.65, 10, labels = "Nestedness")

# Now, we'll calculate the Jaccard index and its partitions of turnover and nestedness. 
#We can calculate Sorensen index instead by using the argument     index.family="sorensen"    .
dist<-beta.pair(select(dat.pa.4,-Year, -Plot, -Site, -Tree, -ID), index.family="sorensen")
dist<-bray.part(dat)
# To get the pairwise Jaccard index turnover partition between communities, type: dist[[1]]. 
#To get nestedness partition, type: dist[[2]]. To get all beta diversity: dist[[3]].
groups <- data$Site
bd<-betadisper(dist[[1]],groups)
plot(bd)
boxplot(bd)
anova(bd)


plot(hclust(Ga1p$beta.sor, method="average"), hang=-1)
#exploratory PCA####
group <- data4$Site
dat <- select(data4, -c(ID, Site, Year, Plot, Tree, TInd, RInd, pH, id2, ID3, sp))

fit <- rda(dat, scale = TRUE)
pl <- ordiplot(fit, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "species", col="blue", cex=0.9)
ordihull(fit, groups = data4$Site, label = TRUE)






