require(rjags)
require(MCMCglmm)
library(MCMCvis)
library(tidyverse)

n.plots<-length(unique(dat$ID_site))#184
n.obs<-as.vector(table(factor(dat$ID_site)))
precipmean<-as.vector(tapply(dat$mean_precip, factor(dat$ID_site), mean))
latitudemean <- as.vector(tapply(dat$latitude, factor(dat$ID_site), mean))
canopymean <- as.vector(tapply(dat$sum_canopy, factor(dat$ID_site), mean))
NH4mean <- as.vector(tapply(dat$n_nh4, factor(dat$ID_site), mean))
NO3mean <- as.vector(tapply(dat$n_no3, factor(dat$ID_site), mean))
#Trecord <- as.vector(tapply(scale(dat$Trecord), factor(dat$ID_site), mean))

precip <-dat$mean_precip
latitude <- dat$latitude
NO3 <- dat$n_no3
NH4 <- dat$n_nh4
Canopy <- dat$sum_canopy


model <- "model{
for(i in 1:n) {
y[i]~dnorm(y.hat[i], tau.plot[plot[i]])
y.hat[i]<- mu + l[plot[i]] + ye[year[i]] + precip[i]*p + latitude[i]*lat +
Canopy[i]*canopyness + NO3[i]*NO3p + NH4[i]*NH4p 
}

mu ~ dnorm (0, .0001)

for (j in 1:L){
l[j] ~ dnorm (0, tau.l)
}

for (k in 1:Y){
ye[k] ~ dnorm (0, tau.y)
}


for(j in 1:L){
tau.plot[j]<-pow(sigma.plot[j], -2)
sigma.plot[j]~dlnorm(mu.plot.sigma[j], res.tau.plot)
mu.plot.sigma[j] <- mu2 + precipmean[j]*pmean + latitudemean[j]*latmean +
canopymean[j]*rmean + NH4mean[j]*nh4mean + NO3mean[j]*no3mean 
}

res.tau.plot <- pow(res.sigma.plot, -2)
res.sigma.plot ~ dunif(0, 10)
tau.l <- pow(sigma.l, -2)
sigma.l ~ dunif(0, 10)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 10)
mu2 ~ dnorm (0, .0001)
p ~ dnorm (0, .0001)
pmean ~ dnorm (0, .0001)
lat ~ dnorm(0, .0001)
latmean ~ dnorm(0, .0001)
NO3p ~ dnorm(0, .0001)
NH4p ~ dnorm(0, .0001)
canopyness ~ dnorm(0, .0001)
rmean ~ dnorm(0, .0001)
nh4mean ~ dnorm(0, .0001)
no3mean ~ dnorm(0, .0001)
}"

##PC1
writeLines(model, "JAGSmodelbelowPC1.txt")
belowPC1 <- list(y =  dat$PC1,  n = length(dat$PC1), L=n.plots, Y=length(unique(dat$survey_year)),
                  plot=rep(1:n.plots, n.obs), year=dat$survey_year-1993,
                  precip=as.vector(precip), latitude=as.vector(latitude), NO3=as.vector(NO3),
                  NH4=as.vector(NH4), Canopy=as.vector(Canopy), precipmean=precipmean, latitudemean=latitudemean,
                  canopymean=canopymean, NH4mean=NH4mean, NO3mean=NO3mean)
jags1 <- jags.model("JAGSmodelbelowPC1.txt" , data=belowPC1, n.chains=2, n.adapt=1000)
smp1 <- coda.samples(jags1, c("mu2","res.sigma.plot","p", "lat", "NO3p", "NH4p",
                              "canopyness", "pmean", "latmean", "rmean", "nh4mean", "no3mean",
                              "sigma.l", "sigma.y"), n.iter = 30500, thin=30, burnin=50,
                     n.chains=2, inits=inits) #, "sigma.plot"

# smp2<-as.mcmc(smp1[1])
# CI <- apply(smp2, 2, quantile, c(0.075, 0.925))
# meansmp <- apply(smp2, 2, mean)
MCMCsummary(smp1)
MCMCplot(smp1, ref_ovl = TRUE, rank=FALSE)
#MCMCtrace(smp1)


##DCA2
writeLines(model, "JAGSmodelbelowDCA2.txt")
belowPC2 <- list(y =  dat$PC2,  n = length(dat$PC2), L=n.plots, Y=length(unique(dat$survey_year)),
                  plot=rep(1:n.plots, n.obs), year=dat$survey_year-1993,
                  precip=as.vector(precip), latitude=as.vector(latitude), NO3=as.vector(NO3),
                  NH4=as.vector(NH4), Canopy=as.vector(Canopy), precipmean=precipmean, latitudemean=latitudemean,
                  canopymean=canopymean, NH4mean=NH4mean, NO3mean=NO3mean)
jags2 <- jags.model("JAGSmodelbelowDCA2.txt" , data=belowPC2, n.chains=2, n.adapt=1000)
smp3 <- coda.samples(jags2, c("mu2","res.sigma.plot","p", "lat", "NO3p", "NH4p",
                              "canopyness", "pmean", "latmean", "rmean", "nh4mean", "no3mean",
                              "sigma.l", "sigma.y"), n.iter = 30500, thin=30, burnin=50,
                     n.chains=2, inits=inits)#, "sigma.plot"

# smp4<-as.mcmc(smp3[1])
# CI2 <- apply(smp4, 2, quantile, c(0.075, 0.925))
# meansmp2 <- apply(smp4, 2, mean)

MCMCsummary(smp3)
MCMCplot(smp3, ref_ovl = TRUE, rank=FALSE)






write.csv2(cbind(paste(names(meansmp[c(4,10,8,1,2,11,69,71)])),
                 paste(round(meansmp[c(4,10,8,1,2,11,69,71)],2), 
                       paste=" (", round(CI[1,c(4,10,8,1,2,11,69,71)],2), 
                       paste=", ", round(CI[2,c(4,10,8,1,2,11,69,71)],2), 
                       paste=")", sep = ""), paste(round(meansmp2[c(4,10,8,1,2,11,69,71)],2), 
                                                   paste=" (",round(CI2[1, c(4,10,8,1,2,11,69,71)],2)
                                                   , paste=", ", round(CI2[2, c(4,10,8,1,2,11,69,71)],2), 
                                                   paste=")", sep = "")), 
           "COMPOSITIONBELOW.csv")



write.csv2(cbind(paste(names(meansmp[c(5,70,9,3,72,73)])),
                 paste(round(meansmp[c(5,70,9,3,72,73)],2), 
                       paste=" (", round(CI[1,c(5,70,9,3,72,73)],2), 
                       paste=", ", round(CI[2,c(5,70,9,3,72,73)],2), 
                       paste=")", sep = ""), paste(round(meansmp2[c(5,70,9,3,72,73)],2), 
                                                   paste=" (",round(CI2[1, c(5,70,9,3,72,73)],2)
                                                   , paste=", ", round(CI2[2, c(5,70,9,3,72,73)],2), 
                                                   paste=")", sep = "")), 
           "STABILITYBELOW.csv")


write.csv2(cbind(paste(unique(factor(dat$ID_site))), 
                 paste(round(meansmp[12:68]^2,2), paste=" (",round(CI[1,12:68]^2,2), 
                       paste=", ", round(CI[2,12:68]^2,2), paste=")", sep = ""), 
                 paste(round(meansmp2[12:68]^2,2), paste=" (",round(CI2[1,12:68]^2,2), 
                       paste=", ", round(CI2[2,12:62]^2,2), paste=")", sep = "")), 
           "HRsBELOW.csv")

##Canopyness PC1
canopymeanBL <- as.vector(tapply(scale(dat$sum_canopy), factor(dat$ID_site), mean))
rBLN <- seq(min(canopymeanBL), max(canopymeanBL), by=0.01)
yR <- matrix(NA, 1016, length(rBLN))
for (i in 1:1016) {
  yR[i,]=exp(smp2[i,6] + smp2[i,9]*rBLN)^2
}
my <- apply(log(yR),2,mean)
ly <- apply(log(yR),2,quantile,0.025)
hy <- apply(log(yR),2,quantile,0.975)

##Canopyness DCA2
canopymeanBL <- as.vector(tapply(scale(dat$sum_canopy), factor(dat$ID_site), mean))
rBLN <- seq(min(canopymeanBL), max(canopymeanBL), by=0.01)
yRR <- matrix(NA, 1016, length(rBLN))
for (i in 1:1016) {
  yRR[i,]=exp(smp4[i,6] + smp4[i,9]*rBLN)^2
}
myR <- apply(log(yRR),2,mean)
lyR <- apply(log(yRR),2,quantile,0.025)
hyR <- apply(log(yRR),2,quantile,0.975)


##NH4inity PC1
NH4meanBL <-  as.vector(tapply(log(dat$n_nh4), factor(dat$ID_site), mean))
NH4BLN <- seq(min(NH4meanBL), max(NH4meanBL), by=0.01)
NH4AB2 <- matrix(NA, 1016, length(NH4BLN))
for (i in 1:1016) {
  NH4AB2[i,]=exp(smp2[i,6] + smp2[i,3]*NH4BLN)^2
}
mx3 <- apply(log(NH4AB2),2,mean)
lx3 <- apply(log(NH4AB2),2,quantile,0.025)
hx3 <- apply(log(NH4AB2),2,quantile,0.975)

##NO3 PC1
NO3meanBL <-  as.vector(tapply(log(dat$n_no3), factor(dat$ID_site), mean))
NO3BLN <- seq(min(NO3meanBL), max(NO3meanBL), by=0.01)
NO3AB2 <- matrix(NA, 1016, length(NO3BLN))
for (i in 1:1016) {
  NO3AB2[i,]=exp(smp2[i,6] + smp2[i,71]*NO3BLN)^2
}
mx4 <- apply(log(NO3AB2),2,mean)
lx4 <- apply(log(NO3AB2),2,quantile,0.025)
hx4 <- apply(log(NO3AB2),2,quantile,0.975)

##precip DCA2
SmeanBL <-  as.vector(tapply(log(dat$mean_precip), factor(dat$ID_site), mean))
SBLN <- seq(min(SmeanBL), max(SmeanBL), by=0.01)
precip2 <- matrix(NA, 1016, length(SBLN))
for (i in 1:1016) {
  precip2[i,]=exp(smp4[i,6] + smp4[i,70]*SBLN)^2
}
ms2 <- apply(log(precip2),2,mean)
ls2 <- apply(log(precip2),2,quantile,0.025)
hs2 <- apply(log(precip2),2,quantile,0.975)

##Temp PC1
TempmeanBL <-  as.vector(tapply(log(dat$Vattentemperatur), factor(dat$ID_site), mean))
TempBLN <- seq(min(TempmeanBL), max(TempmeanBL), by=0.01)
Temp2 <- matrix(NA, 1016, length(TempBLN))
for (i in 1:1016) {
  Temp2[i,]=exp(smp2[i,6] + smp2[i,72]*TempBLN)^2
}
mt3 <- apply(log(Temp2),2,mean)
lt3 <- apply(log(Temp2),2,quantile,0.025)
ht3 <- apply(log(Temp2),2,quantile,0.975)

##Temp DCA2
Temp3 <- matrix(NA, 1016, length(TempBLN))
for (i in 1:1016) {
  Temp3[i,]=exp(smp4[i,6] + smp4[i,72]*TempBLN)^2
}
mt4 <- apply(log(Temp3),2,mean)
lt4 <- apply(log(Temp3),2,quantile,0.025)
ht4 <- apply(log(Temp3),2,quantile,0.975)


varPC1BELOW <- tapply(dat$PC1, factor(dat$ID_site), var)
varDCA2BELOW <- tapply(dat$DCA2, factor(dat$ID_site), var)

CanopyBELOW <- scale(dat$sum_canopy)
NH4BELOW <- log(dat$n_nh4)
precipBELOW <- scale(dat$mean_precip)
