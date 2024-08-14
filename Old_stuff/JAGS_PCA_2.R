library(rjags)
library(MCMCglmm)
library(MCMCvis)
library(tidyverse)
options(mc.cores = parallel::detectCores())
#must first read in dat data
dat <- readRDS("dat.RDS")
#dat <- drop_na(dat)

n.plots<-length(unique(dat$ID_site))#184
n.obs<-as.vector(table(dat$ID_site))
canopymean<-as.vector(tapply(dat$sum_canopy, factor(dat$ID_site), mean))
latitudemean <- as.vector(tapply(scale(dat$latitude), factor(dat$ID_site), mean))
precipmean <- as.vector(tapply(scale(dat$mean_precip), factor(dat$ID_site), mean))
NH4mean <- as.vector(tapply(dat$n_nh4, factor(dat$ID_site), mean))
NO3mean <- as.vector(tapply(dat$n_no3, factor(dat$ID_site), mean))
tempmean <- as.vector(tapply(scale(dat$mean_temp), factor(dat$ID_site), mean))


canopy <-dat$sum_canopy
latitude <- scale(dat$latitude)
NO3 <- dat$n_no3
NH4 <- dat$n_nh4
precip <- scale(dat$mean_precip)
temp <- scale(dat$mean_temp)

model <- "model{
for(i in 1:n) {
y[i]~dnorm(y.hat[i], tau.plot[plot[i]])
y.hat[i]<- mu + l[plot[i]] + ye[year[i]] + canopy[i]*s + latitude[i]*lat +
precip[i]*precipitation + NO3[i]*nitrate + NH4[i]*ammonium + temp[i]*temperature
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
mu.plot.sigma[j] <- mu2 + canopymean[j]*cmean + latitudemean[j]*latmean +
precipmean[j]*pmean + NH4mean[j]*nh4_mean + NO3mean[j]*no3_mean + tempmean[j]*tmean
}

res.tau.plot <- pow(res.sigma.plot, -2)
res.sigma.plot ~ dunif(0, 10)
tau.l <- pow(sigma.l, -2)
sigma.l ~ dunif(0, 10)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 10)
mu2 ~ dnorm (0, .0001)
s ~ dnorm (0, .0001)
cmean ~ dnorm (0, .0001)
lat ~ dnorm(0, .0001)
latmean ~ dnorm(0, .0001)
nitrate ~ dnorm(0, .0001)
ammonium ~ dnorm(0, .0001)
precipitation ~ dnorm(0, .0001)
pmean ~ dnorm(0, .0001)
nh4_mean ~ dnorm(0, .0001)
no3_mean ~ dnorm(0, .0001)
tmean ~ dnorm(0, .0001)
temperature ~ dnorm(0, .0001)
}"

##PCA1
writeLines(model, "JAGSmodelbelowPCA1.txt")
belowPCA1 <- list(y =  dat$PC1,  n = length(dat$PC1), L=n.plots, Y=length(unique(dat$survey_year)),
                  plot=rep(1:n.plots, n.obs),
                  year=as.numeric(levels(as.factor(dat$survey_year)))[as.factor(dat$survey_year)]-1993,
                  canopy=as.vector(canopy), latitude=as.vector(latitude), NO3=as.vector(NO3), NH4=as.vector(NH4),
                  precip=as.vector(precip), temp=as.vector(temp), canopymean=canopymean, latitudemean=latitudemean,
                  precipmean=precipmean, NH4mean=NH4mean, NO3mean=NO3mean, tempmean=tempmean)
jags1 <- jags.model("JAGSmodelbelowPCA1.txt" , data=belowPCA1, n.chains=2, n.adapt=1000)
smp1 <- coda.samples(jags1, c("mu2", "sigma.plot","res.sigma.plot","s", "lat",
                              "nitrate", "ammonium","precipitation", "cmean", "latmean", "pmean",
                              "NH4mean", "NO3mean", "tempmean","temperature","sigma.l", "sigma.y"),
                     n.iter = 30500, thin=30, burnin=500, n.chains=2, inits=inits)

smp2<-as.mcmc(smp1[1])
CI <- apply(smp2, 2, quantile, c(0.075, 0.925))
meansmp <- apply(smp2, 2, mean)


MCMCsummary(smp1)
MCMCplot(smp1, ref_ovl = TRUE, rank=FALSE)

##PCA2
writeLines(model, "JAGSmodelbelowPCA2.txt")
belowPCA2 <- list(y = dat$PC2,  n = length(dat$PC2), L=n.plots, Y=length(unique(dat$survey_year)),
                  plot=rep(1:n.plots, n.obs),
                  year=as.numeric(levels(as.factor(dat$survey_year)))[as.factor(dat$survey_year)]-1993,
                  canopy=as.vector(canopy),
                  latitude=as.vector(latitude), NO3=as.vector(NO3), NH4=as.vector(NH4),
                  precip=as.vector(precip), temp=as.vector(temp), canopymean=canopymean,
                  latitudemean=latitudemean, precipmean=precipmean, NH4mean=NH4mean,
                  NO3mean=NO3mean, tempmean=tempmean)
jags2 <- jags.model("JAGSmodelbelowPCA2.txt" , data=belowPCA2, n.chains=2, n.adapt=1000)
smp3 <- coda.samples(jags2, c("mu2", "sigma.plot","res.sigma.plot","s", "lat", "TotP",
                              "Alkal","richness", "smean", "latmean", "rmean", "NH4mean",
                              "NO3mean", "tempmean","temperature","sigma.l", "sigma.y"), 
                     n.iter = 30500, thin=30, burnin=500, n.chains=2, inits=inits)

smp4<-as.mcmc(smp3[1])
CI2 <- apply(smp4, 2, quantile, c(0.075, 0.925))
meansmp2 <- apply(smp4, 2, mean)

MCMCsummary(smp3)
MCMCplot(smp3, ref_ovl = TRUE, rank=FALSE)
