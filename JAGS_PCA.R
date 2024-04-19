require(rjags)
require(MCMCglmm)
library(MCMCvis)
dat <- all_data2
#must first read in dat data

n.lakes<-241
n.obs<-as.vector(table(factor(dat$ID_site)))
#sizemean<-as.vector(tapply(log(dat$Size), factor(dat$ID_site), mean))
latitudemean <- as.vector(tapply(dat$latitude, factor(dat$ID_site), mean))
tempmean <- as.vector(tapply(scale(dat$mean_temp), factor(dat$ID_site), mean))
precipmean <- as.vector(tapply(scale(dat$mean_precip), factor(dat$ID_site), mean))
NO3mean <- as.vector(tapply(dat$n_no3, factor(dat$ID_site), mean))
NH4mean <- as.vector(tapply(dat$n_nh4, factor(dat$ID_site), mean))
canopymean <- as.vector(tapply(dat$sum_canopy, factor(dat$ID_site), mean))

lat <- dat$latitude
temp <- dat$mean_temp
precip <- dat$mean_precip
NO3 <- dat$n_no3
NH4 <- dat$n_nh4
canopy <- dat$sum_canopy

Alk <- log(dat$`Alk/Acid`)
Rich <- scale(dat$Richness)
Temp <- scale(dat$WaterTemp)

model <- "model{
for(i in 1:n) {
y[i]~dnorm(y.hat[i], tau.ID_site[ID_site[i]])
y.hat[i]<- mu + l[ID_site[i]] + ye[year[i]] + size[i]*s + latitude[i]*lat +
Rich[i]*richness + TP[i]*TotP + Alk[i]*Alkal + Temp[i]*temperature
}

mu ~ dnorm (0, .0001)

for (j in 1:L){
l[j] ~ dnorm (0, tau.l)
}

for (k in 1:Y){
ye[k] ~ dnorm (0, tau.y)
}


for(j in 1:L){
tau.ID_site[j]<-pow(sigma.ID_site[j], -2)
sigma.ID_site[j]~dlnorm(mu.ID_site.sigma[j], res.tau.ID_site)
mu.ID_site.sigma[j] <- mu2 + sizemean[j]*smean + latitudemean[j]*latmean +
richmean[j]*rmean + ALKmean[j]*alkmean + TPmean[j]*tpmean + Tmean[j]*tempmean
}

res.tau.ID_site <- pow(res.sigma.ID_site, -2)
res.sigma.ID_site ~ dunif(0, 10)
tau.l <- pow(sigma.l, -2)
sigma.l ~ dunif(0, 10)
tau.y <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 10)
mu2 ~ dnorm (0, .0001)
s ~ dnorm (0, .0001)
smean ~ dnorm (0, .0001)
lat ~ dnorm(0, .0001)
latmean ~ dnorm(0, .0001)
TotP ~ dnorm(0, .0001)
Alkal ~ dnorm(0, .0001)
richness ~ dnorm(0, .0001)
rmean ~ dnorm(0, .0001)
alkmean ~ dnorm(0, .0001)
tpmean ~ dnorm(0, .0001)
tempmean ~ dnorm(0, .0001)
temperature ~ dnorm(0, .0001)
}"

##DCA1
writeLines(model, "JAGSmodelbelowDCA1.txt")
belowDCA1 <- list(y =  dat$DCA1,  n = length(dat$DCA1), L=n.ID_sites, Y=23,
                  ID_site=rep(1:n.ID_sites, n.obs), year=as.numeric(levels(as.factor(dat$Year)))[as.factor(dat$Year)]-1994,
                  size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk),
                  Rich=as.vector(Rich), Temp=as.vector(Temp), sizemean=sizemean, latitudemean=latitudemean,
                  richmean=richmean, ALKmean=ALKmean, TPmean=TPmean, Tmean=Tmean)
jags1 <- jags.model("JAGSmodelbelowDCA1.txt" , data=belowDCA1, n.chains=2, n.adapt=1000)
smp1 <- coda.samples(jags1, c("mu2", "sigma.ID_site","res.sigma.ID_site","s", "lat",
                              "TotP", "Alkal","richness", "smean", "latmean", "rmean",
                              "alkmean", "tpmean", "tempmean","temperature","sigma.l", "sigma.y"),
                     n.iter = 3050000, thin=3000, burnin=50000, n.chains=1, inits=inits)

smp2<-as.mcmc(smp1[1])
CI <- apply(smp2, 2, quantile, c(0.075, 0.925))
meansmp <- apply(smp2, 2, mean)

MCMCsummary(smp1, round=2)
MCMCtrace(smp1)
MCMCplot(smp1, ref_ovl = TRUE)
##DCA2
writeLines(model, "JAGSmodelbelowDCA2.txt")
belowDCA2 <- list(y = dat$DCA2,  n = length(dat$DCA2), L=n.ID_sites, Y=23, ID_site=rep(1:n.ID_sites, n.obs), year=as.numeric(levels(dat$Year))[dat$Year]-1994,
                  size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk),
                  Rich=as.vector(Rich), Temp=as.vector(Temp), sizemean=sizemean, latitudemean=latitudemean,
                  richmean=richmean, ALKmean=ALKmean, TPmean=TPmean, Tmean=Tmean)
jags2 <- jags.model("JAGSmodelbelowDCA2.txt" , data=belowDCA2, n.chains=2, n.adapt=1000)
smp3 <- coda.samples(jags2, c("mu2", "sigma.ID_site","res.sigma.ID_site","s", "lat", "TotP", "Alkal","richness", "smean", "latmean", "rmean", "alkmean", "tpmean", "tempmean","temperature","sigma.l", "sigma.y"), n.iter = 3050000, thin=3000, burnin=50000, n.chains=1, inits=inits)

smp4<-as.mcmc(smp3[1])
CI2 <- apply(smp4, 2, quantile, c(0.075, 0.925))
meansmp2 <- apply(smp4, 2, mean)

MCMCsummary(smp3, round=2)
MCMCtrace(smp3)
MCMCplot(smp3, ref_ovl = TRUE)
