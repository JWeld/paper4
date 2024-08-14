require(rjags)
require(MCMCglmm)
library(MCMCvis)

#must first read in InvBel data

n.lakes<-57
n.obs<-as.vector(table(factor(InvBel$Lake)))
sizemean<-as.vector(tapply(log(InvBel$Size), factor(InvBel$Lake), mean))
latitudemean <- as.vector(tapply(scale(InvBel$Latfix), factor(InvBel$Lake), mean))
richmean <- as.vector(tapply(scale(InvBel$Richness), factor(InvBel$Lake), mean))
ALKmean <- as.vector(tapply(log(InvBel$`Alk/Acid`), factor(InvBel$Lake), mean))
TPmean <- as.vector(tapply(log(InvBel$`Tot-P`), factor(InvBel$Lake), mean))
Tmean <- as.vector(tapply(scale(InvBel$WaterTemp), factor(InvBel$Lake), mean))

n.plots<-length(unique(dat$ID_site))#184
n.obs<-as.vector(table(dat$ID_site))
canopymean<-as.vector(tapply(dat$sum_canopy, factor(dat$ID_site), mean))
latitudemean <- as.vector(tapply(scale(dat$latitude), factor(dat$ID_site), mean))
precipmean <- as.vector(tapply(scale(dat$mean_precip), factor(dat$ID_site), mean))
NH4mean <- as.vector(tapply(dat$n_nh4, factor(dat$ID_site), mean))
NO3mean <- as.vector(tapply(dat$n_no3, factor(dat$ID_site), mean))
Tmean <- as.vector(tapply(scale(dat$mean_temp), factor(dat$ID_site), mean))

size <-log(InvBel$Size)
latitude <- scale(InvBel$Latfix)
TP <- log(InvBel$`Tot-P`)
Alk <- log(InvBel$`Alk/Acid`)
Rich <- scale(InvBel$Richness)
Temp <- scale(InvBel$WaterTemp)

canopy <-dat$sum_canopy
latitude <- scale(dat$latitude)
NO3 <- dat$n_no3
NH4 <- dat$n_nh4
precip <- scale(dat$mean_precip)
temp <- scale(dat$mean_temp)
model <- "model{
for(i in 1:n) {
y[i]~dnorm(y.hat[i], tau.lake[lake[i]])
y.hat[i]<- mu + l[lake[i]] + ye[year[i]] + size[i]*s + latitude[i]*lat +
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
tau.lake[j]<-pow(sigma.lake[j], -2)
sigma.lake[j]~dlnorm(mu.Lake.sigma[j], res.tau.lake)
mu.Lake.sigma[j] <- mu2 + sizemean[j]*smean + latitudemean[j]*latmean +
richmean[j]*rmean + ALKmean[j]*alkmean + TPmean[j]*tpmean + Tmean[j]*tempmean
}

res.tau.lake <- pow(res.sigma.lake, -2)
res.sigma.lake ~ dunif(0, 10)
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
belowDCA1 <- list(y =  InvBel$DCA1,  n = length(InvBel$DCA1), L=n.lakes, Y=23,
                  lake=rep(1:n.lakes, n.obs), year=as.numeric(levels(as.factor(InvBel$Year)))[as.factor(InvBel$Year)]-1994,
                  size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk),
                  Rich=as.vector(Rich), Temp=as.vector(Temp), sizemean=sizemean, latitudemean=latitudemean,
                  richmean=richmean, ALKmean=ALKmean, TPmean=TPmean, Tmean=Tmean)
jags1 <- jags.model("JAGSmodelbelowDCA1.txt" , data=belowDCA1, n.chains=2, n.adapt=1000)
smp1 <- coda.samples(jags1, c("mu2", "sigma.lake","res.sigma.lake","s", "lat",
                              "TotP", "Alkal","richness", "smean", "latmean", "rmean",
                              "alkmean", "tpmean", "tempmean","temperature","sigma.l", "sigma.y"),
                     n.iter = 3050, thin=3, burnin=50, n.chains=2, inits=inits)

smp2<-as.mcmc(smp1[1])
CI <- apply(smp2, 2, quantile, c(0.075, 0.925))
meansmp <- apply(smp2, 2, mean)
MCMCplot(smp1)

##DCA2
writeLines(model, "JAGSmodelbelowDCA2.txt")
belowDCA2 <- list(y = InvBel$DCA2,  n = length(InvBel$DCA2), L=n.lakes, Y=23, lake=rep(1:n.lakes, n.obs),
                  year=as.numeric(levels(InvBel$Year))[InvBel$Year]-1994, size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk), Rich=as.vector(Rich), Temp=as.vector(Temp), sizemean=sizemean, latitudemean=latitudemean, richmean=richmean, ALKmean=ALKmean, TPmean=TPmean, Tmean=Tmean)
jags2 <- jags.model("JAGSmodelbelowDCA2.txt" , data=belowDCA2, n.chains=2, n.adapt=1000)
smp3 <- coda.samples(jags2, c("mu2", "sigma.lake","res.sigma.lake","s", "lat", "TotP", "Alkal","richness", "smean", "latmean", "rmean", "alkmean", "tpmean", "tempmean","temperature","sigma.l", "sigma.y"), n.iter = 3050000, thin=3000, burnin=50000, n.chains=1, inits=inits)

smp4<-as.mcmc(smp3[1])
CI2 <- apply(smp4, 2, quantile, c(0.075, 0.925))
meansmp2 <- apply(smp4, 2, mean)

