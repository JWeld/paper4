require(rjags)
require(MCMCglmm)

n.lakes<-57
n.obs<-as.vector(table(factor(InvBel$Lake)))
sizemean<-as.vector(tapply(log(InvBel$Size), factor(InvBel$Lake), mean))
latitudemean <- as.vector(tapply(scale(InvBel$Latfix), factor(InvBel$Lake), mean))
richmean <- as.vector(tapply(scale(InvBel$Richness), factor(InvBel$Lake), mean))
ALKmean <- as.vector(tapply(log(InvBel$`Alk/Acid`), factor(InvBel$Lake), mean))
TPmean <- as.vector(tapply(log(InvBel$`Tot-P`), factor(InvBel$Lake), mean))
#Trecord <- as.vector(tapply(scale(InvBel$Trecord), factor(InvBel$Lake), mean))

size <-log(InvBel$Size)
latitude <- scale(InvBel$Latfix)
TP <- log(InvBel$`Tot-P`)
Alk <- log(InvBel$`Alk/Acid`)
Rich <- scale(InvBel$Richness)


model <- "model{
for(i in 1:n) {
y[i]~dnorm(y.hat[i], tau.lake[lake[i]])
y.hat[i]<- mu + l[lake[i]] + ye[year[i]] + size[i]*s + latitude[i]*lat + Rich[i]*richness + TP[i]*TotP + Alk[i]*Alkal 
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
mu.Lake.sigma[j] <- mu2 + sizemean[j]*smean + latitudemean[j]*latmean + richmean[j]*rmean + ALKmean[j]*alkmean + TPmean[j]*tpmean 
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
}"

##DCA1
writeLines(model, "JAGSmodelbelowDCA1.txt")
belowDCA1 <- list(y =  InvBel$DCA1,  n = length(InvBel$DCA1), L=n.lakes, Y=23, lake=rep(1:n.lakes, n.obs), year=as.numeric(levels(InvBel$Year))[InvBel$Year]-1994, size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk), Rich=as.vector(Rich), sizemean=sizemean, latitudemean=latitudemean, richmean=richmean, ALKmean=ALKmean, TPmean=TPmean)
jags1 <- jags.model("JAGSmodelbelowDCA1.txt" , data=belowDCA1, n.chains=1, n.adapt=1000)
smp1 <- coda.samples(jags1, c("mu2", "sigma.lake","res.sigma.lake","s", "lat", "TotP", "Alkal","richness", "smean", "latmean", "rmean", "alkmean", "tpmean", "sigma.l", "sigma.y"), n.iter = 30500, thin=30, burnin=50, n.chains=1, inits=inits)

smp2<-as.mcmc(smp1[1])
CI <- apply(smp2, 2, quantile, c(0.075, 0.925))
meansmp <- apply(smp2, 2, mean)

##DCA2
writeLines(model, "JAGSmodelbelowDCA2.txt")
belowDCA2 <- list(y = InvBel$DCA2,  n = length(InvBel$DCA2), L=n.lakes, Y=23, lake=rep(1:n.lakes, n.obs), year=as.numeric(levels(InvBel$Year))[InvBel$Year]-1994, size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk), Rich=as.vector(Rich), sizemean=sizemean, latitudemean=latitudemean, richmean=richmean, ALKmean=ALKmean, TPmean=TPmean)
jags2 <- jags.model("JAGSmodelbelowDCA2.txt" , data=belowDCA2, n.chains=1, n.adapt=1000)
smp3 <- coda.samples(jags2, c("mu2", "sigma.lake","res.sigma.lake","s", "lat", "TotP", "Alkal","richness", "smean", "latmean", "rmean", "alkmean", "tpmean", "sigma.l", "sigma.y"), n.iter = 30500, thin=30, burnin=50, n.chains=1, inits=inits)

smp4<-as.mcmc(smp3[1])
CI2 <- apply(smp4, 2, quantile, c(0.075, 0.925))
meansmp2 <- apply(smp4, 2, mean)

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


write.csv2(cbind(paste(unique(factor(InvBel$Lake))), 
                 paste(round(meansmp[12:68]^2,2), paste=" (",round(CI[1,12:68]^2,2), 
                       paste=", ", round(CI[2,12:68]^2,2), paste=")", sep = ""), 
                 paste(round(meansmp2[12:68]^2,2), paste=" (",round(CI2[1,12:68]^2,2), 
                       paste=", ", round(CI2[2,12:62]^2,2), paste=")", sep = "")), 
           "HRsBELOW.csv")

##Richness DCA1
richmeanBL <- as.vector(tapply(scale(InvBel$Richness), factor(InvBel$Lake), mean))
rBLN <- seq(min(richmeanBL), max(richmeanBL), by=0.01)
yR <- matrix(NA, 1016, length(rBLN))
for (i in 1:1016) {
  yR[i,]=exp(smp2[i,6] + smp2[i,9]*rBLN)^2
}
my <- apply(log(yR),2,mean)
ly <- apply(log(yR),2,quantile,0.025)
hy <- apply(log(yR),2,quantile,0.975)

##Richness DCA2
richmeanBL <- as.vector(tapply(scale(InvBel$Richness), factor(InvBel$Lake), mean))
rBLN <- seq(min(richmeanBL), max(richmeanBL), by=0.01)
yRR <- matrix(NA, 1016, length(rBLN))
for (i in 1:1016) {
  yRR[i,]=exp(smp4[i,6] + smp4[i,9]*rBLN)^2
}
myR <- apply(log(yRR),2,mean)
lyR <- apply(log(yRR),2,quantile,0.025)
hyR <- apply(log(yRR),2,quantile,0.975)


##Alkalinity DCA1
AlkmeanBL <-  as.vector(tapply(log(InvBel$`Alk/Acid`), factor(InvBel$Lake), mean))
AlkBLN <- seq(min(AlkmeanBL), max(AlkmeanBL), by=0.01)
AlkAB2 <- matrix(NA, 1016, length(AlkBLN))
for (i in 1:1016) {
  AlkAB2[i,]=exp(smp2[i,6] + smp2[i,3]*AlkBLN)^2
}
mx3 <- apply(log(AlkAB2),2,mean)
lx3 <- apply(log(AlkAB2),2,quantile,0.025)
hx3 <- apply(log(AlkAB2),2,quantile,0.975)

##TP DCA1
TPmeanBL <-  as.vector(tapply(log(InvBel$`Tot-P`), factor(InvBel$Lake), mean))
TPBLN <- seq(min(TPmeanBL), max(TPmeanBL), by=0.01)
TPAB2 <- matrix(NA, 1016, length(TPBLN))
for (i in 1:1016) {
  TPAB2[i,]=exp(smp2[i,6] + smp2[i,71]*TPBLN)^2
}
mx4 <- apply(log(TPAB2),2,mean)
lx4 <- apply(log(TPAB2),2,quantile,0.025)
hx4 <- apply(log(TPAB2),2,quantile,0.975)

##Size DCA2
SmeanBL <-  as.vector(tapply(log(InvBel$Size), factor(InvBel$Lake), mean))
SBLN <- seq(min(SmeanBL), max(SmeanBL), by=0.01)
size2 <- matrix(NA, 1016, length(SBLN))
for (i in 1:1016) {
  size2[i,]=exp(smp4[i,6] + smp4[i,70]*SBLN)^2
}
ms2 <- apply(log(size2),2,mean)
ls2 <- apply(log(size2),2,quantile,0.025)
hs2 <- apply(log(size2),2,quantile,0.975)

##Temp DCA1
TempmeanBL <-  as.vector(tapply(log(InvBel$Vattentemperatur), factor(InvBel$Lake), mean))
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


varDCA1BELOW <- tapply(InvBel$DCA1, factor(InvBel$Lake), var)
varDCA2BELOW <- tapply(InvBel$DCA2, factor(InvBel$Lake), var)

RichBELOW <- scale(InvBel$Richness)
AlkBELOW <- log(InvBel$`Alk/Acid`)
SizeBELOW <- scale(InvBel$Size)
