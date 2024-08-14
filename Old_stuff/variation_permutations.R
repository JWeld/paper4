require(rjags)
require(MCMCglmm)

n.lakes<-57
n.obs<-as.vector(table(factor(BIGud$Lake)))
sizemean<-as.vector(tapply(log(BIGud$Size), factor(BIGud$Lake), mean))
latitudemean <- as.vector(tapply(scale(BIGud$Latfix), factor(BIGud$Lake), mean))
richmean <- as.vector(tapply(scale(BIGud$Richness), factor(BIGud$Lake), mean))
ALKmean <- as.vector(tapply(scale(BIGud$`Alk/Acid`), factor(BIGud$Lake), mean))
TPmean <- as.vector(tapply(log(BIGud$`Tot-P`), factor(BIGud$Lake), mean))

size <-log(BIGud$Size)
latitude <- scale(BIGud$Latfix)
TP <- log(BIGud$`Tot-P`)
Alk <- scale(BIGud$`Alk/Acid`)
Rich <- scale(BIGud$Richness)

model <- "model{
for(i in 1:n) {
y[i]~dnorm(y.hat[i], tau.lake[lake[i]])
y.hat[i]<- mu + l[lake[i]] + size[i]*s + latitude[i]*lat + Rich[i]*richness + TP[i]*TotP + Alk[i]*Alkal 
}

mu ~ dnorm (0, .0001)

for (j in 1:L){
l[j] ~ dnorm (0, tau.l)
}

for(j in 1:L){
tau.lake[j]<-pow(sigma.lake[j], -2)
sigma.lake[j]~dlnorm(mu.Lake.sigma[j], res.tau.lake)
mu.Lake.sigma[j] <- mu2 + sizemean[j]*smean + latitudemean[j]*latmean + richmean[j]*rmean + ALKmean[j]*alkmean + TPmean[j]* tpmean
}

res.tau.lake <- pow(res.sigma.lake, -2)
res.sigma.lake ~ dunif(0, 10)
tau.l <- pow(sigma.l, -2)
sigma.l ~ dunif(0, 10)
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

writeLines(model, "C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Scripts/JAGSmodelbelowDCA1.txt")


###Null distribution
n.sim<-10
null.res<-rep(NA, n.sim)

for(i in 1:n.sim){
  dat <- list(y = sample(BIGud$DCA1, length(BIGud$DCA1), replace=FALSE), n = length(BIGud$DCA1), L=n.lakes, lake=rep(1:n.lakes, times=n.obs), size=as.vector(size), latitude=as.vector(latitude), TP=as.vector(TP), Alk=as.vector(Alk), Rich=as.vector(Rich), sizemean=sizemean, latitudemean=latitudemean, richmean=richmean, ALKmean=ALKmean, TPmean=TPmean)
  jagsmodel <- jags.model("C:/Users/hhen0001/Documents/Stability_BMI_56Lakes/Scripts/JAGSmodelbelowDCA1.txt", data=dat, n.chains=1, n.adapt=1000)
  samples.null <- coda.samples(jagsmodel, c("res.sigma.lake"), n.iter = 3000, n.chains=2, inits=inits)
  null.res[i]<-mean(unlist(samples.null))
}

hist(null.res)

p <- null.res>meansmp[7]
table(p)


  