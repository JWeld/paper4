# GAMS --------------------------------------------------------------------
library(tidyverse)
library(mgcv)
library(parallel)
library(mgcViz)
#on p values and gam summary/ gam check difference!####
#The p values relate to two entirely different tests: in summary.gam the p values are of the null hypothesis
#of a zero effect of the indicated spline. There values relate to the F statistic in the table produced by summary.gam,
#in gam.check the p values are for the test of the null hypothesis that the basis dimension used is of sufficient size.
#I.e. these p values relate to the value labelled k-index in the table produced by gam.check

# Spatial autocorrelation
# Two ways within a gam/gamm:
#   - including interaction of latitude/longitude: fit<-gam(response~….+s(lat,long),…)
#   - Covariance: fit<-gamm(response~….,correlation=corGaus(1,form=~lat+long))

#define data used in models####
dat <- filter(all_data, beta.SOR > 0)
k <- 20
bs <- "tp"

#subset for quick calc and/or test and train####
data <- dat
ind <- sample(2, nrow(data), replace = TRUE, prob=c(0.3, 0.7))
train <- data[ind == 1,]
test <- data[ind == 2,]

dat <- train#
dat <- data

### function to print model summaries
printresults<-function(gobject){
  edf=sum(gobject$edf)
  nobs=length(gobject$y)
  BIC<--2*logLik(gobject)+edf*log(nobs)
  results<-c(sum(gobject$edf), gobject$deviance,summary(gobject)$r.sq,
             gobject$aic,BIC)
  results<-round(results,2)
  names(results)<-c("edf","deviance","adj.R-sq","AIC","BIC")
  results
}


gam1 <- gam(dispersion ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+  s(SO4SM, k = k,bs = bs)+
                te(latitude, longitude, year.i,d = c(2, 1), bs = c("cs", "cs"), k = c(10, 10))+
                s(ID_plot, bs="re"),nthreads=12, data = dat, family = "betar") #k5

gam1 <- gam(PC1 ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+  s(SO4SM, k = k,bs = bs)+
              te(latitude, longitude, year.i,d = c(2, 1), bs = c("cs", "cs"), k = c(10, 10))+
              s(ID_plot, bs="re"),nthreads=12, data = dat) #k5

gam2 <- gam(beta.SNE ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+  s(SO4SM, k = k,bs = bs)+
              te(latitude, longitude, survey_year,d = c(2, 1), bs = c("cs", "cs"), k = c(10, 10))+
              s(ID_plot, bs="re"),nthreads=12, data = dat) #k5

gam2 <- bam(dispersion ~ #s(NH4M, k = k,bs = bs) +
              s(NO3M, k = k,bs = bs)+
              #s(SO4SM, k = k,bs = bs)+
              #s(ID_plot, bs="re")+
              s(latitude)+
              s(longitude)+
              s(survey_year),
             # te(latitude, longitude, survey_year),
            family=betar(link="logit"),
            nthreads=12, data = dat)

summary(gam2)

library(betareg)
library(lmtest)
fit <- betareg(dispersion ~ NH4M + NO3M + SO4SM + PREC +
                 latitude + longitude + survey_year, data = filter(dat, dispersion > 0))

fit <- betareg(dispersion ~ NH4M + NO3M + SO4SM + PREC + (1| ID_plot) +
                 latitude + longitude + survey_year, data = filter(dat, dispersion > 0))

sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
           function(x) logLik(update(fit, link = x)))

library(glmmTMB)
fit <- glmmTMB(dispersion ~ NH4M + NO3M + SO4SM + PREC +
                 latitude + longitude + survey_year,
               family= beta_family(link = "logit"),data = filter(dat, dispersion > 0))

fit <- glmmTMB(beta.SOR ~ NH4M + NO3M + SO4SM + PREC +
                 latitude + longitude + survey_year,
               family= beta_family(link = "logit"),data = filter(dat, beta.SOR > 0))

plot(fit$fitted)
summary(fit)
coeftest(fit)


gamm1 <- gamm(beta.SOR ~ s(NH4M, k = k,bs = bs) +
                  s(NO3M, k = k,bs = bs)+
                  s(SO4SM, k = k,bs = bs)+
                
                  #ti(NH4M, NO3M, k = k, bs = bs)+
                  #s(ID_plot, bs="re"),
                s(latitude) +
                s(longitude),
                correlation = corAR1(form = ~ 1| survey_year),
              family=betar(link="logit"), data = dat)

# gamm1.1.b <- gamm(N.y ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+
#                     te(NH4M, NO3M, k = k, bs = bs) + s(survey_year),
#                     #correlation = corAR1(form = ~ 1| survey_year),
#                     correlation = corGaus(1, form = ~ latitude + longitude),
#                     data = dat) 

gamm1.1.b <- gamm(N.y ~ s(NH4M) + s(NO3M)+
                    ti(NH4M, NO3M) + te(latitude, longitude)+
                    s(ID_site, bs="re"),
                  correlation = corAR1(form = ~ 1| survey_year), nthreads=12,
                  data = dat1) 



gamm1.2.b <- gamm(N.y ~ ti(NH4M, k = k,bs = bs) + ti(NO3M, k = k,bs = bs)+
                    ti(NH4M, NO3M, k = k, bs = bs) + te(latitude, longitude,k = k,bs = bs)+
                    s(ID_site, bs="re"),
                  correlation = corAR1(form = ~ 1| survey_year), nthreads=4,
                  data = dat) 

gamm1.3.b <- gamm(N.y ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+
                    te(NH4M, NO3M, k = k,bs = bs) + te(latitude, longitude,k = k,bs = bs),
                  s(ID_site/plot.x, bs="re")+
                    correlation = corAR1(form = ~ 1| survey_year), niterPQL = 20,
                  #correlation = corGaus(1, form = ~ latitude + longitude),
                  data = dat) 

gamm1.4.b <- gamm(N.y ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+
                    te(NH4M, NO3M, k = k,bs = bs) +
                    s(survey_year/ID_site/plot.x, bs="re"),
                  data = dat) 

gamm1.5.b <- gamm(N.y ~ s(NH4M, k = k,bs = bs) + s(NO3M, k = k,bs = bs)+
                    te(NH4M, NO3M, k = k,bs = bs) +
                    s(survey_year/ID_site/plot.x, bs="re"),
                  correlation = corAR1(form = ~ 1| survey_year),
                  data = dat) 


AIC(gamm1.0.b$lme, gamm1.1.b$lme, gamm1.2.b$lme, gamm1.3.b$lme, gamm1.4.b$lme, gamm1.5.b$lme)
AIC(gam1.b,gamm1.b$gam, gam2.b, gam3.b, gam4.b, gam5.b, gam6.b, gam7.b, gam8.b, gam9.b)


fit <- gamm1$gam
fit <- gam2

summary(fit)
#gam.check(fit)
fitb <- getViz(fit)
check(fitb,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

o <- plot( sm(fitb, 1) )
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


plot(dat$beta.SOR, sign(residuals(fit)) + rnorm(nrow(dat), 0, 0.05), col = alpha(1, 0.4), pch = 16,
     ylab = "Residuals sign", xlab = "x0")

plot.gam(fit)

concurvity(fit, full = TRUE)

layout(matrix(1:2, ncol = 2))

acf(resid(fit), main = "ACF", type = "correlation")
acf(resid(fit), main = "ACF", type = "partial")

layout(1)
#previous best N.y ~ ti(NH4M) + ti(NO3M) + ti(NH4M, NO3M) + te(longitude, latitude, survey_year, L.y) + s(sum_canopy)
AIC(gamm1.b$lme, gamm1.1.b$lme, gamm1.2.b$lme, gamm1.3.b$lme)



library(mgcv)
library(voxel)
library(tidyverse)
library(gridExtra)
library(sp)
library(gratia)
vars <- c("NH4M", "NO3M")

map(vars, function(x){
  p <- plotGAM(gam1,smooth.cov = "NO3M") #plot customization goes here
  g <- ggplotGrob(p)
}) %>%
  {grid.arrange(grobs = (.), ncol = 2, nrow = 2)} 

#variogram for spatial autocorrelation
dat2 <- dat
dat2$resid.gam <- residuals(gam1)
plot(dat2$resid.gam, dat2$survey_year)
plot(dat2$resid.gam, dat2$longitude)
plot(dat2$resid.gam)



summary(fit)
plot(fit, pages = 1, all.terms = TRUE, shade = TRUE, shift = coef(fit)[1])
gam.check(fit)


#BAMS####
dat <- all_data
#N####
bam1<- bam(dispersion ~ s(NH4M) + s(NO3M)+
             ti(NH4M, NO3M) +
             te(latitude, longitude, survey_year)+
             #s(sum_canopy) + 
             s(SO4SM)+
             s(PREC)+
             s(ID_site, bs="re"),
           data = dat) 

bam2 <- bam(dispersion ~ s(NH4M) + s(NO3M)+
              #ti(NH4M, NO3M) +
              te(latitude, longitude, survey_year)+
              #s(sum_canopy) + 
              s(SO4SM)+
              s(PREC),
            data = dat) 

bam3 <- bam(dispersion ~ s(NH4M) + s(NO3M)+
              ti(NH4M, NO3M) +
              te(latitude, longitude, survey_year)+
              #s(sum_canopy) + 
              s(SO4SM)+
              s(PREC)+
              s(ID_plot, bs="re"),
            data = dat) 

bam4 <- bam(dispersion ~ s(NH4M) + s(NO3M)+
              ti(NH4M, NO3M) +
              te(latitude, longitude, survey_year)+
              #s(sum_canopy) + 
              s(SO4SM)+
              s(PREC)+
              s(ID_subplot),
            data = dat) 

bam5 <- bam(dispersion ~ s(NH4M) + s(NO3M)+
              ti(NH4M, NO3M) +
              te(latitude, longitude, survey_year)+
              #s(sum_canopy) + 
              s(SO4SM),
            #s(PREC),
            #s(country, bs="re"),
            data = dat) 
summary(bam1)
b <- getViz(bam1)
print(plot(b, allTerms = T), pages = 1) 



summary(bam2)
b <- getViz(bam2)
print(plot(b, allTerms = T), pages = 1) 
gam.check(b,pch=19,cex=.3)

AIC(bam1,bam2,bam3,bam5)

# Make predictions
predictions <- bam1 %>% predict(test)
predictions <- exp(predictions)
plot(predictions, test$dispersion)



