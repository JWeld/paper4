#MODELS with env variables####
library(sjPlot)
library(mgcv)
library(mgcViz)
library(gratia)
library(gvlma)#checks assumptions for models

#X Dispersion, betapart reponses X####
#simple multiple regression####
ml <- lm(dispersion  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           PREC + survey_year, data = subplot.level.dat)
summary(ml)

ml <- lm(beta.SOR  ~  NH4M * NO3M + SO4SM  + latitude + longitude +
           PREC + survey_year, data = subplot.level.dat)
summary(ml)

ml <- lm(beta.SNE  ~  NTOT + SO4SM  + latitude + longitude +
           PREC + survey_year, data = plot.level.dat)

summary(lm)
#seems more going on with PC2 as response
library(visreg)
visreg(ml, "NO3M")
visreg(ml, "NH4M")
visreg(ml, "NTOT")
visreg(ml, "SO4SM")
visreg(ml, "NO3M", by = "survey_year")
visreg(ml, "NH4M", by = "survey_year")
visreg(ml, "NO3M", by = "year_scaled")
visreg(ml, "NH4M", by = "year_scaled")
visreg(ml, "NH4M", by = "NO3M")
visreg(ml, "NO3M", by = "latitude")
visreg(ml, "SO4SM", by = "survey_year")
visreg(ml, "SO4SM", by = "TEMP")


#Nmle AR regressions####

#index year! ??
library(nlme)

summary(mod<-lme(pc_dist~NTOT + SO4SM + latitude + longitude,
                 random=~1|ID_fine,
                 correlation=corCAR1(form=~year.i|ID_fine),
                 data=plot.level.dat))

summary(mod<-lme(beta.SOR~NH4M + NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_fine,
                 #correlation=corCAR1(form=~year.i|ID_fine),
                 data=drop_na(subplot.level.dat)))

summary(mod<-lme(ID_fine_disp~year.i,
                 random=~1|ID_fine,
                 correlation=corCAR1(form=~year.i|ID_fine),
                 data=subplot.level.dat))
anova(mod)
fixed.effects(mod)

plot(mod, col = as.numeric(factor(subplot.level.dat$year.i,
                                  levels = unique(subplot.level.dat$year.i))),
     pch = 16, main = "resids")


#X PCA distance to base reponses X####
#combine PCA base distance and env variables and model
library(corrplot)
combosites5 %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()
combosites5 %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() 
subplot.level.dat %>% ungroup() %>% select(where(is.numeric)) %>% drop_na() %>% cor() %>% corrplot()

#simple multiple regression####
ml <- lm(pc_dist  ~  NTOT+ SO4SM  + latitude + longitude +
           #beta.SNE +
           #beta.SIM +
           #beta.SOR +
           PREC +
           year.i, data = plot.level.dat)
summary(ml)
library(car)
vif(ml)
#attributes(alias(ml)$Complete)$dimnames[[1]] #beta.SOR is linearly dependant on SOR + SNE

ml <- lm(pc_dist_base  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           #beta.SOR +
           PREC + year.i, data = subplot.level.dat)
summary(ml)

ml <- lm(pc_dist_base  ~  NH4M + NO3M + SO4SM  + latitude + longitude +
           beta.SNE +
           beta.SIM +
           beta.SOR+
           PREC + year.i, data = plot.level.dat)
summary(ml)


glm <- glm(pc_dist_base  ~  NH4M + NO3M + SO4SM + TEMP + latitude + longitude +
             PREC + year.i, family = gaussian, data = plot.level.dat)

summary(glm)
#seems more going on with PC2 as response
library(visreg)
visreg(ml, "NO3M")
visreg(ml, "NH4M")
visreg(ml, "SO4SM")
visreg(ml, "NO3M", by = "year.i")
visreg(ml, "NH4M", by = "year.i")
visreg(ml, "NO3M", by = "year_scaled")
visreg(ml, "NH4M", by = "year_scaled")
visreg(ml, "NH4M", by = "NO3M")
visreg(ml, "NO3M", by = "latitude")
visreg(ml, "SO4SM", by = "year.i")
visreg(ml, "SO4SM", by = "TEMP")


#Nmle AR regressions####

#index year?
library(nlme)
library(lme4)
plot.level.dat$year.i <- I(plot.level.dat$survey_year-1997)
plot.level.dat$year.i <- as.numeric(plot.level.dat$year.i)

barplot(table(plot.level.dat$ID_plot))
barplot(table(plot.level.dat$latitude))
barplot(plot.level.dat$NTOT)
barplot(plot.level.dat$year.i)
barplot(table(plot.level.dat$year.i))
hist(plot.level.dat$pc_dist_base)#right skewed
hist(log(plot.level.dat$pc_dist_base))#log transf?

summary(nmod1 <- lme(pc_dist~
                       NH4M +
                       NO3M +
                       SO4SM +
                       latitude + longitude,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=plot.level.dat))

summary(nmod2 <- lme(pc_dist_base~
                       #NH4M * NO3M +
                       NTOT + 
                       SO4SM +
                       latitude + longitude +
                       year.i,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=plot.level.dat))

summary(nmod3 <- lme(log(pc_dist_base+1)~
                       #NH4M * NO3M +
                       NTOT + 
                       SO4SM +
                       latitude + longitude +
                       year.i,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=plot.level.dat))
AIC(nmod1,nmod2,nmod3)

#basic LM
summary(mod<-lm(pc_dist_base~
                  #NH4M * NO3M +
                  NTOT + 
                  SO4SM +
                  latitude + longitude +
                  year.i,
                data=plot.level.dat))

par(mfrow = c(2,2))
plot(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))
vif(mod)#car

#GLM null model
GLM <- gls(log(pc_dist_base+1)~
             #NH4M * NO3M +
             NTOT + 
             SO4SM +
             latitude + longitude +
             survey_year,
           data=subplot.level.dat,
           method = "ML")
summary(GLM)
plot(GLM)

#LME mixed
lmm0 <- lmer(pc_dist_base~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1.5 <- lmer(log(pc_dist_base)~
                 #NH4M * NO3M +
                 NTOT + 
                 SO4SM +
                 latitude +
                 longitude +
                 survey_year+
                 (1|ID_subplot),
               data=subplot.level.dat)

lmm2 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm3 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               #survey_year+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm4 <- lmer(log(pc_dist_base+1)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               survey_year +
               (1|ID_fine),
             
             data=subplot.level.dat)

AIC(GLM, lmm1, lmm1.5, lmm2,lmm3,lmm4)
plot(lmm1.5)
plot(ranef(lmm1.5)) 
# QQ plots (drawn to the same scale!)
par(mfrow = c(1,2))
lims <- c(-3.5,3.5)
qqnorm(resid(GLM, type = "pearson"),
       xlim = lims, ylim = lims,main = "GLM")
abline(0,1, col = "red", lty = 2)
qqnorm(resid(lmm1.5, type = "pearson"),
       xlim = lims, ylim = lims, main = "lmm1")
abline(0,1, col = "red", lty = 2)

mod <- gam1

summary(lmm1)
#percentage variance explained by random factor ID_plot
0.0007465/(0.0007465 + 0.0026378)#22%

summary(lmm1.5)
#percentage variance explained by random factor ID_plot
0.0009634/(0.0009634 + 0.0035281)#22%

summary(lmm2)
#percentage variance explained by random factor ID
0.002088/(0.002088 + 0.001063)#66%

anova(lmm1, lmm1.5)

qqnorm(resid(mod))#2016_SE14_1 is an outlier on 0.20+
qqline(resid(mod))
#sjPlot
plot_model(lmm1.5, type = "re", show.values = TRUE)
plot_model(lmm1.5, show.values = TRUE)


#Facet plots of model predictions####
#comboplots
(mm_plot <- ggplot(plot.level.dat, aes(x = year.i, y = pc_dist_base, colour = ID_plot)) +
   facet_wrap(~ID_plot, nrow=2) +   # a panel for each mountain range
   geom_point(alpha = 0.5) +
   theme_classic() +
   geom_line(data = cbind(plot.level.dat, pred = predict(mod)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
   theme(legend.position = "none",
         panel.spacing = unit(2, "lines"))  # adding space between panels
)

#subplot.level.dat
(mm_plot <- ggplot(subplot.level.dat, aes(x = year.i, y = pc_dist_base, colour = ID_plot)) +
    facet_wrap(~ID_plot, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(subplot.level.dat, pred = predict(mod)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)


summary(mod<-lme(pc_dist_base~NH4M * NO3M + SO4SM + latitude + longitude,
                 random=~1|ID_subplot,
                 correlation=corCAR1(form=~year.i|ID_subplot),
                 data=drop_na(subplot.level.dat)))


#GAM####
#try same model with bam
k=10
bs="ts"
gam0 <- bam(round(pc_dist*100) ~ s(NH4M) +
              s(NO3M)+
              #s(SO4SM.lag,k=k, bs=bs)+
              s(NTOT.lag,k=k, bs=bs)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              s(survey_year, bs="re")+
              #s(year.i,k=k, bs=bs)+
              #s(ID_site, bs= "re")+
              s(ID_plot, bs= "re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = poisson,
            data = subplot.level.dat)
#data = subplot.level.dat)

gam0f <- bam(round(pc_dist*100) ~ #s(NH4M) +
               #s(NO3M)+
               s(SO4SM.lag,k=k, bs=bs)+
               s(NTOT.lag,k=k, bs=bs)+
               s(PREC)+
               #s(latitude,k=k, bs=bs)+
               #s(longitude,k=k, bs=bs)+
               te(latitude, longitude,k=k, bs=bs)+
               s(survey_year,k=k, bs=bs)+
               s(ID_subplot, bs="re"),
             nthreads=8, 
             family = poisson,
             data = subplot.level.dat)
#data = subplot.level.dat)

gam1 <- bam(round(pc_dist*100) ~ #s(NH4M) +
              #s(NO3M)+
              s(SO4SM.lag,k=k, bs=bs)+
              s(NTOT.lag,k=k, bs=bs)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              s(survey_year,k=k, bs=bs)+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = poisson,
            data = subplot.level.dat)
#data = subplot.level.dat)

gam1.5 <- bam(round(pc_dist*100) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT.lag,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(survey_year,k=k, bs=bs)+
                s(ID_plot, bs="re")+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = poisson,
              data = subplot.level.dat)
#data = subplot.level.dat)

gam2 <- bam(log(pc_dist) ~ #s(NH4M) +
              #s(NO3M)+
              #s(SO4SM,k=k, bs=bs)+
              s(NTOT,k=k, bs=bs)+
              #s(beta.SOR)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              #year.i+
              s(survey_year,k=k, bs=bs)+
              s(ID_plot,bs="re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = gaussian,
            method = "REML",
            #data = all_data_dca)
            data = subplot.level.dat)
#data = plot.level.dat)

gam2.1 <- bam(round(pc_dist*100) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT.lag,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(survey_year,k=k, bs=bs)+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = poisson,
              data = subplot.level.dat)


gam2.2 <- bam(sqrt(pc_dist) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT.lag,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(survey_year,k=k, bs=bs)+
                s(ID_plot,bs="re")+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = betar,
              data = subplot.level.dat)

gam2.3 <- bam(log(pc_dist) ~ #s(NH4M) +
                #s(NO3M)+
                #s(SO4SM.lag,k=k, bs=bs)+
                s(NTOT,k=k, bs=bs)+
                #s(PREC)+
                #s(latitude,k=k, bs=bs)+
                #s(longitude,k=k, bs=bs)+
                #te(latitude, longitude,k=k, bs=bs)+
                s(survey_year,k=k, bs=bs)+
                s(ID_plot, bs="re")+
                s(ID_subplot, bs="re"),
              nthreads=8, 
              family = gaussian,
              data = subplot.level.dat)
#data = plot.level.dat)

gam3 <- bam(log(pc_dist_base) ~ s(NH4M) +
              s(NO3M)+
              #s(NTOT)+
              s(SO4SM)+
              ti(SO4SM,NO3M,NH4M)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              s(survey_year,k=k, bs=bs)+
              s(ID_plot, bs= "re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            select = TRUE,
            family = gaussian,
            data = subplot.level.dat)

gam4 <- bam(log(pc_dist) ~ #s(NH4M) +
              #s(NO3M)+
              #s(SO4SM,k=k, bs=bs)+
              #s(NTOT,k=k, bs=bs)+
              te(NTOT, SO4SM)+
              #s(beta.SOR)+
              #s(PREC)+
              #s(latitude,k=k, bs=bs)+
              #s(longitude,k=k, bs=bs)+
              #te(latitude, longitude,k=k, bs=bs)+
              #year.i+
              s(survey_year,k=k, bs=bs)+
              s(ID_plot,bs="re")+
              s(ID_subplot, bs="re"),
            nthreads=8, 
            family = gaussian,
            method = "REML",
            data = subplot.level.dat.no)

# qdat <- select(subplot.level.dat, NTOT, SO4SM, year.i, ID_plot, ID_subplot, latitude, longitude,
#                pc_dist) %>% drop_na()

b <- getViz(gam3)
summary(b)
draw(b, parametric = FALSE, residuals = TRUE)
appraise(b)
acf(residuals(b))
concurvity(b, full = FALSE)
library(sjPlot)
tab_model(b, show.ci = FALSE, show.stat = TRUE, show.r2 = FALSE,
          show.obs = FALSE, string.stat = "Chi.sq", string.est = "edf", file ="GAM1.html")#on dat 1493.023
library(dsm)#vis.concurvity

vis_concurvity(b, type = "estimate")
vis_concurvity(b, type = "worst")

AIC(gam0, gam1,gam1.5,gam2.2,gam2.3, gam3)

#X PCA sequential distance reponses X####

#lm
summary(lm1 <- lm(log(pc_dist)~
                    NTOT.lag +
                    #NH4M +
                    #NO3M +
                    SO4SM.lag +
                    #TEMP +
                    PREC +
                    #SO4SM.lag +
                    latitude +
                    longitude +
                    year.i, 
                  data=subplot.level.dat))
# data=plot.level.dat))
gvlma(lm1)
acf(resid(lm1))
mean(lm1$residuals)
#nmle
hist(subplot.level.dat$pc_dist)
hist(log(subplot.level.dat$pc_dist))#log transf?


summary(nmod1 <- lme(pc_dist~
                       NTOT.lag +
                       year.i ,
                     #NH4M +
                     #NO3M +
                     #SO4SM.lag ,
                     #latitude + longitude,
                     random=~1|ID_plot,
                     correlation=corCAR1(form=~year.i|ID_plot),
                     data=drop_na(select(plot.level.dat, - TEMP))))

summary(nmod2 <- lme(pc_dist~
                       NTOT.lag +
                       year.i + 
                       #NH4M +
                       #NO3M +
                       SO4SM.lag +
                       latitude + longitude,
                     random=~1|ID_plot/ID_subplot,
                     correlation=corCAR1(form=~year.i|ID_plot/ID_subplot),
                     data=drop_na(select(subplot.level.dat, - TEMP))))


summary(nmod3 <- lme(log10(pc_dist)~
                       NTOT.lag +
                       #year.i +
                       #NH4M * NO3M +
                       SO4SM.lag +
                       latitude + longitude,
                     random=~1|ID_plot/ID_subplot,
                     correlation=corCAR1(form=~year.i|ID_plot/ID_subplot),
                     data=drop_na(select(subplot.level.dat, - TEMP))))

summary(nmod4 <- lme(sqrt(pc_dist)~
                       NTOT.lag +
                       #year.i +
                       #NH4M * NO3M +
                       SO4SM.lag +
                       latitude + longitude,
                     random=~1|ID_plot/ID_subplot,
                     correlation=corCAR1(form=~year.i|ID_plot/ID_subplot),
                     data=drop_na(select(subplot.level.dat, - TEMP))))

AIC(nmod1,nmod2,nmod3,nmod4)
plot(nmod1)
anova(nmod1)
fixed.effects(nmod1)
acf(resid(nmod1))
nmod1
plot(nmod1, col = c(1:nlevels(as.factor(subplot.level.dat$ID_plot))), pch = 16)
vif(nmod1)

#Lmer####
#GLM null model
GLM <- gls(log(pc_dist)~
             #NH4M * NO3M +
             NTOT + 
             SO4SM +
             latitude + longitude +
             year.i,
           data=subplot.level.dat,
           # family = "poisson",
           method = "ML")
summary(GLM)
plot(GLM)

#GLM null model
GLM2 <- glm(pc_dist~
              #NH4M * NO3M +
              NTOT + 
              SO4SM +
              latitude +
              longitude +
              year.i,
            data=subplot.level.dat,
            family = "poisson")
summary(GLM2)
plot(GLM2)
plot_model(GLM2, show.values = TRUE) #sjplot

#LME mixed
lmm0 <- lmer(pc_dist~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm1.5 <- lmer(log(pc_dist)~
                 #NH4M * NO3M +
                 NTOT + 
                 #SO4SM +
                 latitude +
                 longitude +
                 year.i+
                 (1|ID_subplot),
               data=subplot.level.dat)

lmm2 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT.lag + 
               #SO4SM +
               latitude +
               longitude +
               year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm3 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               #year.i+
               (1|ID_subplot),
             data=subplot.level.dat)

lmm4 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               SO4SM +
               latitude +
               longitude +
               (1|year.i) +
               (1|ID_fine),
             
             data=subplot.level.dat)

lmm5 <- lmer(log(pc_dist)~
               #NH4M * NO3M +
               NTOT + 
               #SO4SM +
               latitude +
               longitude +
               (1|ID_plot),
             #(1|year.i),
             #(1|ID_fine),
             
             data=subplot.level.dat)

AIC(GLM, lmm1, lmm1.5, lmm2,lmm3,lmm4,lmm5)

mod <- lmm5

plot(mod)
plot(ranef(mod)) 
# QQ plots (drawn to the same scale!)
par(mfrow = c(1,2))
lims <- c(-3.5,3.5)
qqnorm(resid(GLM, type = "pearson"),
       xlim = lims, ylim = lims,main = "GLM")
abline(0,1, col = "red", lty = 2)
qqnorm(resid(mod, type = "pearson"),
       xlim = lims, ylim = lims, main = "lmm1")
abline(0,1, col = "red", lty = 2)
summary(mod)
#percentage variance explained by random factor ID_plot
0.0007465/(0.0007465 + 0.0026378)#22%

summary(lmm1.5)
#percentage variance explained by random factor ID_plot
0.0009634/(0.0009634 + 0.0035281)#22%

summary(lmm2)
#percentage variance explained by random factor ID
0.002088/(0.002088 + 0.001063)#66%

anova(lmm1, lmm1.5)

qqnorm(resid(mod))#2016_SE14_1 is an outlier on 0.20+
qqline(resid(mod))
#sjPlot
plot_model(mod, type = "re", show.values = TRUE)
plot_model(mod, show.values = TRUE)

acf(resid(mod))

library(ggeffects)
ggpredict(mod, terms = c("NTOT")) %>% plot()
ggpredict(mod, terms = c("NTOT", "ID_plot"), type = "re") %>% plot() +
  theme(legend.position = "bottom")

#MCMCglmm####
library(MCMCglmm)
library(MCMCvis)
mod <- MCMCglmm(round(pc_dist*100) ~ NH4M.lag + year.i, random = ~ID_subplot,
                family = "poisson", data = subplot.level.dat)

summary(mod)
plot(mod$Sol)
plot(mod$VCV)
MCMCplot(mod$Sol)
MCMCplot(mod$VCV)


#glmmTMB####
library(glmmTMB)
library(broom)
library(DHARMa)
mod0 <- glmmTMB(log(pc_dist) ~ NTOT + (1|year.i) + #(1|year.i/ID_subplot) +
                  (1|ID_site/ID_plot/ID_subplot),
                #(1|ID_fine),
                family = "gaussian", data = subplot.level.dat)

mod0 <- glmmTMB(log(pc_dist) ~ NTOT + SO4SM +(1|year.i) + #(1|year.i/ID_subplot) +
                  (1|ID_site)+ (1|ID_plot)+ (1|ID_subplot),
                #(1|ID_fine),
                family = "gaussian", data = subplot.level.dat)

mod1 <- glmmTMB(round(pc_dist*100)~ NTOT +
                  #SO4SM +
                  year.i+
                  (1|year.i) + #(1|year.i/ID_subplot) +
                  (1|ID_site)+
                  (1|ID_plot)+
                  (1|ID_subplot),
                #(1|ID_fine),
                family = "poisson", data = subplot.level.dat)
mod <- mod1
summary(mod)
res = simulateResiduals(mod)
plot(res, rank = T)
acf(resid(mod))

mod1 <- glmmTMB(round(pc_dist*100) ~ NTOT + year.i +(1|ID_subplot) +
                  (1|ID_plot),
                family = "poisson", data = subplot.level.dat)


mod2 <- glmmTMB(log(pc_dist) ~ NTOT + year.i +
                  #beta.SIM +
                  (1|ID_subplot) +
                  (1|ID_plot),
                family = "gaussian", data = subplot.level.dat.no)

mod3 <- glmmTMB(log(pc_dist) ~ NTOT +(1|year.i) +(1|ID_subplot) +
                  (1|ID_plot),
                family = "gaussian", data = subplot.level.dat.no)

mod4 <- glmmTMB(sqrt(pc_dist) ~ NTOT + (1|year.i) +(1|ID_subplot) +
                  (1|ID_plot),
                family = "gaussian", data = subplot.level.dat)

mod5 <- glmmTMB(pc_dist ~ NTOT + year.i +(1|ID_subplot) +
                  (1|ID_plot), beta_family(),
                data = subplot.level.dat)

mod6 <- glmmTMB(pc_dist ~ NTOT +(1|ID_subplot/year.i) +
                  (1|ID_plot/year.i), beta_family(),
                data = subplot.level.dat)

AIC(mod1, mod2,mod3,mod4,mod5,mod6)

BIC(mod1, mod2,mod3,mod4,mod5,mod6,GLM, lmm1, lmm1.5,
    lmm2,lmm3,lmm4,lmm5,gam0, gam0f, gam1,gam1.5, gam2,gam2.2,
    gam2.3, gam3) %>% arrange(BIC)

mod <- mod2
summary(mod)
res = simulateResiduals(mod)
#res = simulateResiduals(mod, form = subplot.level.dat$year.i)
plot(res, rank = T)

out <- boxplot.stats(sqrt(subplot.level.dat$pc_dist))$out
out_ind <- which(sqrt(subplot.level.dat$pc_dist) %in% c(out))
out_ind

out <- boxplot.stats(log(subplot.level.dat$pc_dist))$out
out_ind <- which(log(subplot.level.dat$pc_dist) %in% c(out))
out_ind
outliers <- subplot.level.dat[out_ind, ]

subplot.level.dat.no <- subplot.level.dat[-out_ind, ]

acf(resid(mod))

testDispersion(res)
testUniformity(res)
testOutliers(res)
testZeroInflation(res)
testTemporalAutocorrelation(res.t)
testSpatialAutocorrelation(res, x= subplot.level.dat$longitude, y = subplot.level.dat$latitude)

res.t <- recalculateResiduals(res, group = subplot.level.dat$year.i)
plot(res, quantreg = FALSE)



#MAPS#####

library(ggplot2)  # ggplot() fortify()
library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)  # readOGR() spTransform()
library(raster)  # intersect()
library(ggsn)  # north2() scalebar()
library(rworldmap)
world <- getMap(resolution = "low")

(site_check <- ggplot(all_data, mapping = aes(x = longitude, y = latitude)) + 
    geom_point(alpha = 0.5))

clipper_europe <- as(extent(3, 32, 50, 72), "SpatialPolygons")

proj4string(clipper_europe) <- CRS(proj4string(world))

world_clip <- raster::intersect(world, clipper_europe)

world_clip_f <- fortify(world_clip)

(site_map <- ggplot() + 
    geom_polygon(data = world_clip_f, 
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "black") + 
    geom_point(colour = "black", size = 2.5, shape = 21, fill = "grey70",
               stroke = 1,
               aes(x = longitude, y = latitude),
               data = all_data) +
    # labs(shape="Forest type")+
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") + 
    coord_quickmap())

(site_map <- ggplot() + 
    geom_polygon(data = world_clip_f, 
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "black") + 
    geom_point(colour = "black", size = 1.5,
               aes(x = longitude, y = latitude),
               data = all_data) +
    # labs(shape="Forest type")+
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") + 
    coord_quickmap())


