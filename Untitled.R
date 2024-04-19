hannah <- as.data.frame(matrix(rnorm(4*100, 5, 2), ncol=4))
colnames(hannah) <- c("food_quality", "beer_quality", "family_dist","study_vis" )
hannah$hannahs <- rbinom(100,1,.2)
mod <- stan_glm(hannahs ~ food_beer_qual + family_dist + study_vis,
                   chains = 4, cores = 4, iter= 6000, data = hannah, family = "binomial", adapt_delta = 0.95)
summary(mod)
launch_shinystan(mod)


library(MASS)
library(tidyverse)
library(GGally)
set.seed(5)
# create the variance covariance matrix
sigma<-rbind(c(1,0.8,0.7), c(0.8,1, 0.2), c(0.7,0.2,1))
# create the mean vector
mu<-c(7, 5, 6) 
# generate the multivariate normal distribution
df<-as.data.frame(mvrnorm(n=100, mu=mu, Sigma=sigma))
ggpairs(df)

df<-df%>%mutate(hannahs = ifelse(V1>median(V1), sample(c(0,1),n(), replace = TRUE, p=c(0.25, 0.75)) ,
                                       sample(c(0,1),n(), replace = TRUE, p=c(0.75, 0.25))))
colnames(df) <- c("food_beer_qual", "study_vis", "family_dist","hannahs" )
mod <- stan_glm(hannahs ~ food_beer_qual + family_dist + study_vis,
                chains = 4, cores = 4, iter= 6000, data = df, family = "binomial", adapt_delta = 0.95)
summary(mod)
launch_shinystan(mod)
pp_check(mod, "dist", nreps=30)
