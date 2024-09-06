library(gllvm)
library(tidyverse)

#Reduced size version for testing, too slow with full dataset
#Use 20% of dataset as training set
ordispe <- ordispe[, colSums(ordispe != 0) > 0]

sample <- sample(c(TRUE, FALSE), nrow(ordispe), replace=TRUE, prob=c(0.8,0.2))
sample.l  <- ordispe[sample, ]
sample.s   <- ordispe[!sample, ]

sample.s <- sample.s[, colSums(sample.s != 0) > 0]
sample.l <- sample.s[, colSums(sample.l != 0) > 0]

fit_base <- gllvm(sample.s, num.lv = 2, starting.val ="zero", family = "gaussian")
summary(fit_base)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(fit_base, var.colors = 1)
#not looking good! Other family values...
fit_base <- gllvm(sample.s, num.lv = 2, starting.val ="zero", family = "negative.binomial")
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(fit_base, var.colors = 1)

par(mfrow = c(1,1))
ordiplot(fit_base, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(fit_base, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)

y <- as.matrix(ordispe)
X <- as.matrix(ordi_extras)

#try with full dataset
fit_base <- gllvm(ordispe, num.lv = 2, starting.val ="zero", family = "negative.binomial")
summary(fit_base)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(fit_base, var.colors = 1)
par(mfrow = c(1,1))
ordiplot(fit_base, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(fit_base, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)


# Fit model with environmental covariates
fit <- gllvm(y, X, formula = ~ NH4M + NO3M, num.lv = 2, starting.val ="zero",
             family = "negative.binomial")
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(fit_base, var.colors = 1)
par(mfrow = c(1,1))
ordiplot(fit_base, biplot = TRUE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Biplot")
ordiplot(fit_base, biplot = FALSE, ind.spp = 15, xlim = c(-3, 3), ylim = c(-3, 3), 
         main = "Ordination plot", predict.region = TRUE)
coefplot(fit)


