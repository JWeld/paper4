#Random Forest

#random forest####
library(randomForest)
library(randomForestExplainer)
library(pdp)
library(vip)
library(caret)
library(visreg)
#define data for models
df <- subplot.level.dat
#df <- combosites5 %>% ungroup()
df$TEMP <- NULL
df <- df %>% drop_na()
#N
RFmod <- randomForest(pc_dist_base ~ NH4M+
                                NO3M+
                                latitude +
                                longitude +
                                #NTOT  +
                                SO4SM +
                                PREC +
                                #TEMP +
                                #beta.SNE +
                                year.i,importance = TRUE,
                                mtry = 2, data = df)

train(pc_dist_base ~ 
        NH4M+
        NO3M+
        latitude +
        longitude +
        #NTOT  +
        SO4SM +
        PREC +
        #TEMP +
        #beta.SNE +
        year.i, data = df, method ="rf") #optimal 2

RFmod
importance(RFmod)
varImpPlot(RFmod) 
visreg(RFmod)







#beta.SOR
RFmod <- randomForest(pc_dist_base ~ latitude + longitude + NTOT +
                        SO4SM + PREC +
                        year.i, mtry = 2, data = df)

RFmod <- randomForest(pc_dist_base ~ latitude + longitude + NH4M + NO3M + 
                      SO4SM + PREC + 
                        survey_year, mtry = 2, data = df)

train(pc_dist_base ~ latitude + longitude + NH4M + NO3M +
        SO4SM + PREC + 
        survey_year, data = dat, method ="rf") #optimal 2

RFmod
importance(RFmod)
importance(RFmod,type=1)

varImpPlot(RFmod) 
visreg(RFmod)

fit <- RFmod
fit
varImpPlot(fit)
plot(fit)
visreg(fit)



#plot_importance_ggpairs(fit)
vip(fit, horizontal = TRUE)
partialPlot(fit$coefs, pred.data = df, x.var = "NH4M") 
partial(fit, pred.var = "NTOT", plot = TRUE, plot.engine = "ggplot2") #xy plot of curve
partial(fit, pred.var = "NH4M", plot = TRUE, plot.engine = "ggplot2") #xy plot of curve
partial(fit, pred.var = "NO3M", plot = TRUE, plot.engine = "ggplot2")
partial(fit, pred.var = "latitude", plot = TRUE, plot.engine = "ggplot2")
partial(fit, pred.var = "survey_year", plot = TRUE, plot.engine = "ggplot2")


fit %>% 
  partial(pred.var = "NTOT") %>%
  autoplot(smooth = TRUE) +
  theme_light() +
  ggtitle("ggplot2-based PDP")

# Compute partial dependence data for explanatory variables
pd <- partial(fit, pred.var = c("longitude", "latitude"), chull = TRUE)
pd <- partial(fit, pred.var = c("NH4M", "NO3M"), chull = TRUE)

# Default PDP
plotPartial(pd)

# Add contour lines and use a different color palette
plotPartial(pd, contour = TRUE)

# 3-D surface
plotPartial(pd, levelplot = FALSE, colorkey = TRUE)

#everything
data <- drop_na(temp2_scaled)
ind <- sample(2, nrow(data), replace = TRUE, prob=c(0.7, 0.3))
train <- data[ind == 1,]
test <- data[ind == 2,]

RFmodel <- randomForest(rich.y ~ n_nh4 + n_no3 + latitude + longitude + survey_year + 
                          #N.y + N.x + L.y + L.x + #T.x + T.y +
                          #div.x + div.y + rich.y
                          grp_tree_species + mean_temp +
                          mean_winter_temp + mean_precip + N.P +
                          mean_summer_temp + sum_canopy, mtry = 4, data = train)#train

fit <- RFmodel
randomForest::importance(fit)
varImpPlot(fit, main = "rich.x")
#div.x most important but highly correlated with rich.x, mean temps super correlated too.
#remove divs, temps except mean_temp. 80.4% explained
RFmodel <- randomForest(rich.x ~ n_nh4 + n_no3 + latitude + longitude + survey_year + 
                          N.x + L.x + T.x +
                          grp_tree_species + rich.y + mean_temp +
                          mean_precip + N.P +
                          sum_canopy, mtry = 4, data = train)
#take top 10 features 79.68%
RFmodel <- randomForest(rich.x ~ n_nh4 + n_no3 + latitude + longitude + 
                          N.x + T.x + L.x + sum_canopy +
                          mean_temp +
                          mean_precip, mtry = 4, data = train)

#take top 5 features 75.17%
RFmodel <- randomForest(rich.x ~ latitude + N.x + n_nh4 + mean_precip + mean_temp,
                        mtry = 2, data = train)

#top 6 Boruta 76.44%
RFmodel <- randomForest(rich.x ~ n_nh4 + latitude + longitude + mean_temp +
                          N.x + mean_precip, mtry = 4, data = train)



#variable selection for Random Forest (note:masks importance function)
library(Boruta)
b.train <- select(df, -c(dispersion, beta.SIM, beta.SNE, c(1:8)))

boruta.train <- Boruta(beta.SOR ~ . , data = b.train, doTrace = 2)
print(boruta.train)

plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)


# 
# # Split into Train and Validation sets
# # Training Set : Validation Set = 70 : 30 (random)
# set.seed(100)
# train <- sample(nrow(X), 0.7*nrow(X), replace = FALSE)
# TrainSet <- X[train,]
# ValidSet <- X[-train,]
# 
# TrainSet <- drop_na(TrainSet)
# ValidSet <- drop_na(ValidSet)
# 
# summary(TrainSet)
# summary(ValidSet)


# Create a Random Forest model with default parameters
model1 <- randomForest(N.x ~ n_nh4 + n_no3 + latitude + longitude + survey_year + rich.x +
                         rich.y + div.x + div.y + L.x + L.y +
                         sum_canopy, data = TrainSet, importance = TRUE)
model1
importance(model1)

# Predicting on train set
predTrain <- predict(model1, TrainSet)
# Checking classification accuracy
table(predTrain, TrainSet$N.x)  

# Predicting on Validation set
predValid <- predict(model1, ValidSet)
# Checking classification accuracy
mean(predValid == ValidSet$N.x)                    
table(predValid,ValidSet$N.x)

trControl = trainControl(method = 'cv', # Use cross-validation
                         number = 5) # Use 5 folds for cross-validation

# To check important variables
importance(model1)        
varImpPlot(model1)   

# Using For loop to identify the right mtry for model
a=c()
i=5
for (i in 3:8) {
  model3 <- randomForest(N.x ~ n_nh4 + n_no3 + latitude + longitude + survey_year + rich.x +
                           rich.y + div.x + div.y + L.x + L.y +
                           sum_canopy, data = TrainSet, ntree = 500, mtry = i, importance = TRUE)
  predValid <- predict(model3, ValidSet, type = "class")
  a[i-2] = mean(predValid == ValidSet$N.x)
}

a

plot(3:8,a)

