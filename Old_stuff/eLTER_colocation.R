elter_col <- elter_colocations_april_2024
elter_col$ICOS <- as.factor(elter_col$ICOS)
elter_col$SITES <- as.factor(elter_col$SITES)
elter_col$ICP_Forests <- as.factor(elter_col$`ICP Forests`)
elter_col$ICP_IM <- as.factor(elter_col$`ICP Integrated Monitoring`)
elter_col$ACTRIS <- as.factor(elter_col$ACTRIS)
elter_col$Lifewatch <- as.factor(elter_col$Lifewatch)
elter_col$Global_Cryosphere_Watch <- as.factor(elter_col$`Global Cryosphere Watch (GCW)`)

elter_col$`ICP Forests` <- NULL
elter_col$`ICP Integrated Monitoring` <- NULL
elter_col$`Global Cryosphere Watch (GCW)` <- NULL

table(elter_col$ICOS) #no 587 yes 38
table(elter_col$SITES) #no 615 yes 10
table(elter_col$ICP_Forests) #no 517 yes 108
table(elter_col$ICP_IM) #no 610 yes 15
table(elter_col$ACTRIS) #no 622 yes 5
table(elter_col$Lifewatch) #no 620 yes 5
table(elter_col$Global_Cryosphere_Watch) #no 613 yes 12

as.numeric(levels(elter_col$ICOS))[elter_col$ICOS] 

#correlations
elter_cor2 <- elter_col[,4:10] %>% mutate_if(is.character, ~ ifelse(. == "yes", 1, 0))
elter_cor2$sum <- rowSums(elter_cor2)
elter_coloc <- filter(elter_cor2, sum > 0)
elter_coloc$sum <- NULL
elter_cor_mat <-  cor(elter_coloc)
library(corrplot)
corrplot(elter_cor_mat, type = "upper")

