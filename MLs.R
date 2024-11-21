
#Set the folder for data analyses
setwd("/Users/____/Desktop/Data")

#For ELA (Fig.3) of ML (Fig.S10)
#Store the raw data for ELA within the folder "Data" 
#ELA_raw.csv (stored in the "2.Binarized data" of the sheet name "Figs.3 and S11" of DataFCPCfinal241115.xlsx)
#ELA_grouplist.csv (stored in the "3.Group list for ELA" of the sheet name "Figs.3 and S11" of DataFCPCfinal241115.xlsx)

#Make th folder for ELA
dir.create("./Result_ELA", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Analyze based on the following protocol
#(https://github.com/kecosz/ela; https://github.com/kecosz/rELA) and the DOIs (10.5281/zenodo.5492161 and 10.5281/zendo.7979838). 

#Store the calculated data within the folder "Result_ELA" 

#For AA (Fig.4a) 
#memory clear
rm(list=ls(all=TRUE))
invisible(replicate(20, gc()))
library(dplyr)

#Set the folder for data analyses
setwd("/Users/____/Desktop/Data")

#Set the R library
if (!require("arules")) install.packages("arules")
if (!require("arulesViz")) install.packages("arulesViz")
require(arules) 
require(arulesViz)

dir.create("./Result_AA_RF_XG", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#Read the binarized file for AA (stored "Bibarized data for association analysis" in the sheet name "Fig.4a" of DataFCPCfinal241115.xlsx)
readfile="AA_binarized_raw_data.csv"

x=read.csv(readfile, header=T, colClasses="factor")
rulesAp1 <- apriori(x, parameter=list(support=0.25,confidence=0.5,maxlen=2))
rule_Ap1<-rulesAp1[quality(rulesAp1)$lift>1.2]

write(rule_Ap1,"./Result_AA_RF_XG/AA_calculated_data.csv",sep=",")

#For RF and XGBoost (Fig.4a) 
#Read the raw data file for AA-selected components
#For XGBoost
library(xgboost)
FCPC.train <- read.csv("RF_XG_raw_data.csv") # tibble::glimpse() 
FCPC.test <- read.csv("RF_XG_raw_data.csv") #RF_XG_raw_data.csv
dim(FCPC.train) 
train.x <- FCPC.train[, 2:24] #dim(FCPC.train)[1] 36 24
x <- rbind(train.x,FCPC.test[,-1]) 

y_0 <- c(FCPC.train$Species)
y_1 <- as.factor(y_0)
y <- as.integer(y_1)-1

x <- as.matrix(x)
trind <- 1:length(y) 
teind <- (nrow(train.x)+1):nrow(x) 
set.seed(131) #fixed seed

param <- list("objective" = "multi:softprob", 
              "eval_metric" = "mlogloss", 
              "num_class" = 8 # class (Fish_Con; Fish_Comp; Chicken_Con; Chicken_Comp; Pig_Con; Pig_Comp_ThB; Cattle_Before; Cattle_After)
)

k<-round(1+log2(nrow(train.x)))
cv.nround <- 100 #search
bst.cv <- xgb.cv(param=param, data = x[trind,], label = y,  nfold = k, nrounds=cv.nround)

set.seed(131)
nround <- 27

bst <- xgboost(param=param, data = x[trind,], label = y, nrounds=nround)
pred <- predict(bst,x[teind,]) 
pred <- matrix(pred,8,length(pred)/8)　# class (Fish_Con; Fish_Comp; Chicken_Con; Chicken_Comp; Pig_Con; Pig_Comp_ThB; Cattle_Before; Cattle_After)　　
pred <- t(pred)
colnames(pred)<-c("FE_Control","FE_Test","CW_Control","CW_Test","PK_Control","PK_Test","KC_Control","KC_Test") 

head(pred,8) #8

x_1_f <- FCPC.test[,1]
for(i in 1:length(pred)){
  if(pred[i]==0) {pred[i]="Control"}　#Species Control Test
  else if(pred[i]==1) {pred[i]="Test"}
}

table(x_1_f,pred)

sink('./Result_AA_RF_XG/XGboost_pre_x_1_f.txt', append = TRUE)
print (table(x_1_f,pred))
sink()

write.csv(table(x_1_f,pred),"./Result_AA_RF_XG/XGboost_pre_x_1_f.csv")

imp<-xgb.importance(names(y_1),model=bst)
print(imp)
xgb.plot.importance(imp) 

pdf ("./Result_AA_RF_XG/XGBoostgraph.pdf") 
xgb.plot.importance(imp) 
dev.off()

write.csv(print(imp),"./Result_AA_RF_XG/XGBoostgraph_raw.csv")


#For RF
library(randomForest)
set.seed(131)
train.x<- FCPC.train[,2:24] #dim(FCPC.train) [1]  8 74
train.y<-as.factor(FCPC.train[,1])
model.rf<-tuneRF(train.x,train.y,doBest=T)
pred<-predict(model.rf,FCPC.test[,2:24]) #dim(FCPC.train) [1]  8 74
table(FCPC.test[,1],pred)
print(model.rf$importance /sum(model.rf$importance))

write.csv(print(model.rf$importance /sum(model.rf$importance)),"randomForest_raw.csv")

rf_pred <- table(FCPC.test[,1],pred)
write.csv(rf_pred,"./Result_AA_RF_XG/randomForest_pred.csv")

FCPC.test_n <- cbind(y_1, train.x)

set.seed(22)
model = randomForest(y_1 ~ ., data = FCPC.test_n, importance = TRUE, proximity = TRUE)
print(model)
print(varImpPlot(model))

write.csv(print(varImpPlot(model)),"./Result_AA_RF_XG/randomForest_pred_importance_Gini.csv")

par(mar=c(100, 20, 30, 40)) #par(oma = c(3, 3, 3, 2))
rpp2 <- varImpPlot(model)
varFileName <- paste("./Result_AA_RF_XG/randomforest_tree_var.png",sep="") #フォルダの位置確認
par(mar=c(100, 20, 30, 600)) 
png(file=varFileName, res=125, w=750, h=750)
rpp2 <- varImpPlot(model)
dev.off()

#Make the file (ML_mix.xlsx) containing the values of components selected by AA,RF, and XGBoost
#Use python to illustrate the Bubble chart(Fig.4a) 
#Bubblechart.py (load ML_mix.xlsx)