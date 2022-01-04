# Libraries needed

library(MASS)
library(tidyverse)
library(randomForest)
library(caret)
library(nnet)
library(pROC)
library(vsn)
library(RColorBrewer)
library(wesanderson)
library(ggsci)
library(factoextra)
library(circlize)
library(ComplexHeatmap)
library(ggpubr)
library(pheatmap)
library(viridis)

plotDir <- "/Users/ronjaadam/projects/miRNA_mCRC/CMS-miRaCl/analyses/plots/"



##################################################################################################################################      
####################### RF training with Cross validation #######################

########################################################################

#### testing accuracy of RF ranger in caret ####

##create train vs test splits beforehand
# Create custom indices: myFolds

trainingContr <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = F, #no class probabilities for held-out samples during resample
  verboseIter = F,
  sampling = "down", 
  search="grid",
  savePredictions = F
)
tuneGrid <- expand.grid("mtry" = c(2,5,10,25),
                       "splitrule" = "gini",
                       "min.node.size" = c(2,5,10) )
miniGrid <- expand.grid("mtry" = c(2,4,6,10),
            "splitrule" = "gini",
            "min.node.size" = c(2,5,10) )
miniParam <- expand.grid("mtry" = 10,
                         "splitrule" = "gini",
                         "min.node.size" = 2 )

miR.arr.scaled$CMS <- factor(miR.arr.scaled$CMS)
miR.arr.scaled$Sample.ID <- row.names(miR.arr.scaled)
outlierM <- rownames(miR.arr.scaled[is_outlier(rowSds( as.matrix(miR.arr.scaled[, 
                                                                         grep("hsa",
                                                                              colnames(miR.arr.scaled))]))),])#outliers with low mean AND high SD (across all samples)
miR.arr.scaled <- miR.arr.scaled[setdiff(rownames(miR.arr.scaled), outlierM),]

## identify highly correlating (redundant) features
nzv <- nearZeroVar(miR.arr.scaled, saveMetrics=TRUE)
which(nzv$nzv)
descrCor <-  cor(miR.arr.scaled[, setdiff(
  grep("hsa", colnames(miR.arr.scaled)),
                                     which(nzv$nzv))
]) #exclude zero var
highlyCorDescr <- c(colnames(descrCor[,findCorrelation(descrCor, 
                                                              cutoff = .75, exact=T)]), #107
                           colnames(miR.arr.scaled[,nzv$nzv])) #there's no overlap between them
colnames(miR.arr.scaled[,nzv$nzv]) #35

## exclude higly correlating(redundant) features
miR.arr.scaled <- miR.arr.scaled[,setdiff(colnames(miR.arr.scaled),
                                               highlyCorDescr)]
dim(miR.arr.scaled)

## plot correlating features
library(corrplot)
pdf("/Users/ronjaadam/projects/miRNA_mCRC/CMS-miRaCl/analyses/plots/miRaCl20A_correlations_training_20.pdf",
    useDingbats = F)
corrplot(descrCor[match(rownames(importances_df_arr_nzv[1:20,]), rownames(descrCor)),
                  c( ##which features of miRaCL could it correspond to? 
                     match(rownames(importances_df[1:20,]),
                           colnames(descrCor))) ], 
         col=rev(brewer.pal(n=8, name="RdYlBu")), 
         na.label = "o", na.label.col = "#EFEFEF",
           addgrid.col = NA, tl.cex = 1.2,
         tl.col = "#999999")
dev.off()

## these features were 
intersect(rownames(importances_df[1:20,]), highlyCorDescr)#,
      #colnames(descrCor))
## "hsa.mir.218" "hsa.mir.143" "hsa.mir.99a" "hsa.mir.141"


set.seed(5678) 
out_Test <- createDataPartition(miR.arr.scaled$CMS, p = .2, # this adapts to original fractions of classes
                                list = FALSE, 
                                times = 100)

# Make empty variables
store_tests_20 <- list()
store_mtry_20 <- list()
store_nnode_20 <- list()
store_acc_20 <- list()
store_kap_20 <- list()
store_pval_20 <- list()
store_tests <- list()
store_mtry <- list()
store_nnode <- list()
store_acc <- list()
store_kap <- list()
store_pval <- list()


####################### Make the RF model ####################### 
### with  caret, add external CV to get average accuracy and average importances ###

set.seed(5678) 
for(i in c(1:100)){
  
  ##make the out of training data
  rc_vst_out <- as.data.frame(miR.arr.scaled[out_Test[,i], ])
  #toss the out of training samples
  miR.arr.scaled_tr <- miR.arr.scaled[-out_Test[,i],]  
  
  ##identify best parameters 
  set.seed(5678) 
  model_RF_train <- caret::train(CMS~., 
                                 data=miR.arr.scaled_tr[
                                   c(grep("hsa", colnames(miR.arr.scaled), value = T)
                                     , "CMS")],
                                 method="ranger", 
                                 importance = 'impurity',
                                 metric="Kappa", 
                                 tuneGrid= tuneGrid,#miniGrid,# 
                                 num.trees = 1000, #2000 tested best
                                 trControl= trainingContr)
  best_mtry <- model_RF_train$bestTune$mtry
  best_nnode <- model_RF_train$bestTune$min.node.size
  bestParam <- expand.grid("mtry" = best_mtry, # after identifying by grid search representative 10fCV, test different ntrees and classweight combis
                           "splitrule" = "gini",
                           "min.node.size" = best_nnode)
  set.seed(5678) 
  myFolds <- createFolds(miR.arr.scaled_tr$CMS, k = 10, returnTrain = T) #representative samples
  testingContr <- trainControl(
    index = myFolds, #use updated myFolds 
    classProbs = T, #class probabilities for held-out samples during resample
    verboseIter = T,
    savePredictions = TRUE)
  set.seed(5678)
  
  ##use best parameters 
  model_RF_best <- caret::train(CMS~., 
                                data=miR.arr.scaled_tr[,
                                                  c(grep("hsa", colnames(miR.arr.scaled), value = T), 
                                                    "CMS")],
                                method="ranger", 
                                importance = 'impurity',
                                metric="Kappa", 
                                tuneGrid= bestParam, 
                                num.trees = 2000,
                                trControl= testingContr)
  #### get stats from best RF
  ## create the confusion matrix for test data
  pred_iter_valid <- predict(model_RF_best, 
                             newdata = rc_vst_out[,grep("hsa",colnames(rc_vst_out))],
                             type = "raw")
  cmat_iter_RF <- confusionMatrix(pred_iter_valid, na.omit(rc_vst_out$CMS)) 
  ## save the data
  key <- paste( toString(i) )
  store_tests[[key]] <- model_RF_best
  store_mtry[[key]] <- best_mtry
  store_nnode[[key]] <- best_nnode
  store_acc[[key]] <- cmat_iter_RF$overall[["Accuracy"]]
  store_kap[[key]] <- cmat_iter_RF$overall[["Kappa"]]
  store_pval[[key]] <- cmat_iter_RF$overall[["AccuracyPValue"]]
  
  #### reduce to only 20 features
  N_imp_RF <- varImp(model_RF_best)$importance
  N_imp_RF_df <- tibble(variable = rownames(N_imp_RF), 
                        importance = N_imp_RF$Overall ) %>% 
    dplyr::arrange(-importance) 
  model_RF_best_20 <- caret::train(CMS~., 
                                   data=miR.arr.scaled_tr[, c(N_imp_RF_df$variable[1:20],
                                                         "CMS")],
                                   method="ranger", 
                                   importance = 'impurity',
                                   metric="Kappa", 
                                   tuneGrid= miniParam, #fixed param
                                   num.trees = 2000,
                                   trControl= testingContr)
  ## create the confusion matrix for test data
  pred_iter_valid <- predict(model_RF_best_20, 
                             newdata = rc_vst_out[,grep("hsa",colnames(rc_vst_out))],
                             type = "raw")
  cmat_iter_RF <- confusionMatrix(pred_iter_valid, na.omit(rc_vst_out$CMS)) 
  
  ## save the data
  store_tests_20[[key]] <- model_RF_best_20
  store_mtry_20[[key]] <- best_mtry
  store_nnode_20[[key]] <- best_nnode
  store_acc_20[[key]] <- cmat_iter_RF$overall[["Accuracy"]]
  store_kap_20[[key]] <- cmat_iter_RF$overall[["Kappa"]]
  store_pval_20[[key]] <- cmat_iter_RF$overall[["AccuracyPValue"]]
  print(paste("finishing loop", i))
  
}


####### extract average results #######
apply(pred_cptec_RF, 2, function(x) summary(as.factor(x)))
results_final <- resamples(store_tests_20)
summary(results_final)
summary(unlist(store_acc_20)[1:75])
hist(unlist(store_acc_20), breaks=10)
plot(unlist(store_acc_20))
plot(data.frame("Acc"=unlist(store_acc), 
                "Mtry"=unlist(store_mtry), 
                "Node"=unlist(store_nnode))) 
plot(data.frame("Acc"=unlist(store_acc_20), 
                "Mtry"=unlist(store_mtry_20)))
summary(unlist(store_pval))

# write.csv(data.frame("Acc"=unlist(store_acc),
#                      "pval"=unlist(store_pval),
#                      "Kappa"=unlist(store_kap),
#                      "Mtry"=unlist(store_mtry),
#                      "Node"=unlist(store_nnode)),
#           row.names = F,
#           "analyses/tables/results_arrayFFPE_nzv_1-75x5cv_10x10cv_10cv.csv")


######## store CV importances ########
importances_df_arr_nzv <- data.frame("1"=varImp(store_tests[[2]])$importance)
for(i in 1:75){
  importances_df_arr_nzv[,as.character(i)] <- varImp(store_tests[[as.character(i)]])$importance[
    rownames(importances_df_arr_nzv),"Overall"]
}
importances_df_arr_nzv$Overall <- rowMeans(importances_df_arr_nzv[,-1])
importances_df_arr_nzv <- importances_df_arr_nzv[order(-importances_df_arr_nzv$Overall),]
 # write.csv(importances_df_arr_nzv, row.names = T,
 #           "analyses/tables/val_1-75x5cv_10x10cv_importances_arrayFFPE_nzv.csv")

importances_df_arr_nzv <- read.csv("analyses/tables/VAL_1-75x5cv_10x10cv_importances_arrayFFPE_nzv.csv",
                           row.names = 1, sep="\t")
plot(importances_df_arr_nzv$Overall[1:50])
rownames(importances_df_arr_nzv[1:50,])


############ big RF training on FFPE data ###############
#### create parameters for all samples ####
set.seed(5678) 
bigContr <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 100,
  classProbs = F, #no class probabilities for held-out samples during resample
  verboseIter = F,
  sampling = "down", 
  search="grid",
  savePredictions = F
)
set.seed(5678) 
best_mtry <- 25 #model_RF_train_all$bestTune$mtry #
best_nnode <- 10 #model_RF_train_all$bestTune$min.node.size #
bestParam <- expand.grid("mtry" = best_mtry, # after identifying by grid search representative 10fCV, test different ntrees and classweight combis
                         "splitrule" = "gini",
                         "min.node.size" = best_nnode)
####### small RF 20 #######
set.seed(5678) 
best_mtry_20 <- 5 #model_RF_20_train_all$bestTune$mtry #
best_nnode_20 <- 10 #model_RF_20_train_all$bestTune$min.node.size #
bestParam_20 <- expand.grid("mtry" = best_mtry_20, # after identifying by grid search representative 10fCV, test different ntrees and classweight combis
                         "splitrule" = "gini",
                         "min.node.size" = best_nnode_20)

set.seed(5678) 
myFolds <- createFolds(miR.arr.scaled$CMS, k = 100, returnTrain = T) #representative samples
testingContr <- trainControl(
  index = myFolds, #use updated myFolds 
  classProbs = T, #class probabilities for held-out samples during resample
  verboseIter = T,
  savePredictions = TRUE)

###### training on all samples ######
set.seed(5678) 
model_RF_best_ffpe <- caret::train(CMS~., 
                              data=miR.arr.scaled[c(grep("hsa", colnames(miR.arr.scaled), value = T), 
                                                  "CMS")],
                              method="ranger", 
                              importance = 'impurity',
                              metric="Kappa", 
                              tuneGrid= bestParam, 
                              num.trees = 2000, #2000 is best with 5 and 5 
                              trControl= testingContr)
set.seed(5678) 
model_RF_20_best_ffpe <- caret::train(CMS~., 
                                  data=miR.arr.scaled[
                                    c(rownames(importances_df_arr_nzv[1:20,]),
                                      "CMS")],
                                  method="ranger", 
                                  importance = 'impurity',
                                  metric="Kappa", 
                                  tuneGrid= bestParam_20, 
                                  num.trees = 2000,  
                                  trControl= testingContr)

