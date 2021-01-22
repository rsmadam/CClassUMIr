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
tuneGrid <- expand.grid("mtry" = c(25,40,50,60),
                       "splitrule" = "gini",
                       "min.node.size" = c(2,5,10) )
miniGrid <- expand.grid("mtry" = c(2,4,6,10),
            "splitrule" = "gini",
            "min.node.size" = c(1,2,3,5,10) )
miniParam <- expand.grid("mtry" = 10,
                         "splitrule" = "gini",
                         "min.node.size" = 2 )

rc_vst_BR$CMS <- droplevels(rc_vst_BR$CMS)
rc_vst_BR$Sample.ID <- row.names(rc_vst_BR)

set.seed(5678) 
out_Test <- createDataPartition(rc_vst_BR$CMS, p = .2, # this adapts to original fractions of classes
                                list = FALSE, 
                                times = 200)

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
  rc_vst_out <- as.data.frame(rc_vst_BR[out_Test[,i], ])
  #toss the out of training samples
  rc_vst_BR_tr <- rc_vst_BR[-out_Test[,i],]  
  
  ##identify best parameters 
  set.seed(5678) 
  model_RF_train <- caret::train(CMS~., 
                                 data=rc_vst_BR_tr[
                                   c(grep("hsa", colnames(rc_vst_BR), value = T)
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
  myFolds <- createFolds(rc_vst_BR_tr$CMS, k = 10, returnTrain = T) #representative samples
  testingContr <- trainControl(
    index = myFolds, #use updated myFolds 
    classProbs = T, #class probabilities for held-out samples during resample
    verboseIter = T,
    savePredictions = TRUE)
  set.seed(5678)
  
  ##use best parameters 
  model_RF_best <- caret::train(CMS~., 
                                data=rc_vst_BR_tr[, #c(N_imp_RF_df$variable[1:20],
                                                  c(grep("hsa", colnames(rc_vst_BR), value = T), 
                                                    "CMS")],
                                method="ranger", 
                                importance = 'impurity',
                                metric="Kappa", 
                                tuneGrid= bestParam, 
                                num.trees = 2000, #2000 tested best
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
                                   data=rc_vst_BR_tr[, c(N_imp_RF_df$variable[1:20],
                                                         "CMS")],
                                   method="ranger", 
                                   importance = 'impurity',
                                   metric="Kappa", 
                                   tuneGrid= miniParam, #fixed param
                                   num.trees = 2000, #2000 tested best
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
summary(unlist(store_acc_20)[1:50])
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
#           "analyses/tables/results_100x5cv_10x10cv_10cv.csv")


######## store CV importances ########
importances_df <- data.frame("1"=varImp(store_tests[[2]])$importance)
for(i in 81:100){
  importances_df[,as.character(i)] <- varImp(store_tests[[as.character(i)]])$importance[rownames(importances_df),"Overall"]
}
importances_df$Overall <- rowMeans(importances_df[,-1])
importances_df <- importances_df[order(-importances_df$Overall),]
 # write.csv(importances_df, row.names = T,
 #           "analyses/tables/val_100x5cv_10x10cv_importances.csv")

importances_df <- read.csv("analyses/tables/val_100x5cv_10x10cv_importances.csv",
                           row.names = 1)



############ big RF training on all COAD data ###############
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
# model_RF_train_all <- caret::train(CMS~., 
#                                data=rc_vst_BR[ ,c(grep("hsa", colnames(rc_vst_BR), value = T)
#                                    , "CMS")],
#                                method="ranger", 
#                                importance = 'impurity',
#                                metric="Kappa", 
#                                tuneGrid= tuneGrid, 
#                                num.trees = 2000, #2000 is best with 5 and 5 
#                                trControl= bigContr)
best_mtry <- 25 #model_RF_train_all$bestTune$mtry #
best_nnode <- 10 #model_RF_train_all$bestTune$min.node.size #
bestParam <- expand.grid("mtry" = best_mtry, # after identifying by grid search representative 10fCV, test different ntrees and classweight combis
                         "splitrule" = "gini",
                         "min.node.size" = best_nnode)
####### small RF 20 #######
set.seed(5678) 
# model_RF_20_train_all <- caret::train(CMS~., 
#                                       data=rc_vst_BR[
#                                         ,c(rownames(importances_df[1:20,])
#                                            , "CMS")],
#                                       method="ranger", 
#                                       importance = 'impurity',
#                                       metric="Kappa", 
#                                       tuneGrid= miniGrid,
#                                       num.trees = 2000, 
#                                       trControl= bigContr)
best_mtry_20 <- 2 #model_RF_20_train_all$bestTune$mtry #
best_nnode_20 <- 10 #model_RF_20_train_all$bestTune$min.node.size #
bestParam_20 <- expand.grid("mtry" = best_mtry_20, # after identifying by grid search representative 10fCV, test different ntrees and classweight combis
                         "splitrule" = "gini",
                         "min.node.size" = best_nnode_20)

set.seed(5678) 
myFolds <- createFolds(rc_vst_BR$CMS, k = 100, returnTrain = T) #representative samples
testingContr <- trainControl(
  index = myFolds, #use updated myFolds 
  classProbs = T, #class probabilities for held-out samples during resample
  verboseIter = T,
  savePredictions = TRUE)

###### training on all samples ######
set.seed(5678) 
model_RF_best_all <- caret::train(CMS~., 
                              data=rc_vst_BR[c(grep("hsa", colnames(rc_vst_BR), value = T), 
                                                  "CMS")],
                              method="ranger", 
                              importance = 'impurity',
                              metric="Kappa", 
                              tuneGrid= bestParam, 
                              num.trees = 2000, #2000 is best with 5 and 5 
                              trControl= testingContr)
set.seed(5678) 
model_RF_20_best_all <- caret::train(CMS~., 
                                  data=rc_vst_BR[
                                    c(rownames(importances_df[1:20,]),
                                      "CMS")],
                                  method="ranger", 
                                  importance = 'impurity',
                                  metric="Kappa", 
                                  tuneGrid= bestParam_20, 
                                  num.trees = 2000,  
                                  trControl= testingContr)


##############################################################################################################
######### CMS from the test data

######################################
### predictions on EGAS1127 (VU)-data?
VU_rc_vst$Sample.ID <- row.names(VU_rc_vst)
VU_valid <- VU_rc_vst[,]
pred_rc_vst_VU_RF <- predict(model_RF_best_all, newdata = VU_valid, type = "prob")

pred_rc_vst_VU_RF$CMS <- predict(model_RF_best_all, 
                                    newdata = VU_valid, type = "raw")
pred_rc_vst_VU_RF$CMS_20 <- predict(model_RF_20_best_all, 
                                 newdata = VU_valid, type = "raw")

rownames(pred_rc_vst_VU_RF) <- rownames(VU_valid)
pred_rc_vst_VU_RF$Sample.ID <-row.names(pred_rc_vst_VU_RF)
summary(pred_rc_vst_VU_RF$CMS)

### test if result differs when predicting only the primary samples
pred_VU_prim_CMS <- data.frame("CMS.prim"=predict(model_RF_best_all, 
                            newdata = VU_rc_vst[grep("primary",clinVU$sampleType),], 
                            type = "raw")) #predict only primary
pred_VU_prim_CMS$CMS.all <- pred_rc_vst_VU_RF$CMS[match(clinVU[grep("primary",clinVU$sampleType),
                                                              "sampleID"], #match to all
                                                     pred_rc_vst_VU_RF$Sample.ID)]
summary(pred_VU_prim_CMS$CMS.prim==pred_VU_prim_CMS$CMS.all) # no they don't differ


### add VU predictions to clinical annotation table
clinVU$CMS <- 
  pred_rc_vst_VU_RF[match(clinVU$sampleID,
   rownames(pred_rc_vst_VU_RF)),"CMS"] 
clinVU$CMS_20 <- 
  pred_rc_vst_VU_RF[match(clinVU$sampleID,
                          rownames(pred_rc_vst_VU_RF)),"CMS_20"] 

# write.csv2(clinVU, row.names = F, 
#            "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/tables/VUclin-predictedCMS_RF-all_RF-20.csv")
#clinVU <- read.csv2( "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/tables/VUclin-predictedCMS_RF-all_RF-20.csv" )


#######################
##### READ data 
library(RColorBrewer)

miR_READ_vst_BR <- read.csv(file=
                              "/Data/TCGA-miR_READ_vst_BR_CMS-labels.csv",
                            row.names = 1)
miR_READ_vst
selREAD <- 1:length(miR_READ_vst_BR$CMS.cl.rf)#which(  
  #grepl("CMS", miR_READ_vst_BR$CMS.lv ) &
    #miR_READ_vst_BR$CMS.lv == miR_READ_vst_BR$CMS.cl.rf)
pred_read_RF <- predict(model_RF_best_all,
                        newdata = miR_READ_vst_BR[selREAD,],
                                                  #which(miR_READ_vst_BR$CMS.lv == miR_READ_vst_BR$CMS.cl.rf)]),],
                        type = "prob")
summary(pred_read_RF)
###better method to censor by post prob
diffSecond <- function(x, output) {
  ordered = sort(x, decreasing=T)
  censor = abs(ordered[1] - ordered[2])
  censor
} #juts get the absolute difference between first and second class
theSecond <- function(x, output) {
  ordered = order(x, decreasing=T)
  second = ordered[2]
  second
} #juts get the absolute difference between first and second class

pred_read_RF$d2ndProb <- apply(pred_read_RF, 1, diffSecond)
pred_read_RF$d2ndClass <- apply(pred_read_RF[,1:4], 1, theSecond)

pred_read_RF_20 <- predict(model_RF_20_best_all, #get distances also for RF_20
                        newdata = miR_READ_vst_BR[selREAD,],
                        type = "prob")
pred_read_RF$d2ndProb_20 <- apply(pred_read_RF_20, 1, diffSecond)
pred_read_RF$d2ndClass_20 <- apply(pred_read_RF_20, 1, theSecond)

pred_read_RF$CMS <- predict(model_RF_best_all, 
                               newdata = miR_READ_vst_BR[selREAD
                                                         ,grep("hsa", colnames(miR_READ_vst_BR))],
                               type = "raw")
pred_read_RF$CMS_20 <- predict(model_RF_20_best_all, 
                            newdata = miR_READ_vst_BR[selREAD
                                                   ,grep("hsa", colnames(miR_READ_vst_BR))],
                            type = "raw")
pred_read_RF$CMS_tumiR <- predict(model_RF_best_100_TumiRs, 
                               newdata = miR_READ_vst_BR[selREAD
                                                         ,grep("hsa", colnames(miR_READ_vst_BR))],
                               type = "raw")
rownames(pred_read_RF) <- rownames(miR_READ_vst_BR[selREAD,])
pred_read_RF$Sample.ID <-row.names(pred_read_RF)
summary(pred_read_RF$CMS_20)
confusionMatrix(pred_read_RF$CMS, factor(miR_READ_vst_BR$CMS.lv[selREAD]))
confusionMatrix(pred_read_RF$CMS_20, factor(miR_READ_vst_BR$CMS.lv[selREAD]))

###censoring low confidence helps:
selConfid <- which(pred_read_RF$d2ndProb > quantile(pred_read_RF$d2ndProb, prob=.25))
confusionMatrix(pred_read_RF$CMS[selConfid], factor(miR_READ_vst_BR$CMS.lv[selConfid]))
selConfid <- which(pred_read_RF$d2ndProb_20 > quantile(pred_read_RF$d2ndProb_20, prob=.25))
confusionMatrix(pred_read_RF$CMS_20[selConfid], factor(miR_READ_vst_BR$CMS.lv[selConfid]))


##########################################
### colors CMS
#Guinney:
paletteCMS=c( "#E79E1B", 
              "#0071B1", 
             "#C45597",#Guinney:"#C977A4", 
              "#009C74")

paletteCMSn=c("CMS1" ="#E79E1B", 
              "CMS2"= "#0071B1", 
              "CMS3"= "#C45597",#Guinney:"#C977A4", 
              "CMS4"="#009C74",
              "NOLBL"="#d3d3d3")
