# Libraries needed
#BiocManager::install("vsn")
#BiocManager::install("ConsensusClusterPlus")

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


#install.packages('e1071', dependencies=TRUE)
set.seed(67)
plotDir <- "/Users/ronjaadam/projects/miRNA_mCRC/CMS-miRaCl/analyses/plots/"

# ##################################################################################################################################      


tuMirs <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/TuSpecUp_Neerincx-Oncgs15_S3.csv",
                    as.is=T, header = F)
tuMirs <- sub("miR", "mir", gsub("-", ".", sub("-.p$","",tuMirs$V1)))
tuMirsDn <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/TuSpecDn_Neerincx-Oncgs15_S4.csv",
                    as.is=T, header = T)  
tuMirsDn <- sub("miR", "mir", gsub("-", ".", sub("-.p$","",tuMirsDn$miRNA)))
tuMirsDn <- sub("\\.[1,2]$", "", tuMirsDn)
unique(tuMirsDn[which(tuMirsDn %in% colnames(miR_COAD_vst))])
unique(tuMirsDn[which(tuMirsDn %in% rownames(importances_df[1:50,]))])#potentially useful 
tuMirs <- unique( c(tuMirs, tuMirsDn) )

unique(tuMirs[which(tuMirs %in% colnames(miR_COAD_vst))])#only new miRs and mir.92a.1 but mir.92a is also there


#rc_vst <- miR_COAD_vst[grep("CMS", miR_COAD_vst$CMS), 
#                       unique(tuMirs[which(tuMirs %in% colnames(miR_COAD_vst) &
#                     tuMirs %in% colnames(VU_rc_vst))])]## CAVE this still contains non-expressed genes, which should drop out in correlated/nonvariance test



##################################################################################################################################      
####################### RF training with Cross validation ####################### 

##### preprocessing: visualize ######


### PCA plot ###
#install.packages("factoextra")
prin_comp <- prcomp(rc_vst_BR[,-c(grep("S", colnames(rc_vst_BR)))], scale. = T)
pr_var <- round( prin_comp$sdev^2, 5 )
prop_varex <- round( pr_var/sum(pr_var), 5)
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
var.coord <- t(apply(prin_comp$rotation, 1, var_coord_func, prin_comp$sdev)) 
head(var.coord[order(var.coord[,"PC1"]), "PC1"]) #the top features
tail(var.coord[order(var.coord[,"PC1"]), "PC1"])

trans <- preProcess(rc_vst_BR[ #-grep("AA", rownames(rc_vst_BR))
                        ,grep("hsa",colnames(rc_vst_BR))], 
                    method=c(#"BoxCox",
                             "center","scale", 
                             "pca"),
                    thresh = list(thresh = 0.50))
# head(trans$rotation[order(trans$rotation[,"PC1"]), "PC1"]) #the top neg features
# tail(trans$rotation[order(trans$rotation[,"PC1"]), "PC1"]) #the top features
# featPC134 <- unique( names( c( head(trans$rotation[order( 
#   abs( trans$rotation[,"PC3"] ) ), "PC3"], 20), #the top abs features
#   head(trans$rotation[order( 
#     abs( trans$rotation[,"PC4"] ) ), "PC4"], 20), 
#   head(trans$rotation[order( 
#     abs( trans$rotation[,"PC1"] ) ), "PC1"], 20) ) ) )

PC <- predict(trans, rc_vst_BR[#-grep("AA", rownames(rc_vst_BR))
                         ,grep("hsa",colnames(rc_vst_BR))])

# PC1 shows clear separartion of 2 groups independent of CMS, 
# is other clinical data associated with PC1? 
# nzv <- nearZeroVar(clinical, saveMetrics= TRUE)
# summary(nzv$nzv)
# clinical <- clinical[,-which(nzv$nzv)] #toss some uselessclin variables
# 
# clinVar <- clinical$submitter_id[ match( sub("-01A.*","", 
#                                                row.names( PC[ order(PC$PC1), ] ) ),
#                                            clinical$submitter_id ) ]
# clinVar <- factor(as.factor(clinVar), labels=c("right", "right", 
#                                                "Colon,NOS", "left",
#                                                "right", "left",
#                                                "left", "left", 
#                                                "transv"))
# clinVar <- sub("[A,B,C]", "", sub("Stage ", "", clinVar) )
# clinVar <- sub("-.*", "", sub("TCGA-", "", clinVar) )
##some colors for the different variables: 
#scale_color_aaas(alpha=0.8)#lancet() #jco # npg # aaas
#scale_color_manual(values = wes_palette("Darjeeling2", 5, type = c("discrete")))
#scale_color_manual(values = wes_palette("Darjeeling1", 15, 
#                                        type = c( "continuous")) ) ## used this for the TSS
#scale_colour_viridis_d()
#scale_color_brewer(palette = "Dark2", direction = -1)
## results of the plotting dfferent clinical variables on top of the PC1/2:
#-> the TSS, tissue source site, i.e. "AA" from sample bank "indivumed" accounted batch effect
#-> put this variable into batch effect remover

## plot PCA
library(ggsci)
pdf(paste0(plotDir, "miRNA_vst_limBR-TSS2_-outl-lowVarMir-PCA.pdf"),
    onefile = T)
fviz_eig(prin_comp, type="lines" ) #plot variances
ggplot(PC,aes(x=PC1,y=PC2, 
              colour= rc_vst_BR$CMS)) +
  geom_point(na.rm = F) +
  theme_minimal() +
  scale_color_manual(values=paletteCMS)
ggplot(PC,aes(x=PC3,y=PC4,
              colour= rc_vst_BR$CMS)) +
  geom_point(na.rm = F) +
  theme_minimal() +
  scale_color_manual(values=paletteCMS)
ggplot(PC,aes(x=PC5,y=PC6,
              colour= rc_vst_BR$CMS)) +
  geom_point(na.rm = F) +
  theme_minimal() +
  scale_color_manual(values=paletteCMS)
dev.off()

#### plot cluster dendrogram ####
hc <- hclust(dist(rc_vst_BR[,grep("hsa",colnames(rc_vst_BR))]), "ward.D") 
colrc_vst.CMS <- factor(rc_vst_BR$CMS, labels = c("darkgoldenrod1", "blue4",
                                                  "deeppink3",
                                                  "mediumseagreen"))#make color code
colLab <- function(n) {#function to get color labels
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- as.character(colrc_vst.CMS[which(rownames(rc_vst_BR) == a$label)])
    attr(n, "nodePar") <- c(a$nodePar,pch="o",col = labCol) }
  return(n) }
hcd <- as.dendrogram(hc)
clusDendro <- dendrapply(hcd, colLab)# using dendrapply to get color-annotation for each node
plot(clusDendro)

#### plot kmeans class labels on PCA #####
kM <- kmeans(rc_vst_BR[, -grep("S", colnames(rc_vst_BR))], 
             4, nstart = 10)
kM$cluster <- as.factor(kM $cluster)
summary(kM$cluster)
ggplot(rc_vst_BR, aes(PC$PC1, 
                   PC$PC2, 
                   shape = kM$cluster,
                   color=rc_vst_BR$CMS))+
         theme_minimal()+
         scale_color_manual(values=paletteCMS)+
         geom_point() #ok
ggsave(paste0(plotDir, "miRNA_vst_limBR_-outl-lowVarMir-PCA-kMeans.pdf"))

#### tSNE plot ####
set.seed(9)  
library(Rtsne)
library(ggplot2)
library(ggsci)
tsne_model_1 <- Rtsne(as.matrix(rc_vst_BR[,grep("hsa",colnames(rc_vst_BR))]), 
                      check_duplicates=FALSE, 
                     pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_1 <- as.data.frame(tsne_model_1$Y) 
ggplot(d_tsne_1, aes(x=V1, y=V2, colour= rc_vst_BR$CMS))+#rc_vst_BR$CMS)) + #colour=kM$cluster))+# 
  geom_point() +
  xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme_minimal()+
  # scale_color_npg()+
  # scale_color_manual(values = wes_palette("Darjeeling1", 22,
  #                                         type = c( "continuous")) ) ## used this for the TSS
  scale_colour_manual(values = paletteCMS)
ggsave(paste0(plotDir, "miRNA_vst_limBR-TSS2_-outl-lowvarMir_tSNE_CMS.pdf"))

## from the kmeans it looks like there's a good concordance between 
# CMS1 and kCl2, CMS2 and kCl1, CMS4 and kCl3, 
# whereas CMS3 and kCl4 is more outliers and intermediate samples
# leave out CMS3 samples for everything? or combined class CMS2+CMS3? No
# only concordant samples in the training helps with PCA 5&6 now separating CMS3



```


```{r}
########################################################################

#### testing accuracy of RF ranegr in caret ####


##create train vs test splits beforehand
# Create custom indices: myFolds

trainingContr <- trainControl(
  #index = myFolds, #either use fixed myfolds or repeated CV
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
#TODO use metric "ROC" for model comparison,
#TODO: keep nzv but use preprocess pca instead
#TODO: caret::train on whole dataset because it uses the trainControl function as training recipe and to compute stats on cv out of trainin
#TODO: try bootstrap aggregation instead of cv
#TODO: compare different models in model_list using e.g. summary() or dotplot(resamples(model_list), metric="ROC)


###### Classifier with RandomForest or caret ######
rc_vst_BR$CMS <- droplevels(rc_vst_BR$CMS)
rc_vst_BR$Sample.ID <- row.names(rc_vst_BR)

# Make an empty variable
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
set.seed(5678) 
out_Test <- createDataPartition(rc_vst_BR$CMS, p = .2, # this adapts to original fractions of classes
                                list = FALSE, 
                                times = 200)
set.seed(5678) 
for(i in c(61:80)){
  ####################### Make the RF model ####################### 
  ##with  caret

  ##make the out of training data
  rc_vst_out <- as.data.frame(rc_vst_BR[out_Test[,i], ])
  #toss the out of training samples
  rc_vst_BR_tr <- rc_vst_BR[-out_Test[,i],]  
  
  ##identify best parameters 
  set.seed(5678) 
  model_RF_train <- caret::train(CMS~., 
                           data=rc_vst_BR_tr[
                                          #c(N_imp_RF_df$variable[1:20]
                             c(grep("hsa", colnames(rc_vst_BR), value = T)
                                                 , "CMS")],
                           method="ranger", 
                           importance = 'impurity',
                           metric="Kappa", 
                           tuneGrid= tuneGrid,#miniGrid,# 
                           num.trees = 1000, #2000 is best with 5 and 5 
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
                                num.trees = 2000, #2000 is best with 5 and 5 
                                trControl= testingContr)
  ## get stats from best RF
  # create the confusion matrix for test data
  pred_iter_valid <- predict(model_RF_best, 
                             newdata = rc_vst_out[,grep("hsa",colnames(rc_vst_out))],
                             type = "raw")
  cmat_iter_RF <- confusionMatrix(pred_iter_valid, na.omit(rc_vst_out$CMS)) 
  
  ## save all data
  key <- paste( toString(i) )
  store_tests[[key]] <- model_RF_best
  store_mtry[[key]] <- best_mtry
  store_nnode[[key]] <- best_nnode
  store_acc[[key]] <- cmat_iter_RF$overall[["Accuracy"]]
  store_kap[[key]] <- cmat_iter_RF$overall[["Kappa"]]
  store_pval[[key]] <- cmat_iter_RF$overall[["AccuracyPValue"]]
  
  ## reduce to only 20 features
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
                                 num.trees = 2000, #2000 is best with 5 and 5 
                                 trControl= testingContr)
  ## create the confusion matrix for test data
  pred_iter_valid <- predict(model_RF_best_20, 
                             newdata = rc_vst_out[,grep("hsa",colnames(rc_vst_out))],
                             type = "raw")
  cmat_iter_RF <- confusionMatrix(pred_iter_valid, na.omit(rc_vst_out$CMS)) 
 
 store_tests_20[[key]] <- model_RF_best_20
 store_mtry_20[[key]] <- best_mtry
 store_nnode_20[[key]] <- best_nnode
 store_acc_20[[key]] <- cmat_iter_RF$overall[["Accuracy"]]
 store_kap_20[[key]] <- cmat_iter_RF$overall[["Kappa"]]
 store_pval_20[[key]] <- cmat_iter_RF$overall[["AccuracyPValue"]]
 print(paste("finishing loop", i))
 
}

apply(pred_cptec_RF, 2, function(x) summary(as.factor(x)))
#
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
#           "analyses/tables/results_100x5cv_51-80_10x10cv_10cv.csv")


### make boxplot for Accuracy and Kappa
results_table <- as.data.frame(matrix(nrow = 100,ncol=2))
results_table$V1 <- factor(rep(c("Accuracy", "Kappa") , each=50))
results_table$V2 <- c(unlist(store_acc_20), unlist(store_kap_20))
colnames(results_table) <- c("variable", "measure")
ggplot(results_table, aes(variable, measure) ) + 
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values=c("dodgerblue4", "cadetblue")) +
  theme_minimal()
ggsave(paste0(plotDir, "results_50x_10x10cv_RF20-Acc-boxplot.pdf"),
       width=3 ,height=4)
store_tests[[1]]$finalModel

#### store importances ####
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



#### create RF for all samples ####
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
### big RF 
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
### small RF 20 
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

model_RF_best_100_TumiRs <- caret::train(CMS~., 
                                         data=rc_vst_BR[c(rownames(importances_df[1:100,])[
                                           which(rownames(importances_df[1:100,]) %in% tuMirs)],#grep("hsa", colnames(rc_vst_BR), value = T), 
                                                          "CMS")],
                                         method="ranger", 
                                         importance = 'impurity',
                                         metric="Kappa", 
                                         tuneGrid= bestParam_20, 
                                         num.trees = 2000, #2000 is best with 5 and 5 
                                         trControl= testingContr)
model_RF_all_TumiRs <- caret::train(CMS~., 
                                         data=rc_vst_BR[colnames(rc_vst_BR) %in% c(tuMirs, "CMS")],
                                         method="ranger", 
                                         importance = 'impurity',
                                         metric="Kappa", 
                                         tuneGrid= bestParam_20, 
                                         num.trees = 2000, #2000 is best with 5 and 5 
                                         trControl= testingContr)

##############################################################################################################
# Predictions
######### CMS from the test data
pred_rc_vst_BR_valid_RF <- predict(model_RF, 
                                   newdata = rc_vst_out[,grep("hsa",colnames(rc_vst_out))],
                                   type = "prob")
pred_rc_vst_BR_valid_RF$CMS <-  predict(model_RF, 
                                        newdata = rc_vst_out[,grep("hsa",colnames(rc_vst_out))],
                                        type = "raw")
# Print the confusion matrix for test data
cmat_RF <- confusionMatrix(pred_rc_vst_BR_valid_RF$CMS, na.omit(rc_vst_out$CMS))
summary(pred_rc_vst_BR_valid_RF)
cmat_RF

# ROC curve for test data, each CMS vs other
for (cms in 1:4){
refLabel <- factor(ifelse(na.omit(rc_vst_out$CMS)==paste0("CMS", cms), paste0("CMS", cms), 
                   "other"), levels=c("other",paste0("CMS", cms)))
result.roc <- roc(refLabel, pred_rc_vst_BR_valid_RF[,paste0("CMS", cms)]) # Draw ROC curve.
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft", 
     main=paste0("CMS", cms, " prob-thresh. (AUC)"))
result.coords <- coords(result.roc, "best", transpose = TRUE,
                        best.method="closest.topleft", 
                        ret=c("threshold", "accuracy"))
print(result.coords)#to get threshold and accuracy
}


######################################
### predictions on VU-data?
VU_rc_vst$Sample.ID <- row.names(VU_rc_vst)
VU_valid <- VU_rc_vst[,]#-grep("", row.names(VU_rc_vst)),
                     # grep("hsa", colnames(VU_rc_vst))]
pred_rc_vst_VU_RF <- predict(model_RF_best_all, newdata = VU_valid, type = "prob")
pred_rc_vst_VU_RF_TuMirs <- predict(model_RF_best_100_TumiRs, newdata = VU_valid, type = "prob")

# function to identify where first and second predicted prob are very close together
fCens = function(x, output) {
  ordered = sort(x, decreasing=T)
  censor = abs(ordered[1] - ordered[2]) <0.02
  censor
} 
pred_rc_vst_VU_RF$ordCens <- apply(pred_rc_vst_VU_RF,1, fCens)
pred_rc_vst_VU_RF$CMS <- predict(model_RF_best_all, 
                                    newdata = VU_valid, type = "raw")
pred_rc_vst_VU_RF$CMS_20 <- predict(model_RF_20_best_all, 
                                 newdata = VU_valid, type = "raw")
pred_rc_vst_VU_RF$CMS_TuMirs <- predict(model_RF_best_100_TumiRs, 
                                        newdata = VU_valid, type = "raw")

#pred_rc_vst_VU_RF$CMS[pred_rc_vst_VU_RF$ordCens] <- NA
rownames(pred_rc_vst_VU_RF) <- rownames(VU_valid)
pred_rc_vst_VU_RF$Sample.ID <-row.names(pred_rc_vst_VU_RF)
summary(pred_rc_vst_VU_RF$CMS)

## test if result differs when predicting only the primary samples
pred_VU_prim_CMS <- data.frame("CMS.prim"=predict(model_RF_best_all, 
                            newdata = VU_rc_vst[grep("primary",clinVU$sampleType),], 
                            type = "raw")) #predict only primary
pred_VU_prim_CMS$CMS.all <- pred_rc_vst_VU_RF$CMS[match(clinVU[grep("primary",clinVU$sampleType),
                                                              "sampleID"], #match to all
                                                     pred_rc_vst_VU_RF$Sample.ID)]
summary(pred_VU_prim_CMS$CMS.prim==pred_VU_prim_CMS$CMS.all) # no they don't differ

## plot VU predictions
data.frame("CMS"=pred_rc_vst_VU_RF[grep("[S,P][0-9]", rownames(pred_rc_vst_VU_RF)),
                                   "CMS"]) %>% 
  #[grep("^[P,M][0-9]", t_VU_rc_vst_surv$Sample.ID), "CMS"]) %>% #[!is.na(t_VU_rc_vst_surv$OS), "CMS"]) %>% 
  ggplot(aes(x = CMS, fill = CMS )) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=c(paletteCMS[1:4])) + #"#4F2776"
  ggtitle("CMS counts")

### add VU predictions to clinical annotation table
clinVU$CMS <- 
  pred_rc_vst_VU_RF[match(clinVU$sampleID,
   rownames(pred_rc_vst_VU_RF)),"CMS"] 
clinVU$CMS_20 <- 
  pred_rc_vst_VU_RF[match(clinVU$sampleID,
                          rownames(pred_rc_vst_VU_RF)),"CMS_20"] 

clinVU$Response_CR.y <- factor(clinVU$Response_CR.y, levels=
                                 c("CR", "PR", "SD",  "PD"), 
                               ordered=T)
# write.csv2(clinVU, row.names = F, 
#            "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/tables/VUclin-predictedCMS_RF-all_RF-20.csv")
#clinVU <- read.csv2( "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/tables/VUclin-predictedCMS_RF-all_RF-20.csv" )

library(wesanderson)
### plot: what are the tissue origins from the VU predictions?
ggplot(clinVU[which(grepl("CRC", clinVU$sampleType) & !is.na(clinVU$OS)),]
       %>% dplyr::count(CMS_20, LRcolon), 
       aes( CMS_20, n, fill = LRcolon )) +
  geom_bar(stat="identity", ) +
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Darjeeling1", 12,
                                         type = c( "continuous"))
                    ,na.value="#999999")+#values=c("#4F2776",paletteCMS[c(2,4)])) + #
  ggtitle("CMS counts CRC")
ggsave(paste0(plotDir,"VU_vst_RF_20_all-2000-10-2_4CMS_locColon.pdf"))

## overlap paired samples P&M
library(ComplexHeatmap)
concordPM <- clinVU[grep("^[P,M][0-9]", clinVU$sampleID), 
                     c("CMS", "sampleID", "sampleType", "patient")]
which(duplicated(clinVU$patient) & grepl("^S",clinVU$sampleID)) # no pairs identifiable in OS sample set 
concordPM$PM <- factor(sub("[0-9].*", "",concordPM$sampleID), 
                          levels=c("P", "M"), ordered=T)
concordPM$nCMS <- as.numeric(sub("CMS","", concordPM$CMS))
wConcPM <- concordPM %>% spread(PM,CMS)
wsConcPM <- wConcPM[-which(is.na(wConcPM$M)),-5]
wsConcPM <- left_join(wsConcPM, concordPM[grep("P", concordPM$sampleID),
                          c("patient", "CMS")])
wsConcPM$patient <- as.numeric(wsConcPM$patient)
wsConcPM$nCMS.M <- wsConcPM$nCMS
wsConcPM$nCMS.P<- as.numeric(sub("CMS","", wsConcPM$CMS))

rownames(wsConcPM) <- make.names(wsConcPM$patient, unique=T)
col_vec <- structure(as.character(wes_palette("Darjeeling1", 12,
                                   type = c( "continuous")))[1:7], 
          names = levels(droplevels(wsConcPM$sampleType)))

haConc <- HeatmapAnnotation(df=data.frame("sampleType"=wsConcPM$sampleType),
                            col = list( "sampleType"= col_vec),
                            na_col = "grey")
## plot concordance P&M
pdf(paste0(plotDir,
           "VU_vst_RF_all_best-2000-25-5_concordancePM_complexheatmap.pdf"),
    width=12, height=3)
draw(Heatmap(t(wsConcPM[,c("nCMS.M", "nCMS.P")]),
        name = "CMS", 
        col = structure(paletteCMS, names = c("1", "2", "3", "4")),
        cluster_rows = F,
        cluster_columns = T,
        column_split = wsConcPM$sampleType,
        top_annotation = haConc
        ))
dev.off()





##############################
##### on FFPE data 
pred_mir22_RF <- predict(model_RF_best_30,
                        newdata = scale(miR22_qn[,]),
                        type = "raw")
summary(pred_mir22_RF)
miR22.clin$CMSmir <- pred_mir22_RF
mRNA21.clin$CMSmRNA==miR22.clin$CMSmir
confusionMatrix(miR22.clin$CMSmir, mRNA21.clin$CMSmRNA)

###TODO: what to do with unmeasurable 
##TODO: PCA miR22 and mRNA21 datasets
miR22.clin$OS <- as.numeric(sub("overall survival (os): ","", fixed = T,
                     as.character(miR22.clin$Sample_characteristics_ch1_OS)))
miR22.clin$Event_OS <- ifelse( grepl("dead", as.character(miR22.clin$Sample_characteristics_ch1_OSevent)),
                                     1, 0)
miR22.clin$CMSmRNA <- mRNA21.clin$CMSmRNA

fit <- survfit(Surv(OS, Event_OS) ~ CMSmRNA,
               data = miR22.clin[,]) 

# visualize with survminer
pdf(paste0(plotDir,
           "GSE29623-mRNA21_RF_survOS.pdf"),
    width=6, height=7)
print( ggsurvplot(fit, data = miR22.clin[
  c( "OS", "Event_OS", "CMSmRNA")], #"OS", "Event_OS",
  risk.table = TRUE, pval = T,
  palette = paletteCMS[c(1,2,3,4)] ))
dev.off()

draw(Heatmap(t(data.frame(mRNA21.clin$CMSmRNA, miR22.clin$CMSmir)),
             name = "CMS", 
             col = structure(paletteCMS, 
                             names = c("CMS1", "CMS2", "CMS3", "CMS4")),
             cluster_rows = F,
             cluster_columns = T,
))



##############################
#### CPTAC-2
#miR_Cptac_vst_scale <- data.frame(t(scale(t(miR_Cptac_vst[,]))))
pred_cptec_RF <- predict(model_RF_20_best_all,
                        newdata = miR_Cptac_vst, 
                        type = "prob")
rownames(pred_cptec_RF) <- rownames(miR_Cptac_vst)
pred_cptec_RF$d2ndProb <- apply(pred_cptec_RF, 1, diffSecond)
pred_cptec_RF$d2ndClass <- apply(pred_cptec_RF[,1:4], 1, theSecond)
pred_cptec_RF$CMS <- predict(model_RF_20_best_all,
                             newdata = miR_Cptac_vst, #rpm_vst slightly better than raw_qn with RF20 (RF20 better than RF)
                             type = "raw")
summary(pred_cptec_RF$CMS)

clinSupp.cptac$CMSmir <- pred_cptec_RF$CMS[match(clinSupp.cptac$sample.label, 
                                                 rownames(pred_cptec_RF))]
clinSupp.cptac <- clinSupp.cptac[!is.na(clinSupp.cptac$sample.label), ]

confusionMatrix(factor(clinSupp.cptac$CMSmir), factor(clinSupp.cptac$CMS))

col_type <- structure(c("grey", "cadetblue", "bisque"),
                      names = levels(clinSupp.cptac$Mucinous))
clinSupp.cptac$LRsite <- factor(clinSupp.cptac$Subsite, 
                                labels=c( "Right", "Right", "Left",
                                         "Flex/Transv", "Left",
                                         "Flex/Transv", "Flex/Transv"))
clinSupp.cptac$LRsite <- factor(clinSupp.cptac$LRsite, 
                                levels=c("Right", 
                                         "Flex/Transv", "Left"), ordered=T)
clinSupp.cptac$TumorPurity
col_site <- structure(c("bisque","cadetblue", "dodgerblue4"),
                      names = levels(clinSupp.cptac$LRsite))
col_msi <- structure(c("grey","darkgoldenrod2", "bisque", "cadetblue"),
                      names = levels(clinSupp.cptac$MSI_PCR_Result))
col_ums <- structure(c("grey", "dodgerblue4", "cadetblue", "darkgoldenrod2"),
                     names = levels(clinSupp.cptac$UMS))
col_dec <- structure(c("grey","dodgerblue4", "white"),
                     names = levels(clinSupp.cptac$Vital.Status))

haCptac <- HeatmapAnnotation(df=clinSupp.cptac[order(clinSupp.cptac$CMSmir),
                                               c(  "Stage", "LRsite",
                                                   "Vital.Status",
                                       "TumorPurity", "MSI_PCR_Result")],
                                       col = list( Stage=col_stage,
                                                   LRsite=col_site,
                                                   Vital.Status=col_dec,
                                                   TumorPurity=col_fun, 
                                                   MSI_PCR_Result=col_msi,
                                                   UMS=col_ums),
                                       na_col = "grey")
clinSupp.cptac$conc <- ifelse(as.character(clinSupp.cptac$CMSmir)==
                                as.character(clinSupp.cptac$CMS),
                              "concordant", "discordant")
pdf(paste0(plotDir,
           "CPTAC2_vst_scale_RF_all_scale_mini20-2000-10-2_concordancePM_complexheatmap.pdf"),
    width=20, height=3)
draw(Heatmap(t(clinSupp.cptac[order(clinSupp.cptac$CMSmir),c("CMSmir", "CMS")]),
             name = "CMS", 
             col = structure(paletteCMS, names = c("CMS1", "CMS2", "CMS3", "CMS4")),
             cluster_rows = F,
             cluster_columns = T,
             column_split = clinSupp.cptac$conc[order(clinSupp.cptac$CMSmir)],
             top_annotation = haCptac
))
dev.off()
summary(clinSupp.cptac$Vital.Status)
clinSupp.cptac$Event_OS <- ifelse(clinSupp.cptac$Vital.Status=="Deceased", 1, 0)








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

# confusionM.READ <- data.frame("CMS1.mRNA"=c(4,0,0,0),
#                               "CMS2.mRNA"=c(0,31,1,14),
#                               "CMS3.mRNA"=c(1,2,8,1),
#                               "CMS4.mRNA"=c(0,2,0,25))
# confusionM.CPTAC <- data.frame("CMS1.mRNA"=c(7,1,3,2),
#                                "CMS2.mRNA"=c(0,25,0,8),
#                                "CMS3.mRNA"=c(2,4,6,4),
#                                "CMS4.mRNA"=c(1,4,0,17))
#rownames(confusionM.CPTAC)<-paste0("CMS", 1:4, "miR") 
confusionM.READ <- confusionMatrix(pred_read_RF$CMS, 
                                   factor(miR_READ_vst_BR$CMS.lv[selREAD]))[["table"]]
rownames(confusionM.READ) <- paste0("CMS", 1:4, "miR") 
colnames(confusionM.READ) <- paste0("CMS", 1:4, "mRNA") 

pheatmap(as.matrix(confusionM.READ), 
         cluster_rows = F,
         cluster_cols = F, 
         color= brewer.pal(n = 9, name = "Purples"),
)

## plot: what are the READ predictions
data.frame("CMS"=pred_read_RF$CMS) %>%#
  ggplot(aes(x = CMS, fill = CMS )) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=c(paletteCMS[]),na.value="#999999") + #"#4F2776"
  ggtitle("CMS counts")
data.frame("CMS_20"=pred_read_RF$CMS_20) %>%#
  ggplot(aes(x = CMS_20, fill = CMS_20 )) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=c(paletteCMS[]),na.value="#999999") + #"#4F2776"
  ggtitle("CMS counts")

##plot: what is the READ concordance
concREAD <- data.frame("CMS"=factor(c(miR_READ_vst_BR[selREAD , "CMS.lv"],
                       as.character(pred_read_RF$CMS), as.character(pred_read_RF$CMS_20))),
           "pat"=c(rownames(miR_READ_vst_BR[selREAD,]),
           rownames(miR_READ_vst_BR[selREAD,]),
           rownames(miR_READ_vst_BR[selREAD,])),
           "type"= factor(c(rep("mRNA", length(miR_READ_vst_BR[selREAD,1])),
                     rep("miR", length(miR_READ_vst_BR[selREAD,1])),
                         rep("miR_RF_20", length(miR_READ_vst_BR[selREAD,1])))))
# ggplot(concREAD, aes(x=pat, y=type, fill= CMS)) + 
#   geom_tile()+
#   theme_bw()+
#   scale_fill_manual(values=c(paletteCMS[],"#999999"),na.value="#999999")+
#   theme(axis.text=element_text(angle = 90), 
#         axis.title= element_text(angle = 90))
# ggsave(paste0(plotDir,
#            "READ_RFkappaCaret_mini20_concordancePM.pdf"),
#     width=12, height=3)

#### overlap READ samples mRNA vs mirNA
concREAD$CMS <- as.numeric(sub("CMS", "", concREAD$CMS))
concREAD$CMS <- factor(paste0("CMS", concREAD$CMS))

wideREAD <- reshape(concREAD, direction="wide", 
                    idvar="pat", timevar = "type", v.names = "CMS")

library("ggalluvial")
ggplot(as.data.frame(wideREAD),
       aes(axis1 = CMS.mRNA, axis3 =CMS.miR, axis2 = CMS.miR_RF_20)) +
  geom_alluvium(aes(fill = CMS.mRNA), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  scale_x_discrete(limits = c("mRNA", "miR_RF20", "miR"), expand = c(.05, .05)) +
  scale_fill_manual(values = paletteCMSn) +
  theme_minimal() + 
  ggtitle("READ mRNA vs. miR CMS")
ggsave(paste0(plotDir, "READ_alluvial_RF20.pdf"))

##make dataframe for heatmap 
concREAD$CMS <- as.numeric(sub("CMS", "", concREAD$CMS))
wideREAD <- reshape(concREAD, direction="wide", 
                    idvar="pat", timevar = "type", v.names = "CMS")
wideREAD$d2ndProb <- pred_read_RF$d2ndProb
wideREAD$d2ndProb_20 <- pred_read_RF$d2ndProb_20
wideREAD$d2ndClass <- pred_read_RF$d2ndClass_20

wideREAD <- wideREAD[order(wideREAD$CMS.mRNA),]
wideREAD$type <- factor(clinical.read$primary_diagnosis[match(wideREAD$pat, 
                                                                   clinical.read$submitter_id)])
wideREAD$stage <- factor(sub("[A,B,C]$","" ,clinical.read$ajcc_pathologic_stage[match(wideREAD$pat, 
                                                         clinical.read$submitter_id)]))
wideREAD$conc <- ifelse(wideREAD$CMS.mRNA==wideREAD$CMS.miR,
                       "concordant", "discordant")
rownames(wideREAD) <- make.names(wideREAD$pat, unique=T)

wideREAD$conc[which(is.na(wideREAD$conc))] <- "NA"
wideREAD$purity <- purity.ABS.read[match(rownames(wideREAD),
                                        gsub("-",".",purity.ABS.read$TCGA.patID)),
                                  "purity"]
wideREAD$Ca.DNA <- purity.ABS.read[match(rownames(wideREAD),
                                        gsub("-",".",purity.ABS.read$TCGA.patID)),
                                  "Cancer.DNA.fraction"]
wideREAD$mRNApredPval <- log10(data.frame("mRNApredPval"=CMS_samples$min_Pval[match(rownames(wideREAD),
                                         gsub("-",".",CMS_samples$SampleId))])[,"mRNApredPval"])
summary(wideREAD$d2ndProb)
summary(wideREAD$d2ndProb_20)
wideREAD$uncertainty <- factor(ifelse(wideREAD$d2ndProb < 0.086 | wideREAD$d2ndProb_20 < 0.12, 
                        "two calls", "clear winner"))
wideREAD$d2ndClass[which(wideREAD$uncertainty=="clear winner")] <- NA

col_stage <- structure(c("bisque", "cadetblue", "dodgerblue4","chocolate4"),
                      names = levels(wideREAD$stage))
col_type <- structure(c("dodgerblue4", "chocolate4", 
                        "cadetblue", "darkgoldenrod2", "bisque"),
                      names = levels(wideREAD$type))
col_bin <- structure(c("white", "grey"),
                      names = levels(wideREAD$uncertainty))
col_fun <- colorRamp2(c(0, 1), c("white", "dodgerblue4"))
col_fun_neg <- colorRamp2(c(0, -5), c("white", "dodgerblue4"))

haConc <- HeatmapAnnotation(df=wideREAD[,c( "type", "stage",
                                          "purity", "d2ndProb",
                                          "d2ndProb_20", #"uncertainty",
                                          "mRNApredPval")],
                            col = list( stage=col_stage,
                                        type=col_type,
                                        purity=col_fun,
                                        d2ndProb=col_fun,
                                        d2ndProb_20=col_fun,
                                        #uncertainty=col_bin,
                                        mRNApredPval=col_fun_neg),
                            na_col = "grey")
## plot concordance P&M
pdf(paste0(plotDir,
           "READ_RF_all_20-2000-10-2_wt-2ndCLass_RF_20.pdf"),
    width=22, height=5)
draw(Heatmap(t(wideREAD[,c("CMS.mRNA", "CMS.miR", "CMS.miR_RF_20", "d2ndClass")]),
             col = structure(paletteCMS, names = c("1", "2", "3", "4")),
             cluster_rows = F,
             cluster_columns = F, 
             column_split=wideREAD$uncertainty,
             top_annotation = haConc,
             column_names_gp = gpar(fontsize = 6)
))
dev.off()


ggplot(wideREAD, aes(-mRNApredPval, d2ndProb) ) + 
  geom_point(aes(col=factor(wideREAD$CMS.miR)))+
  facet_grid(wideREAD$CMS.miR) +
  scale_colour_manual(values=paletteCMS, name="CMS") +
  scale_fill_manual(values=paletteCMS, name="CMS") +
  theme_minimal() +
  geom_smooth(aes(color=factor(wideREAD$CMS.miR), 
                  fill=factor(wideREAD$CMS.miR)), 
              method = "lm", show.legend = F) +
  geom_jitter(aes(color=factor(wideREAD$CMS.miR), 
                  fill=factor(wideREAD$CMS.miR)), width = 0.1)+
  labs(y="distance to 2nd class miR", x = "-log10(pval) pred. class mRNA")+
  stat_cor()
ggsave(paste0(plotDir, "bestRF_all_mRNApval-vs-classProbDifference-by-CMS.pdf"),
       width=4, height=6)





#### on independent 3 samples from GSE121842
predict(model_RF, newdata = miRC42_vst)


### store a good model
# selCaretRF <- model_RF
# sele_pred_rc_vst_VU_RF <- pred_rc_vst_VU_RF



##################################################################################################################################      
####################### Check importance of features ####################### 

##for caret R:
N_imp_RF <- importance(model_RF_best$finalModel)[,"MeanDecreaseGini"] # for caret RF
N_imp_RF <- varImp(model_RF)
N_imp_RF_df <- tibble(variable = names(N_imp_RF), 
                      importance = N_imp_RF ) %>% 
  dplyr::arrange(-importance) 
##for caret ranger:
N_imp_RF <- varImp(model_RF_best)$importance
N_imp_RF_df <- tibble(variable = rownames(N_imp_RF), 
                    importance = N_imp_RF$Overall ) %>% 
  dplyr::arrange(-importance) 

"Purples"
"YlOrBn"
"PuBu"
"PuRd"
"Greens"
head( N_imp_RF_df, 30 ) %>%
  ggplot(aes(x = reorder(variable, -importance), 
             y = importance, fill = importance)) +
  geom_bar(stat = "identity" ) +
  scale_fill_gradientn(colours=brewer.pal(9, "Purples")[3:9]) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = "none")




##########################################
########### survival effects?? ###########

# get VU suvival data

#for(sel.mir in c("hsa.mir.429", "hsa.mir.552", 
#                 "hsa.mir.200a", "hsa.mir.200b", "hsa.mir.141" )){
sel.mir <- "hsa.mir.552"
cut625 <- median(VU_rc_vst[grep("primary", clinVU$sampleType), sel.mir])
bin.sel.mir <- ifelse(VU_rc_vst[,sel.mir] >=cut625, "low", "high")
clinVU$bin.sel.mir <- factor(bin.sel.mir)
#}

# install.packages('survival')
# install.packages('survminer')
# BiocManager::install("RTCGA.clinical") # data for examples
library(survival)
library(survminer)
library(RTCGA.clinical)
#event: in OS analyisi 0=alive, 1=dead, people with 0, meaning without Death/event will be censored 
clinVU$Event_PFS <- ifelse(clinVU$PFS_1==clinVU$OS & clinVU$Event_OS==0, 0, 1) 
#if OS=PFS patient alive at end of study or dropped out from any cause -> censor that 3 patients
clinVU$Event_PFS[is.na(clinVU$PFS_1)] <- NA

fit_20 <- survfit(Surv(OS, Event_OS) ~ CMS_20, #OS, Event_OS
               data = clinVU[#grepl("Oxaliplatin", clinVU$Scheme_for_analysis),])#
                 which(#clinVU$CMS_20 %in% c("CMS2", "CMS4")
                 #grepl("tachroon", clinVU$Syn_metachroom)
                # &
                 grepl("prim",#allMet #"[S][0-9]", #only survival cohort
                       clinVU$sampleType)),]) #[grep("P[0-9]", t_VU_rc_vst_surv$Sample.ID),]
# visualize with survminer
pdf(paste0(plotDir,
           "VU_vst_RF_all_2000-10-2_noclwt_survOS.pdf"),
    width=6, height=7)
print( ggsurvplot(fit, data = clinVU[#grepl("Oxaliplatin", clinVU$Scheme_for_analysis),#
  which(#clinVU$CMS_20 %in% c("CMS2", "CMS4") &
                                                 grepl("prim", 
                                                      clinVU$sampleType)), #
                                                 #grepl("tachroon", clinVU$Syn_metachroom)),
                              c( "OS", "Event_OS", "CMS")], #"OS", "Event_OS",
           risk.table = TRUE, pval = T, xlim=c(0,1500),
           palette = paletteCMS[c(2,4)] ))#"#4F2776" #c("#2a7886","#79bac1" ))# 
dev.off()

 
### get TCGA survival rectal dataset
rc_vst_surv <- cbind(clinical.read, "CMS"=factor(miR_READ_vst_BR$CMS.lv[
  match(clinical.read$submitter_id, rownames(miR_READ_vst_BR))], ordered=T))
rc_vst_surv$OS <- clinical.read$days_to_last_follow_up
rc_vst_surv$OS[which(is.na(rc_vst_surv$OS))] <- rc_vst_surv$days_to_death[which(is.na(rc_vst_surv$OS))]
rc_vst_surv$Event_OS <- ifelse(clinical.read$vital_status=="Alive", 0, 1) #Surv takes 0 as alive and 1 as dead
rc_vst_surv$tissue_or_organ_of_origin[which(rc_vst_surv$CMS=="CMS1")]
fit <- survfit(Surv(OS, Event_OS) ~ CMS,
               data = rc_vst_surv)

# TCGA.survInfo <- survivalTCGA(READ.clinical, extract.cols = "admin.disease_code") 
# TCGA.survInfo$CMSpred <- pred_read_RF[match(TCGA.survInfo$bcr_patient_barcode, 
#                            rownames(miR_READ_vst_BR)), "CMS"]
# rc_vst_surv <- cbind( miR_READ_vst_BR, TCGA.survInfo[
#   match(rownames(miR_READ_vst_BR), TCGA.survInfo$bcr_patient_barcode),])
# rc_vst_surv$bin.sel.mir <- factor(ifelse(rc_vst_surv$hsa.mir.552 >= median(rc_vst_surv$hsa.mir.552),
#                                                         "low", "high"))
# rc_vst_surv$Event_OS <- ifelse(rc_vst_surv$patient.vital_status==1, 1, 0)
# fit <- survfit(Surv(times, Event_OS) ~ CMSpred,
#                data = rc_vst_surv)
# visualize with survminer
pdf(paste0(plotDir,
           "SURV_READ-cms-lv.pdf"),
    width=6, height=7)
print( ggsurvplot(fit, data = rc_vst_surv, risk.table = TRUE, pval = T,
           palette =  paletteCMS, surv.median.line = "hv", pval.method = T ) )
dev.off()

###TCGA COAD data 
rc_vst_surv <- cbind(TCGA.survInfo, "CMS"=factor(miR_COAD_vst$CMS.lv[
    match(TCGA.survInfo$submitter_id, rownames(miR_COAD_vst))], ordered=T))
TCGA.survInfo <- survivalTCGA(COAD.clinical, extract.cols = "admin.disease_code") 
rpm_surv <- cbind( rpm, TCGA.survInfo[
  match(sub("-01.*","",rownames(rpm)), 
        TCGA.survInfo$bcr_patient_barcode),])
rc_vst_surv$OS <- clinical.coad$days_to_death
rc_vst_surv$OS[which(is.na(rc_vst_surv$OS))] <- rc_vst_surv$days_to_last_follow_up[which(is.na(rc_vst_surv$OS))]
rc_vst_surv$Event_OS <- ifelse(clinical.coad$vital_status=="Alive", 0, 1) #Surv takes 0 as alive and 1 as dead
summary(factor(rc_vst_surv$tissue_or_organ_of_origin[which(rc_vst_surv$CMS=="CMS1")])) #to check if CMS cbind to clin seems correct
fit <- survfit(Surv(OS, Event_OS) ~ CMS,
               data = rc_vst_surv[grep(" iii?[a,b,c]?",rc_vst_surv$tumor_stage),])
# visualize with survminer
pdf(paste0(plotDir,
           "SURV_COAD-cms-rf.cl.pdf"),
    width=6, height=7)
ggsurvplot(fit, data = rc_vst_surv[grep(" iii?[a,b,c]?",rc_vst_surv$tumor_stage),], risk.table = TRUE, pval = T,
           palette =  paletteCMS, surv.median.line = "hv" ) 
dev.off()


# Fit a Cox proportional hazards model
rc_vst_surv$CMS1 <- factor(ifelse(rc_vst_surv$CMS=="CMS1", "CMS1", "other"))
rc_vst_surv$CMS2 <- factor(ifelse(rc_vst_surv$CMS=="CMS2", "CMS2", "other"))
rc_vst_surv$CMS3 <- factor(ifelse(rc_vst_surv$CMS=="CMS3", "CMS3", "other"))
rc_vst_surv$CMS4 <- factor(ifelse(rc_vst_surv$CMS=="CMS4", "CMS4", "other"))

surv_object <- Surv(time = rc_vst_surv$OS, 
                    event = rc_vst_surv$Event_OS)
fit.coxph <- coxph(surv_object ~ #CMS1 + CMS2 + CMS3 + CMS4,
                   hsa.mir.625 +
                   hsa.mir.592 + hsa.mir.218+
                   hsa.mir.552 + hsa.mir.31 +
                     hsa.mir.615 + hsa.mir.92b
                   + hsa.mir.143 + hsa.mir.942
                   + hsa.mir.375 + hsa.mir.141,
                   data = rc_vst_surv)
ggforest(fit.coxph, data = rc_vst_surv)


##plot density distribution
ggplot( rc_vst_BR, aes( x= hsa.mir.31) ) + 
  geom_density( aes( group=CMS, colour=CMS, fill=CMS ), alpha=0.3 ) +
  scale_fill_manual( values = paletteCMS ) +
  scale_color_manual( values = paletteCMS ) +
  theme_minimal( )

### colors CMS
# paletteCMS=c("#DAA520", 
#   "#1874CD", 
#   "#CD2990", 
#   "#458B74")
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




WangXiaoWang <- c("hsa.mir21", "hsa.mir26", "hsa.mir31", "hsa.mir141", "hsa.mir145", "hsa.mir196", 
  "hsa.mir200",
  "hsa.mir222", "hsa.mir29b", "hsa.mir503", "hsa.mir322", "hsa.mir195", "hsa.mir122a", 
  "hsa.mir146a",
  "hsa.mir493",
  "hsa.mir.429", "hsa.mir.9", "hsa.mir.21", "hsa.mir.145", 
  "hsa.mir.155", "hsa.mir.192", "hsa.mir.374")




