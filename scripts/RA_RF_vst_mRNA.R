# Libraries needed
library(MASS)
library(tidyverse)
library(randomForest)
library(caret)
library(nnet)
library(pROC)
#install.packages('e1071', dependencies=TRUE)
set.seed(67)

##################################################################################################################################      
# Open the csv file
rpm <- read.table("data/cleaned/RPM_prob_p005.csv")
rc <- read.table("data/cleaned/RC_prob_p005.csv")

# Only use the features with read count mean > 1
col_means <- colMeans(rc[-c(1,2)])
feature_bigger1 <- min( which( sort( col_means )>1)) # First feature with a mean > 1
feature_names_bigger1 <- names( sort( col_means )[feature_bigger1 :length(sort(col_means))]) 
feature_names <- c( "CMS", feature_names_bigger1 )

rpm <- rpm[,feature_names]
rpm$Sample.ID <- row.names(rpm)
rc_vst$Sample.ID <- row.names(rc_vst)
# keep out of bag samples for performance testing 
out_Test <- sample(1:339, 339/5)
rpm_out <- rpm[out_Test, ]
rpm <- rpm[-out_Test,]
# save(selectedRF,file = "selRF.RData")
# save(model_RF,file = "lastModelRF.RData")
# save(rpm_out, file= "rpm_out.RData")
# save(rpm, file="rpm_input.RData")
selectedRF <- load(file = "selRF.RData")
model_RF <- load(file = "lastModelRF.RData")
rpm_out <- load(file= "rpm_out.RData")
rpm <- load( file="~/projects/miRNA_mCRC/miRNA classifier/rpm_input.RData")
TCGARnaseqRaw <- read.table("~/projects/miRNA_mCRC/TCGA_COADREAD_mRNA_raw_gene_counts.tsv")


### normalize raw counts (rc) with varianceStabilizingTransformation 
library(DESeq2)
rc_CMS <- rc$CMS
rc <- rc[,-1]
rc_vst <- varianceStabilizingTransformation( t( round( rc[, 
                            which( colMeans( rc ) > 1)],0))) #toss low count features
##OR for mRNA:
# rc_vst <- varianceStabilizingTransformation( 
#   as.matrix( round( TCGARnaseqRaw[,which( 
#     colMeans( TCGARnaseqRaw ) > 1)],0)))
rc_vst <- as.data.frame( t( rc_vst[,] ) )
#rc_vst$CMS <- patCMSmRNA
rc_vst$CMS <- rc_CMS
rc_vst <- rc_vst[-out_Test, ] #toss the out of training samples 
  
##################################################################################################################################      
####################### RF training with Cross validation ####################### 

##### RF with caret ######
##preprocessing: visualize

##identifying low info features, plot their correlation
nzv <- nearZeroVar(rc_vst, saveMetrics= TRUE)
summary(nzv$nzv) 
descrCor <-  cor(rc_vst[,-c(which(nzv$nzv),20532)])
library(corrplot)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75, exact=T )
filteredDescr <- rc_vst[,-c(highlyCorDescr,20532)]

## visualization of class distances
centroids <- classDist( filteredDescr, rc_vst$CMS, pca = T, keep = 30 )
distances <- predict(centroids, rc_vst)
distances <- as.data.frame(distances)
head(distances)
xyplot(dist.CMS3 ~ dist.CMS4,
       data = distances, 
       groups = rc_vst$CMS, 
       auto.key = list(columns = 2, col = paletteCMS ), 
       col = paletteCMS)

### PCA plot ###
#install.packages("factoextra")
library(factoextra)
prin_comp <- prcomp(rc_vst[,1:20531], scale. = T)
fviz_eig(prin_comp, type="lines" ) #plot variances
pr_var <- round( prin_comp$sdev^2, 5 )
prop_varex <- round( pr_var/sum(pr_var), 5)
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
var.coord <- t(apply(prin_comp$rotation, 1, var_coord_func, prin_comp$sdev)) 
head(var.coord[order(var.coord[,"PC3"]), "PC3"]) #the top features
tail(var.coord[order(var.coord[,"PC3"]), "PC3"])

library(caret)
trans <- preProcess(rc_vst[
                        ,1:20531], 
                    method=c(#"BoxCox",
                             #"center","scale", #makes no difference in vst
                             "pca"))
PC <- predict(trans, rc_vst[,1:20531])
ggplot(PC,aes(x=PC1,y=PC2,
              colour=rc_vst$CMS[
                ])) +
  geom_point() +
  scale_color_manual(values=paletteCMS)


#### tSNE plot ####
set.seed(9)  
library(Rtsne)
tsne_model_1 <- Rtsne(as.matrix(rc_vst[,1:20531]), check_duplicates=FALSE, 
                     pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_1 <- as.data.frame(tsne_model_1$Y)  
ggplot(d_tsne_1, aes(x=V1, y=V2, colour=rc_vst[, "CMS"])) +  
  geom_point() +
  xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = paletteCMS)


#### plot kmeans class labels on PCA #####
kM <- kmeans(rc_vst[, 1:20531], 4, nstart = 20)
kM$cluster <- as.factor(kM $cluster)
summary(kM$cluster)
ggplot(rc_vst, aes(PC$PC1, 
                 PC$PC2, 
                 color = kM$cluster)) + geom_point()
#### heatmap ####
heatmap(log10(t(rc_vst[,-c(492,493)]+1))) 
#TODO: plot kmeans class labels on heatmap

#TODO: check if the TCGA data was normalized or should be, test everything with rlog,
#test everything with mRNA


#### plot cluster dendrogram ####
hc <- hclust(dist(rc_vst[,1:20531]), "ward.D") #TODO: test other distance methods

# customized plot
# make color code
colrc_vst.CMS <- factor(rc_vst$CMS, labels = c("darkgoldenrod1", 
                                         "blue4", "deeppink3",
                                         "mediumseagreen"))
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- as.character(colrc_vst.CMS[which(rc_vst$Sample.ID == a$label)])
    attr(n, "nodePar") <- c(a$nodePar,
                            pch="o",
                            col = labCol) }
  return(n) }
hcd <- as.dendrogram(hc)
# using dendrapply to get color-annotation for each node
clusDendro = dendrapply(hcd, colLab)
plot(clusDendro)


#### testing accuracy of RF vs SVM ####
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
# SVM
fit.svm <- train(CMS~., data=rc_vst[,-c(highlyCorDescr,493)], method="svmRadial", 
                 metric="Kappa", trControl=control)
# SVM
fit.svm.acc <- train(CMS~., data=rc_vst[,-c(highlyCorDescr,493)], method="svmRadial", 
                 metric="Accuracy", trControl=control)
# Random Forest
fit.rf <- train(CMS~., data=rc_vst[,-c(highlyCorDescr,493)], method="rf", 
                metric="Kappa", trControl=control)
# Random Forest
fit.rf.acc <- train(CMS~., data=rc_vst[,-c(highlyCorDescr,493)], method="rf", 
                metric="Accuracy", trControl=control)
# compare accuracy of models
results <- resamples(list(svm.kappa=fit.svm, svm.acc=fit.svm.acc, 
                          rf.kappa=fit.rf, rf.acc=fit.rf.acc))
summary(results)
dotplot(results)




###### RF with randomForest ######
#test optimal number of features
#mtry is no of Variables randomly chosen at each split
oob.err <- double(11)
test.acc <- double(11)
i=1
for(mtry in seq(22,222, by= 20)) {
  rf=randomForest(CMS ~ . - Sample.ID, data = rc_vst_train , mtry=mtry,
                  ntree=1000, cutoff=c(.2,.4,.2,.2)) 
  oob.err[i] = mean(rf$err.rate) #Error of all Trees fitted
  
  pred<-predict(rf,rc_vst_valid) #Predictions on Test Set for each Tree
  test.acc[i]= confusionMatrix(pred, rc_vst_valid$CMS)[["overall"]][["Accuracy"]]
  
  cat(mtry," ") #printing the output to the console
  i<- i+1
}
matplot( seq(22,222, by= 20) , cbind(oob.err,test.acc), 
         pch=19 , col=c("red","blue"),
         type="b",ylab="Accuracy",
         xlab="Number of Predictors Considered at each Split")
legend("right",legend=c("Out of Bag Error",
                           "Accuracy"),pch=19, 
       col=c("red","blue"))


# Make an empty variable
cv_comparison = NULL
list_RF <- vector("list", 10)
names(list_RF) <- paste0("model_cv_", 1:10)
for (i in 1:10){
  ###added cross validation because due to oversampling we cant assume independence of samples in rc_vst, so oob error rate might be too optimistic
  #https://stats.stackexchange.com/questions/283760/is-cross-validation-unnecessary-for-random-forest
  selectTrain <- rc_vst  %>% group_by(CMS) %>% #or rc_vst[which(rc_vst$CMS != "CMS2"),] #for trying how well separation goes without CMS2
    dplyr::sample_frac(0.7) %>% 
    dplyr::select(Sample.ID)
  rc_vst_train <- rc_vst[which(rc_vst$Sample.ID %in% selectTrain$Sample.ID), ]
  rc_vst_train$CMS <- factor(rc_vst_train$CMS)
  rc_vst_valid <- rc_vst[-which(rc_vst$Sample.ID %in% selectTrain$Sample.ID), ]
  
  # # Do OVERSAMPLING, duplicate all samples from CMS 1, 3 and 4 from the training set
  # rc_vst_train_CMS134 <- rc_vst_train_pre_over %>%
  #   dplyr::filter(CMS == "CMS3" | CMS == "CMS4")
  # rc_vst_train <- dplyr::bind_rows(rc_vst_train_pre_over, rc_vst_train_CMS134)

  ####################### Make the RF model ####################### 
  ##with randomForest
  model_RF <- randomForest(CMS ~ . - Sample.ID, importance = TRUE,
                           replace = F, #less likely to overuse CMS2?
                           proximity = TRUE, data = rc_vst_train, 
                           ntree = 2000,
                           mtry = 100, 
                           keep.inbag=TRUE, 
                           #cutoff=c(0.3,0.3,0.2,0.2))
                           classwt=c(0.3,0.9,0.6,0.15))
                           #sampsize=c(30,15,30,30), strata=rc_vst_train$CMS)
  ##with  caret
  model_RF <- caret::train()
  list_RF[[paste0("model_cv_", i)]] <- model_RF
  model_RF
  
  # Predict CMS from the validation data
  pred_rc_vst_valid_RF <- predict(model_RF, newdata = rc_vst_valid)
  
  # Print the confusion matrix
  cmat_RF <- confusionMatrix(pred_rc_vst_valid_RF, rc_vst_valid$CMS)
  cmat_RF
  
  # ####################### Make the RF_no0 model #######################
  # model_RF_no0 <- randomForest(CMS ~ . - Sample.ID, data = rc_vst_train_no0)
  # #model_RF_no0
  # list_RF <- append( model_RF_no0, list_RF )
  # 
  # # Predict CMS from the validation data
  # pred_rc_vst_valid_RF_no0 <- predict(model_RF_no0, rc_vst_valid_no0)
  # 
  # # Print the confusion matrix
  # cmat_RF_no0 <- confusionMatrix(pred_rc_vst_valid_RF_no0, rc_vst_valid_no0$CMS)
  # cmat_RF_no0
  # 
  # Make a row with the accuracy of the models, how well CMS 3 is predicted, 
  # and the amount of CMS3 in the training set (to check whether that matters, doesnt seem so)
  new_row <- data_frame("model_RF_acc" = cmat_RF$overall[1], 
                        # "model_RF_no0_acc" = cmat_RF_no0$overall[1], 
                        "specifiticy_CMS1" = cmat_RF$byClass[1, 1], 
                        "specifiticy_CMS2" = cmat_RF$byClass[2, 1], 
                        "specifiticy_CMS3" = cmat_RF$byClass[3, 1], 
                        "specifiticy_CMS4" = cmat_RF$byClass[4, 1], 
                        )
  
  # Add this row to cv_comparison
  cv_comparison <- rbind(cv_comparison, new_row)
  ##################################################################################################################################      
}
cv_comparison 
colMeans(cv_comparison)
which(cv_comparison$`model_RF_acc`==max(cv_comparison$`model_RF_acc`))
model_RF <- list_RF$model_cv_9
## on VU-data?
t_VU_rc_vst_df$Sample.ID <- row.names(t_VU_rc_vst_df)
pred_rc_vst_VU_RF <- predict(model_RF, 
                              newdata = t_VU_rc_vst_df[grep("^S",
                                                         row.names(t_VU_rc_vst_df)),
                                                    ])
summary(pred_rc_vst_VU_RF)

### store a good model
# selectedRF <- model_RF
# sele_pred_rc_vst_VU_RF <- pred_rc_vst_VU_RF
# Predict CMS from the out of box data
pred_rc_vst_valid_RF <- predict(model_RF, newdata = rc_vst_out)
# Print the confusion matrix for external validation
cmat_RF <- confusionMatrix(pred_rc_vst_valid_RF, rc_vst_out$CMS)
summary(pred_rc_vst_valid_RF)

cmat_RF

###how well did the trees improve the model
oob.error.data <- data.frame(Trees=rep(1:nrow(model_RF$err.rate), 
                                        times=5),
                             Type=rep(c("OOB",paste0("CMS", c(1,2,3,4))),
                                        each=nrow(model_RF$err.rate)),
                             Error=c(model_RF$err.rate[,"OOB"],
                                     model_RF$err.rate[,"CMS1"],
                                     model_RF$err.rate[,"CMS2"],
                                     model_RF$err.rate[,"CMS3"],
                                     model_RF$err.rate[,"CMS4"]))
ggplot(data=oob.error.data,
       aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type)) +
  scale_colour_manual(values=c(paletteCMS, "#999999"))

MDSplot(model_RF, rc_vst$CMS, palette= paletteCMS, 3)

# ROC metric 
#library(pROC)
pred <- multiclass.roc(model_RF$predicted, model_RF$votes)
rs <- pred[['rocs']]
plot.roc(rs[[6]][[1]], col="#336633")
sapply(2:length(rs),function(i) lines.roc(rs[[i]][[1]],col=
                                            c("#336633", "#FF3333",
                                              "#669900", "#660066",
                                              "#003333", "#330000")[i], 
                                          lty=8-i))

##################################################################################################################################      
####################### Check importance of features ####################### 
N_imp_RF <- importance(model_RF)[,"MeanDecreaseGini"]
summary(importance(model_RF))
N_imp_RF_df <- tibble(variable = names(N_imp_RF), 
                    importance = N_imp_RF ) %>% 
  dplyr::arrange(-importance) 

"Purples"
"YlOrBn"
"PuBu"
"PuRd"
"Greens"


head( N_imp_RF_df, 20 ) %>%
  ggplot(aes(x = reorder(variable, -importance), 
             y = importance, fill = importance)) +
  geom_bar(stat = "identity" ) +
  scale_fill_gradientn(colours=brewer.pal(9, "Purples")[3:9]) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = "none")

##TODO: check out inTrees package to get the rules from the forest 
#TODO: really TMM normalized?
#get best RFmodel and it's importance
#feature selection process?? 
#solve class imbalance problem

### survival effects??
# get suvival data
clinVU <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/191104patientsforCMS_D6response.csv")
#match sample names to patient ID
clinVU$sampleID <- rownames(t_VU_rc_vst_df)[
  match(sub(" .*","",clinVU$ID),
  sub("\\..*","",rownames(t_VU_rc_vst_df)) )]
t_VU_rc_vst_df$OS <- clinVU$OS[match(t_VU_rc_vst_df$Sample.ID, clinVU$sampleID)]
t_VU_rc_vst_df$Event_OS <- clinVU$Event_OS[match(t_VU_rc_vst_df$Sample.ID, 
                                              clinVU$sampleID)]

#for(sel.mir in c("hsa.mir.429", "hsa.mir.552", 
#                 "hsa.mir.200a", "hsa.mir.200b", "hsa.mir.141" )){
cut625 <- median(t_VU_rc_vst_df[!is.na(t_VU_rc_vst_df$OS), "hsa.mir.592"])
t_VU_rc_vst_surv <- t_VU_rc_vst_df %>% dplyr::mutate(bin.sel.mir = 
                                                 ifelse(hsa.mir.592 >=cut625,
                                                        "low", "high"))
t_VU_rc_vst_surv$CMS <- sele_pred_rc_vst_VU_RF[match(t_VU_rc_vst_surv$Sample.ID,
                                               names(sele_pred_rc_vst_VU_RF))]
t_VU_rc_vst_surv <- t_VU_rc_vst_surv[-220,]
t_VU_rc_vst_surv$bin.sel.mir <- factor(t_VU_rc_vst_surv$bin.sel.mir)
#}

# install.packages('survival')
# install.packages('survminer')
# BiocManager::install("RTCGA.clinical") # data for examples
library(survival)
library(survminer)
library(RTCGA.clinical)

fit <- survfit(Surv(OS, Event_OS) ~ bin.sel.mir,
               data = t_VU_rc_vst_surv)
# visualize with survminer
ggsurvplot(fit, data = t_VU_rc_vst_surv, risk.table = TRUE, pval = T,
           palette = c("#2a7886","#79bac1" ), main="mir592" ) 

### get TCGA survivall data 
TCGA.survInfo <- survivalTCGA(COAD.clinical, extract.cols = "admin.disease_code") 
rc_vst_surv <- cbind( rc_vst, TCGA.survInfo[
  match(sub("-01.*","",rownames(rc_vst)), 
  TCGA.survInfo$bcr_patient_barcode),])
rc_vst_surv <- rc_vst_surv %>% dplyr::mutate(bin.sel.mir = 
                                         ifelse(hsa.mir.200a >=
                                                  median(rc_vst_surv$hsa.mir.200a),
                                                        "low", "high"))
rc_vst_surv$bin.sel.mir <- factor(rc_vst_surv$bin.sel.mir)

fit <- survfit(Surv(times, patient.vital_status) ~ bin.sel.mir,
               data = rc_vst_surv)
# visualize with survminer
ggsurvplot(fit, data = rc_vst_surv, risk.table = TRUE, pval = T,
           palette =  c("#2a7886","#79bac1" ), main="out-of-box" ) 

# Fit a Cox proportional hazards model
surv_object <- Surv(time = rc_vst_surv$times, 
                    event = rc_vst_surv$patient.vital_status)
fit.coxph <- coxph(surv_object ~ hsa.mir.625 +
                   hsa.mir.375 + hsa.mir.592+ 
                   hsa.mir.552 + hsa.mir.141 + hsa.mir.200a,  
                   data = rc_vst_surv)
ggforest(fit.coxph, data = rc_vst_surv)


# Check class proportions
t_VU_rc_vst_surv %>% ggplot(aes(x = CMS, fill = CMS)) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=c(#"#DAA520", 
                      "#1874CD", 
                      "#CD2990", 
                      "#458B74")) +
  ggtitle("CMS counts")



##plot density distribution
ggplot( rc_vst, aes( x=hsa.mir.4636 ) ) + 
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
