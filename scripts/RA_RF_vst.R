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
library(wesanderson)
library(ggsci)
#install.packages('e1071', dependencies=TRUE)
set.seed(67)
plotDir <- "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/plots/"

##################################################################################################################################      
# Open the csv file
rpm <- read.table("data/cleaned/RPM_prob_p005.csv")
rc <- read.table("data/cleaned/RC_prob_p005.csv")
rc$CMS <- relevel(rc$CMS, "CMS2")

# Only use the features with raw read count mean > 1
feature_names <- names(which( colMeans(rc[,-grep("CMS",colnames(rc))])>1))
rpm <- rpm[,c(feature_names,"CMS")]
rpm$Sample.ID <- row.names(rpm)

## keep out of bag samples for performance testing 
set.seed(42)
out_Test <- sample(1:274, 274/5)
# out_Test <- c(251,208,135,239,112,194,277,272,117,316,
#               34,310,298,13,266,125,249,96,196,157,
#               197,98,178,287,144,171,244,164,88,14,226,
#               60,324,74,274,66,44,199,182,276,
#               73,129,295,195,156,219,92,314,106,230,
#               110,271,315,27,228,306,111,297,321,114,
#               52,8,283,202,267,339,63)
#rpm_out <- rpm[out_Test, ] #out_Test:
rpm <- rpm[out_Test,"CMS"]


# save(selectedRF,file = "selRF.RData")
# save(model_RF,file = "lastModelRF.RData")
# save(rpm_out, file= "rpm_out.RData")
# save(rpm, file="rpm_input.RData")
load(file = "selRF.RData")
load(file = "lastModelRF.RData")
load(file= "~/projects/miRNA_mCRC/miRNA classifier/rpm_out.RData")
load( file="~/projects/miRNA_mCRC/miRNA classifier/rpm_input.RData")
TCGARnaseqRaw <- read.table("~/projects/miRNA_mCRC/TCGA_COADREAD_mRNA_raw_gene_counts.tsv")


### normalize raw counts (rc) with varianceStabilizingTransformation 
library(DESeq2)
rc_CMS <- rc$CMS
rc_vst <- varianceStabilizingTransformation( t( round( rc[, 
                                                          feature_names],0))) #toss low count features
rc_vst <- as.data.frame( t( rc_vst[,] ) )

###########
### use new download and concordant classifier labels
rc_CMS <- droplevels(miR_COAD_vst$CMS[grep("CMS", miR_COAD_vst$CMS)])
#rc_CMS <- relevel(rc_CMS, "CMS2")
feature_names <- gsub("-",".",names(which(rowMeans(
  miR_COAD_comb[,grep("CMS", miR_COAD_vst$CMS)])>1))) #expressed before vst?
rc_vst <- miR_COAD_vst[grep("CMS", miR_COAD_vst$CMS),feature_names]



#create the TSS annotation for use in some plots: 
# pre Batch effect removal pre outlier removal
dfAnnot <- data.frame("CMS"=rc_CMS)
rownames(dfAnnot) <- rownames(rc_vst)
dfAnnot$TSS <- sub("-.*","", sub( "TCGA-", "", rownames(dfAnnot) ) )
dfAnnot$TSS <- factor(dfAnnot$TSS, labels=c(rep("other",1),"AA",
                                            rep("other", 14 )))
dfAnnot$mean <- rowMeans(rc_vst[,grep("hsa", colnames(rc_vst))])
dfAnnot$sums <- rowSums(rc_vst[,grep("hsa", colnames(rc_vst))])
dfAnnot$sDs <- rowSds(as.matrix(rc_vst[,grep("hsa", colnames(rc_vst))]))
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
} #Tukey's mild outlier definition

dfAnnot$outls <- ifelse(is_outlier(dfAnnot$sDs), sub("-.*", "",
                                                     sub("TCGA-..-","",rownames(dfAnnot))), 
                        as.numeric(NA))
##outliers with low mean AND high SD (across all samples):
outlierS <- rownames(dfAnnot[is_outlier(dfAnnot$sDs),])
              #        ("TCGA-AA-3555-01A", "TCGA-AA-3560-01A",
              # "TCGA-AA-3664-01A", "TCGA-AA-A00Q-01A")

##batch effect removal (AA batch is different)
#limma was suggested also by Michael Love, https://support.bioconductor.org/p/76099/
##TODO: test to remove batch effects with ComBat
#rc_CMS <- relevel(rc_CMS, "CMS2")
rc_vst$CMS <- rc_CMS
designM <- model.matrix( ~ rc_CMS ) ## do i need an intercept? does this compare all to CMS1? 
rc_vst_limBatchRem <- t(limma::removeBatchEffect(t(rc_vst[,grep("hsa", 
                                                            colnames(rc_vst))]), #the entire dataset, 339 samples
                                                 design = designM,
                                                 batch=dfAnnot$TSS) )#only AA vs other
#function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects
rc_vst_limBatchRem <- as.data.frame(rc_vst_limBatchRem)
rc_vst_limBatchRem$CMS <- rc_CMS

## meanSD Plot should show horizontal line
meanSdPlot(t(as.matrix(rc_vst_limBatchRem[,grep("hsa", 
                                    colnames(rc_vst_limBatchRem))])))
ggsave(paste0(plotDir, "rc_vst_BR_AAvsOther_meanRankSDPlot.pdf"))

#toss outliers and the out of training samples
rc_vst_BR <- rc_vst_limBatchRem[-c(out_Test, 
                                   which(rownames(rc_vst_limBatchRem) %in%
                                                  outlierS)),]  
rc_vst_BR$Sample.ID <- row.names(rc_vst_BR)


#create the TSS annotation for use in some plots after BR
dfAnnot_BR <- data.frame("CMS"=rc_vst_BR$CMS)
rownames(dfAnnot_BR) <- rownames(rc_vst_BR)
dfAnnot_BR$TSS <- sub("-.*","", sub( "TCGA-", "", rownames(dfAnnot_BR) ) )
# dfAnnot_BR$TSS <- factor(dfAnnot_BR$TSS, labels=c("other","AA",
#                                             rep("other", 14 )))
dfAnnot_BR$mean <- rowMeans(rc_vst_BR[,grep("hsa", 
                                            colnames(rc_vst_limBatchRem))])
dfAnnot_BR$sums <- rowSums(rc_vst_BR[,grep("hsa", 
                                           colnames(rc_vst_limBatchRem))])
dfAnnot_BR$sDs <- rowSds(as.matrix(rc_vst_BR[,grep("hsa", 
                                                   colnames(rc_vst_limBatchRem))]))

## visualize batch effects based on boxplots of samples or batches:
library(ggrepel)
ggplot(dfAnnot_BR, aes(x=TSS, y=sDs, colour=TSS)) + 
  geom_boxplot(outlier.size=0.8) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.6) +
  ggtitle("SDs of rc_vst_BR") +
  theme_classic() +
  #geom_text_repel(aes(label = outls), na.rm = TRUE) +
  scale_color_manual(values = wes_palette("Darjeeling1", 22,
                                         type = c( "continuous")) ) ## used this for the TSS
ggsave(paste0(plotDir, "rc_vst_BR_designM&AAvsother_boxplotSDs_batchTSS.pdf"))


##make the out of training data
rc_vst_out <- as.data.frame(rc_vst_limBatchRem[out_Test, ])
rc_vst_out$Sample.ID <- row.names(rc_vst_out)
rc_vst_out$CMS <- rc_vst[row.names(rc_vst_out), "CMS"]
rc_vst_out <- rc_vst_out[-which(rownames(rc_vst_out) %in% outlierS),] #make sure no outliers in validation set




##################################################################################################################################      
####################### RF training with Cross validation ####################### 

##### RF with caret ######
##preprocessing: visualize
library(RColorBrewer)
##identifying low info features, plot their correlation
nzv <- nearZeroVar(rc_vst_BR, saveMetrics=TRUE)
summary(nzv$nzv) #no features to remove
descrCor <-  cor(rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))])
library(corrplot)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75, exact=T)
corrplot(descrCor[highlyCorDescr,highlyCorDescr], type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"), addgrid.col = NA, tl.cex = .5,
         tl.col = "#999999")
fairlyCorDescr <- findCorrelation(descrCor, cutoff = .25, exact=T)

summary(descrCor[upper.tri(descrCor)])
filteredDescr <- rc_vst_BR[,-c(highlyCorDescr, grep("S", colnames(rc_vst_BR)))]

## visualization of class distances
centroids <- classDist( filteredDescr, rc_vst_BR$CMS, pca = T, keep = 30 )
distances <- predict(centroids, rc_vst_BR[,1:491])
distances <- as.data.frame(distances)
head(distances)
pdf(paste0(plotDir, "vst_miRNA_limBRem_classCentroidDist_XY_CMS.pdf"),
    onefile = T)
xyplot(dist.CMS1 ~ dist.CMS3,
       data = distances, 
       groups = rc_vst_BR$CMS, 
       auto.key = list(columns = 2, col = paletteCMS ), 
       col = paletteCMS)
xyplot(dist.CMS2 ~ dist.CMS4,
       data = distances, 
       groups = rc_vst_BR$CMS, 
       auto.key = list(columns = 2, col = paletteCMS ), 
       col = paletteCMS)
dev.off()


### PCA plot ###
#install.packages("factoextra")
library(factoextra)
prin_comp <- prcomp(rc_vst_BR[,-c(highlyCorDescr,492,493)], scale. = T)
pr_var <- round( prin_comp$sdev^2, 5 )
prop_varex <- round( pr_var/sum(pr_var), 5)
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
var.coord <- t(apply(prin_comp$rotation, 1, var_coord_func, prin_comp$sdev)) 
head(var.coord[order(var.coord[,"PC1"]), "PC1"]) #the top features
tail(var.coord[order(var.coord[,"PC1"]), "PC1"])


library(caret)
trans <- preProcess(rc_vst_BR[ #-grep("AA", rownames(rc_vst_BR))
                        ,grep("hsa",colnames(rc_vst_BR))], 
                    method=c(#"BoxCox",
                             "center","scale", 
                             "pca"),
                    thresh = list(thresh = 0.50))
head(trans$rotation[order(trans$rotation[,"PC1"]), "PC1"]) #the top neg features
tail(trans$rotation[order(trans$rotation[,"PC1"]), "PC1"]) #the top features
featPC134 <- unique( names( c( head(trans$rotation[order( 
  abs( trans$rotation[,"PC3"] ) ), "PC3"], 20), #the top abs features
  head(trans$rotation[order( 
    abs( trans$rotation[,"PC4"] ) ), "PC4"], 20), 
  head(trans$rotation[order( 
    abs( trans$rotation[,"PC1"] ) ), "PC1"], 20) ) ) )

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
pdf(paste0(plotDir, "vst_miRNA_limBRem_PCA.pdf"),
    onefile = T)
fviz_eig(prin_comp, type="lines" ) #plot variances
ggplot(PC,aes(x=PC1,y=PC2, 
              colour= rc_vst_BR$CMS[])) +
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

#### plot kmeans class labels on PCA #####
kM <- kmeans(rc_vst_BR[, 1:491], 4, nstart = 20)
kM$cluster <- as.factor(kM $cluster)
summary(kM$cluster)
ggplot(rc_vst_BR, aes(PC$PC1, 
                   PC$PC2, 
                   color = kM$cluster)) + geom_point()
#-> CMS2 and 3 are hard to separate
dfAnnot_BR$kM3 <- kM$cluster
dfAnnot_BR$kM3CMS <- factor(dfAnnot_BR$kM3, 
                            labels = c("CMS1", "CMS4", "CMS23"))
dfAnnot_BR$kM3CMS <- relevel(dfAnnot_BR$kM3CMS, "CMS23")
dfAnnot_BR$CMS23 <- factor(dfAnnot_BR$CMS, 
                          labels= c("CMS23", "CMS1", "CMS23", "CMS4"))
rc_vst_BR$CMS23 <- dfAnnot_BR$CMS23
summary(dfAnnot_BR$CMS23 == dfAnnot_BR$kM3CMS)

#### tSNE plot ####
set.seed(9)  
library(Rtsne)
library(ggplot2)
library(ggsci)
tsne_model_1 <- Rtsne(as.matrix(rc_vst_BR[,grep("hsa",colnames(rc_vst_BR))]), 
                      check_duplicates=FALSE, 
                     pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_1 <- as.data.frame(tsne_model_1$Y) 
ggplot(d_tsne_1, aes(x=V1, y=V2, colour= dfAnnot_BR$TSS))+#rc_vst_BR$CMS)) + #colour=kM$cluster))+# 
  geom_point() +
  xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme_minimal()+
  # scale_color_npg()+
  scale_color_manual(values = wes_palette("Darjeeling1", 22,
                                          type = c( "continuous")) ) ## used this for the TSS
  scale_colour_manual(values = paletteCMS)
ggsave(paste0(plotDir, "vst_miRNA_limBRem_tSNE_CMS.pdf"))

## from the kmeans it looks like there's a good concordance between 
# CMS1 and kCl2, CMS2 and kCl1, CMS4 and kCl3, 
# whereas CMS3 and kCl4 is more outliers and intermediate samples
# leave out CMS3 samples for everything? or combined class CMS2+CMS3? No
#TODO: check if using only concordant samples in the training helps?? CMS have different probabilities 

library(circlize)

#### heatmap VU mir Data ####
dfAnnot_BR <- data.frame("CMS"=rc_vst_BR$CMS)
rownames(dfAnnot_BR) <- rownames(rc_vst_BR)
dfAnnot_BR$pval <- as.numeric(CMS_samples$min_Pval[match(sub("-01A","",
                                                             rownames(dfAnnot_BR)), 
                                                         CMS_samples$SampleId)])
col_fun = colorRamp2(c(0, 0.05), c("darkblue", "white"))
ha <- HeatmapAnnotation( df = dfAnnot_BR[,c("pval", "CMS")],

#### heatmap TCGA mir data #### or for featPC134
dfAnnot_BR$pval <- as.numeric(CMS_samples$min_Pval[match(sub("-01A","",
                                                             rownames(dfAnnot_BR)), 
                                              CMS_samples$SampleId)])
col_fun = colorRamp2(c(0, 0.05), c("darkblue", "white"))
ha <- HeatmapAnnotation( df = dfAnnot_BR[,c("pval", "CMS")],
                         col = list("pval" = col_fun,
                                     "CMS" = c("CMS1"="#E79E1B", "CMS2"="#0071B1",
                                                  "CMS3"="#C45597", "CMS4"="#009C74",
                                                  "normal"="#410071", "NA"="#999999")
                             ),
                         na_col = "grey")

imp.miR <- rownames(varImp(model_RF)$importance)[ which(varImp(model_RF)$importance>
                                                          mean(varImp(model_RF)$importance) )] #>mean
impDF <- data.frame("imp"=varImp(model_RF)$importance[imp.miR,])
rownames(impDF) <- imp.miR
colRdYlBl <- circlize::colorRamp2(c(seq(-5, -2, length = 2),seq(-2,2 ,length = 7),
                          seq(2, 5, length = 2)), rev( c('#a50026','#d73027','#f46d43','#fdae61','#fee090', '#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695') ) )

haMRNA <- HeatmapAnnotation( df = data.frame("CMS"=CMS_samples$CMS[match(colnames(TCGARnaseqDF),
                                                    gsub("-",".",CMS_samples$SampleId))],
                                             "pval"=as.numeric(CMS_samples$min_Pval[match(colnames(TCGARnaseqDF),
                                                                     gsub("-",".",CMS_samples$SampleId))])),
                         col = list("pval" = col_fun,
                                    "CMS" = c("CMS1"="#E79E1B", "CMS2"="#0071B1",
                                              "CMS3"="#C45597", "CMS4"="#009C74",
                                              "normal"="#410071", "NA"="#999999")
                         ),
                         na_col = "grey")
mat <- t( scale(as.matrix(t( #mRNA:
  na.exclude(TCGARnaseqDF[match(CMSsymbols,
                     rownames(TCGARnaseqDF)),]))),#rc_vst_BR[,imp.miR]), 
                center=T, scale=T ) )
htFirst <- Heatmap( mat,
                    col = colRdYlBl,
                    cluster_columns = T, cluster_rows = T, 
                    #clustering_distance_columns = "ward.D",#
                    #clustering_distance_rows = "ward.D",
                    top_annotation = haMRNA, 
                    #show_row_names = T, 
                    column_split = CMS_samples$CMS[match(colnames(TCGARnaseqDF),
                                                         gsub("-",".",CMS_samples$SampleId))],#dfAnnot_BR[, "CMS"],
                    show_column_names = T, 
                    column_names_gp = gpar(fontsize = 4),
                    row_names_gp = gpar(fontsize = 4),
                    heatmap_height = unit(1.2, "mm")*nrow(mat),
                    heatmap_width = unit(1.2, "mm")*ncol(mat))
htSecond <- Heatmap( as.matrix(log2(impDF)) ,
                    col =  colorRamp2(c(10, 0), c("darkred", "white")),
                    show_row_names = T,
                    width = unit(3, "mm"),
                    column_title_side ="bottom",
                    column_title = "log2(Impurity)")

print(htFirst)
pdf( file = paste0(  "analyses/plots/rc_vst_heatmap_mRNA.pdf"),
     onefile = TRUE, height=20, width=20 )
print(htFirst)
dev.off() #



#TODO: test everything with mRNA


#### plot cluster dendrogram ####
hc <- hclust(dist(rc_vst_BR[,1:491]), "ward.D") #TODO: test other distance methods

# customized plot
# make color code
colrc_vst.CMS <- factor(rc_vst_BR$CMS23, labels = c(
                                         "blue4","darkgoldenrod1", 
                                         #"deeppink3",
                                         "mediumseagreen"))
# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- as.character(colrc_vst.CMS[which(rownames(rc_vst_BR) == a$label)])
    attr(n, "nodePar") <- c(a$nodePar,
                            pch="o",
                            col = labCol) }
  return(n) }
hcd <- as.dendrogram(hc)
# using dendrapply to get color-annotation for each node
clusDendro = dendrapply(hcd, colLab)
plot(clusDendro)


#### testing accuracy of RF vs SVM ####
library(caret)
##create train vs test splits beforehand
# Create custom indices: myFolds
rc_vst_BR <- rc_vst_BR[-which(rc_vst_BR$CMS %in% "NOLBL"),]
rc_vst_BR$CMS <- droplevels(rc_vst_BR$CMS)
myFolds <- createFolds(rc_vst_BR$CMS, k = 10)

# Create reusable trainControl object: myControl
myControl <- trainControl(
  #summaryFunction = twoClassSummary,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = TRUE,
  savePredictions = TRUE,
  index = myFolds
)
tGrid <- data.frame("mtry" = c(15,50,100,150,200,250,300),
                       "splitrule" = "gini",
                       "min.node.size" = 5 )

#TODO use metric "ROC" for model comparison,
#TODO: keep nzv but use preprocess pca instead
#TODO: caret::train on whole dataset because it uses the trainControl function as training recipe and to compute stats on cv out of trainin
#TODO: try bootstrap aggregation instead of cv
#TODO: compare different models in model_list using e.g. summary() or dotplot(resamples(model_list), metric="ROC)


# SVM
fit.svm <- train(CMS~., data=rc_vst_BR[,-c(highlyCorDescr,493)], 
                 method="svmRadial", metric="Kappa", trControl=myControl)
# SVM
fit.svm.acc <- train(CMS~., data=rc_vst_BR[,-c(highlyCorDescr,493)], 
                     method="ranger", metric="ROC", trControl=myControl)
# Random Forest
fit.rf <- train(CMS~., data=rc_vst_BR[,-c(highlyCorDescr,493)], 
                method="ranger", metric="Kappa", 
                tuneGrid= tGrid,
                trControl=myControl)
# Random Forest
fit.rf.acc <- train(CMS~., data=rc_vst_BR[,-c(highlyCorDescr,493)], 
                    method="rf", metric="Accuracy", trControl=myControl)
# compare accuracy of models
results <- resamples(list(svm.kappa=fit.svm, svm.acc=fit.svm.acc, 
                          rf.kappa=fit.rf, rf.acc=fit.rf.acc))
summary(results)
pdf("analyses/plots/rc_vst_limBR_CMS23_RFvsSVM_AccuracyKappa.pdf",
    height = 4.1, width = 5.8)
dotplot(results)
dev.off()



###### Classisier with RandomForest or caret ######

# Make an empty variable
cv_comparison = NULL
list_RF <- vector("list", 10)
names(list_RF) <- paste0("model_cv_", 1:10)
rc_vst_BR$CMS == dfAnnot_BR$CMS#NULL # forget about CMS when comaring CMS23
#rc_vst_BR$CMS23 <- dfAnnot_BR$CMS23 # forget about CMS23 when comaring CMS

for (i in 1:10){
  ###added cross validation because due to oversampling we cant assume independence of samples in rc_vst_BR, so oob error rate might be too optimistic
  #https://stats.stackexchange.com/questions/283760/is-cross-validation-unnecessary-for-random-forest
  selectTrain <- rc_vst_BR  %>% group_by(CMS) %>% #or rc_vst_BR[which(rc_vst_BR$CMS != "CMS2"),] #for trying how well separation goes without CMS2
    dplyr::sample_frac(0.7) %>% 
    dplyr::select(Sample.ID)
  rc_vst_BR_train <- rc_vst_BR[which(rc_vst_BR$Sample.ID %in% selectTrain$Sample.ID), ]
  rc_vst_BR_valid <- rc_vst_BR[-which(rc_vst_BR$Sample.ID %in% selectTrain$Sample.ID), ]
  
  ####################### Make the RF model ####################### 
  # ##with randomForest
  model_RF <- randomForest(CMS ~ ., importance = TRUE,
                           replace = F, #less likely to overuse CMS2?
                           proximity = TRUE,
                           data = rc_vst_BR_train[,-c(highlyCorDescr, 492)],
                           ntree = 1000,
                           mtry = 50,
                           keep.inbag=TRUE)#,
                           #cutoff=c(0.2,0.4,0.2,0.2))
                           #classwt=c(0.3,0.9,0.6,0.15))
                           #sampsize=c(30,15,30,30), strata=rc_vst_BR_train$CMS)
  # list_RF[[paste0("model_cv_", i)]] <- model_RF
  # 
  # i <- i+1
  ##with  caret
  model_RF <- caret::train(CMS~., 
                           data=rc_vst_BR[,c(grep("hsa", colnames(rc_vst_BR)), 
                                             grep("CMS", colnames(rc_vst_BR)))],
                           method="ranger", 
                           importance = 'impurity',
                           metric="Kappa", 
                           tuneGrid= tGrid,
                           class.weight=c(0.9,1,0.1,2),
                           trControl= myControl)

  list_RF[[paste0("model_cv_", i)]] <- model_RF
}
#save(model_RF, file="modelRF-caret-ranger-jun2020.R")

##############################################################################################################
# Predict CMS from the out of box data
pred_rc_vst_BR_valid_RF <- predict(model_RF, 
                                   newdata = rc_vst_out)
# Print the confusion matrix for external validation
cmat_RF <- confusionMatrix(pred_rc_vst_BR_valid_RF, rc_vst_out$CMS)
summary(pred_rc_vst_BR_valid_RF)
cmat_RF

## on VU-data?
VU_rc_vst$Sample.ID <-row.names(VU_rc_vst)
VU_valid <- VU_rc_vst[grep("", row.names(VU_rc_vst)),
                      grep("hsa", colnames(VU_rc_vst))]
pred_rc_vst_VU_RF <- predict(model_RF, newdata = VU_valid, type = "prob")
pred_rc_vst_VU_RF$CMS <- predict(model_RF, newdata = VU_valid, type = "raw")
pred_rc_vst_VU_RF$maxP <- rowMax(pred_rc_vst_VU_RF[,1:4])

summary(pred_rc_vst_VU_RF)

## plot VU predictions
data.frame("CMS"=pred_rc_vst_VU_RF[grep("[S,P][0-9]", rownames(VU_valid))]) %>% 
  #[grep("^[P,M][0-9]", t_VU_rc_vst_surv$Sample.ID), "CMS"]) %>% #[!is.na(t_VU_rc_vst_surv$OS), "CMS"]) %>% 
  ggplot(aes(x = CMS, fill = CMS )) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=c(paletteCMS[1:4])) + #"#4F2776"
  ggtitle("CMS counts")

t_VU_rc_vst_surv$CMS <- 
  pred_rc_vst_VU_RF[match(t_VU_rc_vst_surv$Sample.ID,
   rownames(VU_valid))] 
t_VU_rc_vst_surv$sampleOrigin <- 
  clinVU[match(t_VU_rc_vst_surv$Sample.ID,
                             clinVU$sampleID), "SampleOrigin"] 

## plot: what are the tissue origins from the VU predictions?
ggplot(t_VU_rc_vst_surv %>% dplyr::count(CMS, sampleOrigin), 
       aes( CMS, n, fill = sampleOrigin )) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Darjeeling1", 16,
                                         type = c( "continuous")))+#values=c("#4F2776",paletteCMS[c(2,4)])) + #
  ggtitle("CMS counts")
ggsave(paste0(plotDir,"VU_vst_RFkappaCaret_4CMS.pdf"))


#### overlap paired samples P&M
concordPM <- t_VU_rc_vst_surv[grep("^[P,M][0-9]", t_VU_rc_vst_surv$Sample.ID), 
                     c("CMS", "Sample.ID", "sampleOrigin")]
concordPM$patient <- factor(sub("_.*","", sub("\\..*","", 
                                       sub("^[P,M]", "", concordPM$Sample.ID))))
concordPM$PM <- factor(sub("[0-9].*", "",concordPM$Sample.ID), 
                          levels=c("P", "M"), ordered=T)
concordPM$nCMS <- as.numeric(sub("CMS","", concordPM$CMS))
wConcPM <- concordPM %>% spread(PM, CMS)
## plot concordance P&M
ggplot(concordPM, aes(patient, PM, fill= CMS)) + 
  geom_tile()+
  theme_bw()+
  scale_fill_manual(values=paletteCMS[c(2:4)])
ggsave("VU-predict-concordancePM.pdf", width=12, height=3)


##### on READ data 
miR_READ_vst
pred_read_RF <- predict(model_RF, newdata = miR_READ_vst[,1:1706])
summary(pred_read_RF)
confusionMatrix(pred_read_RF, factor(miR_READ_vst$CMS.lv))
#confusionMatrix(pred_read_RF, factor(miR_READ_vst$CMS.cl.rf)) # I don't think thiese predicitons are right, too balanced CMS2,3,4
miR_READ_vst$CMS.cl.rf
## plot: what are the READ predictions
data.frame("CMS"=pred_read_RF[]) %>%#
  ggplot(aes(x = CMS, fill = CMS )) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=c(paletteCMS[])) + #"#4F2776"
  ggtitle("CMS counts")
##plot: what is the read concordance
concREAD <- data.frame("CMS"=factor(c(miR_READ_vst[,"CMS.cl.rf"],
                       as.character(pred_read_RF[]))),
           "pat"=c(rownames(miR_READ_vst[,]),
           rownames(miR_READ_vst[,])),
           "type"= factor(c(rep("mRNA", length(miR_READ_vst[,"CMS.cl.rf"])),
                     rep("miR", length(miR_READ_vst[,"CMS.cl.rf"])))))
ggplot(concREAD, aes(x=pat, y=type, fill= CMS)) + 
  geom_tile()+
  theme_bw()+
  scale_fill_manual(values=paletteCMS[])+
  theme(axis.text=element_text(angle = 90), 
        axis.title= element_text(angle = 90))


#### on independent 3 samples from GSE121842
predict(model_RF, newdata = miRC42_vst)


### store a good model
# selCaretRF <- model_RF
# sele_pred_rc_vst_VU_RF <- pred_rc_vst_VU_RF



###how well did the trees improve the model
oob.error.data <- data.frame(Trees=rep(1:nrow(model_RF$finalModel$err.rate), 
                                        times=4),
                             Type=rep(c("OOB",paste0("CMS", c(23,1,4))),
                                        each=nrow(model_RF$finalModel$err.rate)),
                             Error=c(model_RF$finalModel$err.rate[,"OOB"],
                                     model_RF$finalModel$err.rate[,"CMS23"],
                                     model_RF$finalModel$err.rate[,"CMS1"],
                                     model_RF$finalModel$err.rate[,"CMS4"]))
ggplot(data=oob.error.data,
       aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type)) +
  theme_minimal() +
  scale_colour_manual(values=c(paletteCMS[2],"#4F2776",paletteCMS[4], "#999999"))



# ROC metric 
#library(pROC)

caTools::colAUC(predicted_probabilities, actual, plotROC = TRUE)
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

##for caret R:
N_imp_RF <- importance(model_RF$finalModel)[,"MeanDecreaseGini"] # for caret RF
N_imp_RF <- varImp(model_RF)
N_imp_RF_df <- tibble(variable = names(N_imp_RF), 
                      importance = N_imp_RF ) %>% 
  dplyr::arrange(-importance) 
##for caret ranger:
N_imp_RF <- varImp(model_RF)$importance
N_imp_RF_df <- tibble(variable = rownames(N_imp_RF), 
                    importance = N_imp_RF$Overall ) %>% 
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



### survival effects?? ###
# get suvival data
VU_rc_vst$OS <- clinVU$OS[match(VU_rc_vst$Sample.ID, clinVU$sampleID)]
VU_rc_vst$Event_OS <- clinVU$Event_OS[match(VU_rc_vst$Sample.ID, 
                                              clinVU$sampleID)]

#for(sel.mir in c("hsa.mir.429", "hsa.mir.552", 
#                 "hsa.mir.200a", "hsa.mir.200b", "hsa.mir.141" )){
cut625 <- median(VU_rc_vst[!is.na(VU_rc_vst$OS), "hsa.mir.200b"])
t_VU_rc_vst_surv <- VU_rc_vst %>% dplyr::mutate(bin.sel.mir = 
                                                 ifelse(hsa.mir.200b >=cut625,
                                                        "low", "high"))
t_VU_rc_vst_surv$CMS <- pred_rc_vst_VU_RF[match(t_VU_rc_vst_surv$Sample.ID,
                                                rownames(VU_rc_vst[grep("",
                                                 row.names(VU_rc_vst)), ]))]
#predict doesn't create named vector anymore, so match to same as predicted to

t_VU_rc_vst_surv$bin.sel.mir <- factor(t_VU_rc_vst_surv$bin.sel.mir)
#}

# install.packages('survival')
# install.packages('survminer')
# BiocManager::install("RTCGA.clinical") # data for examples
library(survival)
library(survminer)
library(RTCGA.clinical)

fit <- survfit(Surv(OS, Event_OS) ~ CMS,
               data = t_VU_rc_vst_surv[which(t_VU_rc_vst_surv$CMS %in% c("CMS2", "CMS4") &
                                               grepl("[S][0-9]", #only survival cohort
                                            t_VU_rc_vst_surv$Sample.ID)),],) #[grep("P[0-9]", t_VU_rc_vst_surv$Sample.ID),]
# visualize with survminer
ggsurvplot(fit, data = t_VU_rc_vst_surv[which(t_VU_rc_vst_surv$CMS %in% c("CMS2", "CMS4") &
                                                grepl("[S][0-9]", #only survival cohort
                                                      t_VU_rc_vst_surv$Sample.ID)),], 
           risk.table = TRUE, pval = T, xlim=c(0,2000),
           palette = paletteCMS[c(2,4)], main="CMS" )#"#4F2776" #c("#2a7886","#79bac1" ))# 
# ggsave(paste0(plotDir, "selCaretRFKappa_classweight09100120_surv_KaplanM_VUdata_onlyS.pdf"), 
#        height=8, width=7)
c("#2a7886","#79bac1" )

### get TCGA survival data 
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
t_VU_rc_vst_surv 
rc_vst[-grep("-AA-", rownames(rc_vst)),] %>% ggplot(aes(x = CMS, fill = CMS)) +
  geom_bar() +
  theme_minimal() + 
  scale_fill_manual(values=paletteCMS) +
  ggtitle("CMS counts")
ggsave(paste0(plotDir, "vst_miRNA_noAA_proportions.pdf")


##plot density distribution
ggplot( rc_vst_BR, aes( x=hsa.mir.218 ) ) + 
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
