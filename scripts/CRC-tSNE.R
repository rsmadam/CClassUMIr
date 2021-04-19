##### general functions and paramters #####
# Libraries needed
library(tidyverse)
library(readxl)
library("DESeq2")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(TCGAbiolinks)
library(CMSclassifier)
library(ggsci)
library(ggrepel)
library(wesanderson)
library(corrplot)
library(Rtsne)

plotDir <- "/Users/ronjaadam/projects/miRNA_mCRC/CMS-miRaCl/analyses/plots/"
projDir <- "/Users/ronjaadam/projects/miRNA_mCRC/CMS-miRaCl/"
projHome <- "/Users/ronjaadam/projects/miRNA_mCRC/"


########## CMS labels publication #######
##Get CMS labels from Louis table
CMS_samples <- read_excel('Data/raw/CMS4 classification.xls', na = "NA",
                            col_types = c("text", rep("numeric", 4),
                                          "text", "numeric"))
# Rename first column
colnames(CMS_samples)[1] <- "SampleId"
# Only use the data from TCGA
CMS_samples <- CMS_samples[ grep("^tcga", CMS_samples$SampleId) , ]
# Adjust the SampleId so it matches with the miRNA df
CMS_samples$SampleId <- gsub("tcga_rnaseq.", "", CMS_samples$SampleId)
# Filter the ones that have a minimal p_value <0.5: not really needed here, they are censored to NA anyways
# Remove duplicates
CMS_samples$SampleId.2 <- CMS_samples$SampleId
CMS_samples$SampleId <- gsub("ll.", "", CMS_samples$SampleId)
CMS_samples <- CMS_samples[!duplicated(CMS_samples$SampleId),]
# from the paper link https://www.synapse.org/#!Synapse:syn2623706/files/
CMS_samples_official <- read.table(file = "Data/cms_labels_public_all.txt")
CMS_samples$CMS_final_netw_RF <- CMS_samples_official$V5[match(CMS_samples$SampleId, 
                                                               as.character(CMS_samples_official$V1))]
remove(CMS_samples_official)

TCGARnaseqDF.COAD <- read.csv(paste0(projDir,"/Data/TCGA-COAD-mRNA-GA-HiSeq_scaled_estimate_01Aonly_BR.csv"),
                                   row.names = 1)
predictedCMS.COAD <- read.csv(paste0(projDir, "/Data/TCGA-COAD-allRNA-01Aonly_CMS-labels.csv"),
                    row.names = 1)


################### get COAD miRNA from TCGA ################

## fetch TCGA COAD miRNAs:
query.exp <- GDCquery(project = c("TCGA-COAD", "TCGA-READ"),
                      legacy = F,
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification",
                      experimental.strategy = "miRNA-Seq"
)
GDCdownload(query.exp)
TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
TCGAmiR.CRC <- TCGARnaseqSE[,grep("read_count_", colnames(TCGARnaseqSE)) ] #or normalized_count_ for normalized_results
colnames(TCGAmiR.CRC) <- sub("read_count_","",colnames(TCGAmiR.CRC))
rownames(TCGAmiR.CRC) <- TCGARnaseqSE$miRNA_ID
remove(TCGARnaseqSE)
##TCGA data has similar counts for isomirs, summarize by mean up
TCGAmiR.CRC$miRcomb <- sub("-[1-9]$","",row.names(TCGAmiR.CRC))
miR_CRC_comb <- as.data.frame(TCGAmiR.CRC %>% dplyr::group_by(miRcomb) %>%
                                 dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                     .funs = mean)) #affects around 180miRNA
row.names(miR_CRC_comb) <- miR_CRC_comb$miRcomb #afterwards remove column 1
miR_CRC_comb$miRcomb <- NULL
##keep only one sample per patient
colnames(miR_CRC_comb) <- sub("-01A-.*","",colnames(miR_CRC_comb))
miR_CRC_comb <- miR_CRC_comb[,-grep("-.+-.+-.+-.+", colnames(miR_CRC_comb) )] #get rid of 20 normal and extra samples
miR_CRC_comb <- miR_CRC_comb[,!duplicated(colnames(miR_CRC_comb))] #miR_CRC_comb more duplicated patients
# ##normalize, in right direction for RF
# miR_CRC_vst <- as.data.frame( t( varianceStabilizingTransformation(
#   as.matrix( round( miR_CRC_comb, 0 ) ) ) ) )
# colnames(miR_CRC_vst) <- gsub("-", ".", colnames(miR_CRC_vst))
# 


########## CMS labels ############
####combine miR data with labels
##get the CMS labels from CMSclassifier
rownames(miR_CRC_vst) <- sub("//.*","",rownames(miR_CRC_vst))
miR_COAD_vst$CMS.cl.rf <- predictedCMS.COAD[match(rownames(miR_COAD_vst),
                                                        rownames(predictedCMS.COAD)),"RF"]
miR_COAD_vst$CMS.cl.rf <- factor(miR_COAD_vst$CMS.cl.rf)

##get the CMS labels Guinney et al. publication
miR_CRC_vst$CMS.lv <- CMS_samples$CMS[match(rownames(miR_CRC_vst),
                                             CMS_samples$SampleId)]#get the CMS labels from Louis' table
miR_COAD_vst$CMS.lv <- factor(miR_COAD_vst$CMS.lv)
miR_COAD_vst$CMS.guin <- CMS_samples$CMS_final_netw_RF[match(rownames(miR_COAD_vst),
                                             CMS_samples$SampleId)]#get the CMS labels from Louis' table
miR_COAD_vst$CMS.guin[which(miR_COAD_vst$CMS.guin == "NOLBL")] <- NA
miR_COAD_vst$CMS.guin <- droplevels(miR_COAD_vst$CMS.guin)
miR_COAD_vst$CMS <- miR_COAD_vst$CMS.cl.rf
miR_COAD_vst$CMS[-which(miR_COAD_vst$CMS.cl.rf==
                          miR_COAD_vst$CMS.guin & miR_COAD_vst$CMS.guin==
                          miR_COAD_vst$CMS.lv)] <- NA #keep only samples where the predictions were the same
write.table(miR_COAD_vst, paste0(projHome, "miRNA classifier/Data/miR_COAD_vst_CMS.txt"))

miR_COAD_vst <- read.table( "Data/miR_COAD_vst_CMS.txt") #here the outliers and non-classifiable are still in but the normal/extra samples are out

### batch effect removal in miR data (AA batch is different)
#limma was suggested also by Michael Love, https://support.bioconductor.org/p/76099/
tssAnnot <- data.frame("CMS"=miR_COAD_vst$CMS.cl.rf) #use the predicted labels here in order to have labels for all
tssAnnot$CMS <- factor(tssAnnot$CMS, levels=c("CMS1", "CMS2", "CMS3", "CMS4", "NOLBL"))
tssAnnot$CMS[which(is.na(tssAnnot$CMS))] <- "NOLBL"#use the predicted labels here in order to have labels for all

## pre Batch effect removal pre outlier removal
rownames(tssAnnot) <- rownames(miR_COAD_vst)
tssAnnot$TSS <- sub("-.*","", sub( "TCGA-", "", rownames(tssAnnot) ) )
tssAnnot$TSS2 <- factor(tssAnnot$TSS, labels=c(rep("other",5),"AA",
                                            rep("other", 19 )))#only AA vs other
designM <- model.matrix( ~ 0+ tssAnnot$CMS ) ## do i need an intercept? otherwise it compares all to CMS1 or CMS2 if releveled?
coad_vst_limBatchRem <- t(limma::removeBatchEffect(t(miR_COAD_vst[,grep("hsa",
                                                                colnames(miR_COAD_vst))]), #the entire dataset, 339 samples
                                                 design = designM,
                                                 batch=tssAnnot$TSS) )
#function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects
coad_vst_limBatchRem <- as.data.frame(coad_vst_limBatchRem)
coad_vst_limBatchRem$CMS <- miR_COAD_vst$CMS




miR_READ_vst <- read.csv( paste0(projHome, "CMS-miRaCl/Data/TCGA-miR_READ_vst_CMS-labels.csv"),
                          row.names=1)
### batch effect removal in miR data (AA batch is different)
#limma was suggested also by Michael Love, https://support.bioconductor.org/p/76099/
tssAnnot <- data.frame("CMS"=miR_READ_vst$CMS.cl.rf) #use the predicted labels here in order to have labels for all
tssAnnot$CMS <- factor(tssAnnot$CMS, levels=c("CMS1", "CMS2", "CMS3", "CMS4", "NOLBL"))
tssAnnot$CMS[which(is.na(tssAnnot$CMS))] <- "NOLBL"#use the predicted labels here in order to have labels for all

## pre Batch effect removal pre outlier removal
rownames(tssAnnot) <- rownames(miR_READ_vst)
tssAnnot$TSS <- sub("-.*","", sub( "TCGA-", "", rownames(tssAnnot) ) )
designM <- model.matrix( ~ 0+ tssAnnot$CMS ) ## do i need an intercept? otherwise it compares all to CMS1 or CMS2 if releveled?
read_vst_limBatchRem <- t(limma::removeBatchEffect(t(miR_READ_vst[,grep("hsa",
                                                                      colnames(miR_READ_vst))]), #the entire dataset, 339 samples
                                                 design = designM,
                                                 batch=tssAnnot$TSS) )
#function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects
read_vst_limBatchRem <- as.data.frame(read_vst_limBatchRem)
read_vst_limBatchRem$CMS <- miR_READ_vst$CMS



##the two datasets together
miR_CRC_vst <- rbind(miR_COAD_vst[,1:1706], 
                     miR_READ_vst[,1:1706])
crc_vst_limBatchRem <-  rbind(coad_vst_limBatchRem, read_vst_limBatchRem)
miR_CRC_vst$CMS <- c(miR_COAD_vst$CMS.cl.rf, miR_READ_vst$CMS.cl.rf)

miR_CRC_vst$CMS  <- factor(miR_CRC_vst$CMS , labels=paste0("CMS", 1:4))
miR_CRC_vst$CMS <- factor(miR_CRC_vst$CMS , levels=c(levels(miR_CRC_vst$CMS ), "NOLBL"))
miR_CRC_vst$CMS[is.na( miR_CRC_vst$CMS )] <- "NOLBL"
miR_CRC_CMS <- miR_CRC_vst$CMS
names(miR_CRC_CMS) <- rownames(miR_CRC_vst)

colCRC <- factor( sub("-.*","",sub("TCGA-","",rownames(miR_CRC_vst))), 
                  levels = c(levels( factor( sub("-.*","",sub("TCGA-","",rownames(miR_COAD_vst))) ) ), 
                  levels( factor( sub("-.*","",sub("TCGA-","",rownames(miR_READ_vst))) ) )), 
                  ordered=T)
shOrgan <- ifelse(colCRC %in% levels( factor( sub("-.*","",sub("TCGA-","",rownames(miR_COAD_vst))) ) ), 
                  "colon", "rectum")
# 
# designM <- model.matrix( ~ 0 + miR_CRC_CMS ) ## do i need an intercept? otherwise it compares all to CMS1 or CMS2 if releveled?
# 
# CRC_vst_limBatchRem <- t(limma::removeBatchEffect(t(miR_CRC_vst[,grep("hsa",
#                                                                       colnames(miR_CRC_vst))]), #the entire dataset, 339 samples
#                                                  design = designM,
#                                                  batch=colCRC) )



### tSNE raw:
set.seed(5678)
tsne_model_COAD <- Rtsne( as.matrix( t(miR_CRC_comb[grep("hsa", rownames(miR_CRC_comb)),
                                                 match(rownames(miR_CRC_vst), colnames(miR_CRC_comb))])), 
                          check_duplicates=F, 
                          pca=TRUE, perplexity=20, 
                          theta=0.5, dims=2)

ggplot(as.data.frame(tsne_model_COAD$Y), 
       aes(x=V1, y=V2, color = colCRC, shape=shOrgan))+
  geom_point() +
  xlab("") + ylab("") +
  theme_minimal()+
  scale_color_manual(values =c( c(wes_palette("Darjeeling1", 25, 
                                              type = c( "continuous")) ),#[1:25], 
                                c(wes_palette("Moonrise3", 26,  type = c( "continuous"))[1:13] )))

ggsave(paste0(plotDir,"/COAD-READ-tSNE_raw.pdf" ), useDingbats=F )


### tSNE pre and post BR:
set.seed(5678)
tsne_model_COAD <- Rtsne( as.matrix( miR_CRC_vst[,
                                                   grep("hsa", colnames(miR_CRC_vst))]), 
                          check_duplicates=F, 
                          pca=TRUE, perplexity=20, 
                          theta=0.5, dims=2)

ggplot(as.data.frame(tsne_model_COAD$Y), 
       aes(x=V1, y=V2, color = colCRC, shape=shOrgan))+
  geom_point() +
  xlab("") + ylab("") +
  theme_minimal()+
  scale_color_manual(values =c( c(wes_palette("Darjeeling1", 25, 
                                            type = c( "continuous")) ),#[1:25], 
                     c(wes_palette("Moonrise3", 26,  type = c( "continuous"))[1:13] )))

ggsave(paste0(plotDir,"/COAD-READ-tSNE_rc_vst_pre-limBR.pdf" ), useDingbats=F )



set.seed(5678)
tsne_model_COAD <- Rtsne( as.matrix( crc_vst_limBatchRem[,
                                                 grep("hsa", colnames(crc_vst_limBatchRem))]), 
                          check_duplicates=F, 
                          pca=TRUE, perplexity=20, 
                          theta=0.5, dims=2)

ggplot(as.data.frame(tsne_model_COAD$Y), 
       aes(x=V1, y=V2, color = colCRC, shape=shOrgan))+
  geom_point() +
  xlab("") + ylab("") +
  theme_minimal()+
  scale_color_manual(values =c( c(wes_palette("Darjeeling1", 25, 
                                              type = c( "continuous")) ),#[1:25], 
                                c(wes_palette("Moonrise3", 26,  type = c( "continuous"))[1:13] )))

ggsave(paste0(plotDir,"/COAD-READ-tSNE_rc_vst_post-limBR.pdf" ), useDingbats=F )

