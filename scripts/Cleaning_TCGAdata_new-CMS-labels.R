# Libraries needed
library(tidyverse)
library(readxl)
library("DESeq2")
library("AnnotationDbi")
library("org.Hs.eg.db")
plotDir <- "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/plots/"
#function needed 
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
} #davetang

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
} #Tukey's mild outlier definition


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
# # Filter the ones that have a minimal p_value <0.5 #not really needed here, they are censored to NA anyways
# CMS_samples <- CMS_samples[which(CMS_samples$min_Pval <= 0.05),]#throw out NAs
# Remove duplicates
CMS_samples$SampleId.2 <- CMS_samples$SampleId
CMS_samples$SampleId <- gsub("ll.", "", CMS_samples$SampleId)
CMS_samples <- CMS_samples[!duplicated(CMS_samples$SampleId),]
# from the paper link https://www.synapse.org/#!Synapse:syn2623706/files/
CMS_samples_official <- read.table(file = "Data/cms_labels_public_all.txt")
CMS_samples$CMS_final_netw_RF <- CMS_samples_official$V5[match(CMS_samples$SampleId, 
                                                               as.character(CMS_samples_official$V1))]
remove(CMS_samples_official)

# #### get raw mRNA from TCGA ####
# ## fetch TCGA data COAD:
library(TCGAbiolinks)
clinical.coad <- GDCquery_clinic(project = c("TCGA-COAD"), type = "clinical")
# CMS_samples$coad <- clinical.coad$tissue_or_organ_of_origin[match(CMS_samples$SampleId,
#                               as.character(clinical$submitter_id))]
# CMS_samples$read <- clinical.coad$tissue_or_organ_of_origin[match(CMS_samples$SampleId,
#                                                 as.character(clinical$submitter_id))]
oldRowNames <- rownames(read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-allRNA-01Aonly-RSEM-ComBat.csv",
                                   header = T)) ##genes
rownames(clinical.coad) <- clinical.coad$submitter_id
# query.exp <- GDCquery(project = c("TCGA-COAD"), 
#                       legacy = T,
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       file.type = "results", # or normalized_results
#                       experimental.strategy = "RNA-Seq", 
#                       platform = "Illumina GA"
# )
# GDCdownload(query.exp)
# TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
# #after packages update, rownames are gone. using the old rownames from a saved table because the number of rows is the same
# TCGARnaseqRSEM.GA <- as.data.frame(TCGARnaseqSE)[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_
# query.exp <- GDCquery(project = c("TCGA-COAD"), 
#                       legacy = T,
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       file.type = "results", # or normalized_results
#                       experimental.strategy = "RNA-Seq", 
#                       platform = "Illumina HiSeq"#OR: GA (for READ there seem to be duplicates with Illumina HiSeq & Illumina GA)
# )
# GDCdownload(query.exp)
# TCGARnaseqSE <- as.data.frame(GDCprepare(query= query.exp, summarizedExperiment = F, save = F))
# TCGARnaseqRSEM.HiSeq <- as.data.frame(TCGARnaseqSE)[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_
# 
# ##for the CMSclassifier use RSEM normalization, log transform and remove batch effect GA/HiSeq
# unique(c(colnames(TCGARnaseqRSEM.HiSeq), colnames(TCGARnaseqRSEM.GA)))
# TCGARnaseqRaw <- cbind( TCGARnaseqRSEM.HiSeq, TCGARnaseqRSEM.GA )
# colnames( TCGARnaseqRaw ) <- c(paste0("HS.",colnames(TCGARnaseqRSEM.HiSeq)), 
#                                paste0("GA.",colnames(TCGARnaseqRSEM.GA)))
# TCGARnaseqRaw <- TCGARnaseqRaw[,grep("-01A-.*", colnames(TCGARnaseqRaw))] #keep only first primary tu sample
# remove(TCGARnaseqSE) #455 unique COAD patients when combining GA and HiSeq
# colnames(TCGARnaseqRaw) <- sub("scaled_estimate_","",colnames(TCGARnaseqRaw))
# colnames(TCGARnaseqRaw) <- sub("-01A-.*", "",colnames(TCGARnaseqRaw))
# batchInfo <- data.frame("sample"=colnames(TCGARnaseqRaw),
#                         "batch"=sub("\\..*","",colnames(TCGARnaseqRaw)))
# ##batch effect removal: ComBat runs long and outputs all NAs :( -> use limma
# # TCGARnaseqDF.COAD <- ComBat(dat=TCGARnaseqDF.COAD, 
# #                             batch=batchInfo$batch,
# #                             par.prior=FALSE, #distribution needs to be nonparametric
# #                             prior.plot=FALSE)
# TCGARnaseqDF.COAD <- limma::removeBatchEffect(as.matrix(log2(TCGARnaseqRaw+0.000001)), #the entire dataset
#                                               batch=batchInfo$batch)
# colnames(TCGARnaseqDF.COAD) <- sub("^...","",colnames(TCGARnaseqDF.COAD)) #GA./HS. being the batch identifier in the colname, keeps HiSeq source in dupliciates
# TCGARnaseqDF.COAD <- TCGARnaseqDF.COAD[,!duplicated(colnames(TCGARnaseqDF.COAD))]
# rownames(TCGARnaseqDF.COAD) <- oldRowNames
# 
# ## get the classification from CMSclassifier package
# library(CMSclassifier)
# rownames(TCGARnaseqDF.COAD) <- sub(".*\\|", "", rownames(TCGARnaseqDF.COAD))#keep only the ENTREZ for classififer
# CMScoad.rf <-classifyCMS(TCGARnaseqDF.COAD, method="RF")
# CMScoad.ss <- classifyCMS(TCGARnaseqDF.COAD, method="SSP")
# CMSsymbols <- mapIds(org.Hs.eg.db, keys=listModelGenes(method = "RF"), 
#                      'SYMBOL', 'ENTREZID')
# write.csv(data.frame("LV.table"=CMS_samples$CMS_final_netw_RF[match(rownames(CMScoad.rf$predictedCMS), 
#                                                                CMS_samples$SampleId)], 
#                 CMScoad.rf$predictedCMS, 
#                 CMScoad.ss$predictedCMS), row.names = T,
#             "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-allRNA-01Aonly_CMS-labels.csv")
# write.csv(TCGARnaseqDF.COAD,
#             "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-mRNA-GA-HiSeq_scaled_estimate_01Aonly_BR.csv")
TCGARnaseqDF.COAD <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-mRNA-GA-HiSeq_scaled_estimate_01Aonly_BR.csv",
                                   row.names = 1)
predictedCMS.COAD <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-allRNA-01Aonly_CMS-labels.csv",
                    row.names = 1)



############
## get COAD miRNA from TCGA
## fetch TCGA COAD miRNAs:
library(TCGAbiolinks)
# query.exp <- GDCquery(project = c("TCGA-COAD"), 
#                       legacy = F,
#                       data.category = "Transcriptome Profiling",
#                       data.type = "miRNA Expression Quantification",
#                       experimental.strategy = "miRNA-Seq"
# )
# GDCdownload(query.exp)
# TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
# TCGAmiR.COAD <- TCGARnaseqSE[,grep("read_count_", colnames(TCGARnaseqSE)) ] #or normalized_count_ for normalized_results
# colnames(TCGAmiR.COAD) <- sub("read_count_","",colnames(TCGAmiR.COAD))
# rownames(TCGAmiR.COAD) <- TCGARnaseqSE$miRNA_ID
# remove(TCGARnaseqSE)
# ##TCGA data has apparently similar counts for isomirs, summarize by mean up
# TCGAmiR.COAD$miRcomb <- sub("-[1-9]$","",row.names(TCGAmiR.COAD))
# miR_COAD_comb <- as.data.frame(TCGAmiR.COAD %>% dplyr::group_by(miRcomb) %>% 
#                                  dplyr::summarise_if(.predicate = function(x) is.numeric(x),
#                                                      .funs = mean)) #affects around 180miRNA
# row.names(miR_COAD_comb) <- miR_COAD_comb$miRcomb #afterwards remove column 1
# miR_COAD_comb$miRcomb <- NULL
# ##keep only one sample per patient
# colnames(miR_COAD_comb) <- sub("-01A-.*","",colnames(miR_COAD_comb))
# miR_COAD_comb <- miR_COAD_comb[,-grep("-.+-.+-.+-.+", colnames(miR_COAD_comb) )] #get rid of 20 normal and extra samples
# miR_COAD_comb <- miR_COAD_comb[,!duplicated(colnames(miR_COAD_comb))] #no more duplicated patients
# ##normalize, in right direction for RF
# miR_COAD_vst <- as.data.frame( t( varianceStabilizingTransformation( 
#   as.matrix( round( miR_COAD_comb, 0 ) ) ) ) )
# colnames(miR_COAD_vst) <- gsub("-", ".", colnames(miR_COAD_vst))
# 
# miR_COAD_qn <-as.data.frame( t(quantile_normalisation(miR_COAD_comb))) # classifier not improved
# colnames(miR_COAD_qn) <- gsub("-", ".", colnames(miR_COAD_qn))
# 
# 
# ####combine with labels
# ##get the CMS labels from CMSclassifier 
# rownames(miR_COAD_vst) <- sub("//.*","",rownames(miR_COAD_vst))
# miR_COAD_vst$CMS.cl.rf <- predictedCMS.COAD[match(rownames(miR_COAD_vst),
#                                                         rownames(predictedCMS.COAD)),"RF"]
# #miR_COAD_vst$CMS.cl.rf[is.na(miR_COAD_vst$CMS.cl.rf)] <- "NOLBL"
# miR_COAD_vst$CMS.cl.rf <- factor(miR_COAD_vst$CMS.cl.rf)
# 
# ##get the CMS labels from Louis table
# miR_COAD_vst$CMS.lv <- CMS_samples$CMS[match(rownames(miR_COAD_vst), 
#                                              CMS_samples$SampleId)]#get the CMS labels from Louis table
# #miR_COAD_vst$CMS.lv[which(is.na(miR_COAD_vst$CMS.lv))] <- "NOLBL"
# miR_COAD_vst$CMS.lv <- factor(miR_COAD_vst$CMS.lv)
# miR_COAD_vst$CMS.guin <- CMS_samples$CMS_final_netw_RF[match(rownames(miR_COAD_vst), 
#                                              CMS_samples$SampleId)]#get the CMS labels from Louis table
# miR_COAD_vst$CMS.guin[which(miR_COAD_vst$CMS.guin == "NOLBL")] <- NA
# miR_COAD_vst$CMS.guin <- droplevels(miR_COAD_vst$CMS.guin)
# miR_COAD_vst$CMS <- miR_COAD_vst$CMS.cl.rf
# miR_COAD_vst$CMS[-which(miR_COAD_vst$CMS.cl.rf==
#                           miR_COAD_vst$CMS.guin & miR_COAD_vst$CMS.guin==
#                           miR_COAD_vst$CMS.lv)] <- NA #keep only samples where the predictions were the same
# write.table(miR_COAD_vst, "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/miR_COAD_vst_CMS.txt")
# ### batch effect removal (AA batch is different)
# #limma was suggested also by Michael Love, https://support.bioconductor.org/p/76099/
# #rc_CMS <- relevel(rc_CMS, "CMS2")
# tssAnnot <- data.frame("CMS"=miR_COAD_vst$CMS.cl.rf) #use the predicted labels here in order to have labels for all
# tssAnnot$CMS <- factor(tssAnnot$CMS, levels=c("CMS1", "CMS2", "CMS3", "CMS4", "NOLBL"))
# tssAnnot$CMS[which(is.na(tssAnnot$CMS))] <- "NOLBL"#use the predicted labels here in order to have labels for all
# 
# # pre Batch effect removal pre outlier removal
# rownames(tssAnnot) <- rownames(miR_COAD_vst)
# tssAnnot$TSS <- sub("-.*","", sub( "TCGA-", "", rownames(tssAnnot) ) )
# tssAnnot$TSS2 <- factor(tssAnnot$TSS, labels=c(rep("other",5),"AA",
#                                             rep("other", 19 )))#only AA vs other
# designM <- model.matrix( ~ 0+ tssAnnot$CMS ) ## do i need an intercept? otherwise it compares all to CMS1 or CMS2 if releveled? 
# rc_vst_limBatchRem <- t(limma::removeBatchEffect(t(miR_COAD_vst[,grep("hsa", 
#                                                                 colnames(miR_COAD_vst))]), #the entire dataset, 339 samples
#                                                  design = designM,
#                                                  batch=tssAnnot$TSS) )
# 
# #function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects
# rc_vst_limBatchRem <- as.data.frame(rc_vst_limBatchRem)
# rc_vst_limBatchRem$CMS <- miR_COAD_vst$CMS
# 
# ##outlier detection
# tssAnnot$mean <- rowMeans(rc_vst_limBatchRem[,grep("hsa", colnames(rc_vst_limBatchRem))])
# tssAnnot$sums <- rowSums(rc_vst_limBatchRem[,grep("hsa", colnames(rc_vst_limBatchRem))])
# tssAnnot$sDs <- rowSds(as.matrix(rc_vst_limBatchRem[,grep("hsa", colnames(rc_vst_limBatchRem))]))
# tssAnnot$outls <- is_outlier(tssAnnot$sDs)
# outlierS <- rownames(tssAnnot[is_outlier(tssAnnot$sDs),])#outliers with low mean AND high SD (across all samples)
# outlierS
# 
# # ## meanSD Plot should show horizontal line
# # pdf(paste0(plotDir, "meanRankSDPlot_COAD.pdf"), 
# #     onefile = T)
# # #before: 
# # meanSdPlot(t(as.matrix(log2(t(miR_COAD_comb[grep("hsa", 
# #                                                 rownames(miR_COAD_comb)),]+1)))),
# #            ylab= "sd log2 raw")
# # #qn:
# # meanSdPlot(t(as.matrix(log2(miR_COAD_qn[!tssAnnot$outls,grep("hsa", 
# #                                              colnames(miR_COAD_qn))]+1))), 
# #            ylab= "sd log2 QN")
# # #vst:
# # meanSdPlot(t(as.matrix(miR_COAD_vst[!tssAnnot$outls,grep("hsa", 
# #                                           colnames(miR_COAD_vst))])),
# #            ylab= "sd vst") #peak sd is pushed under 1
# # #after vst & BR
# # meanSdPlot(t(as.matrix(rc_vst_limBatchRem[!tssAnnot$outls,grep("hsa", 
# #                                                     colnames(rc_vst_limBatchRem))])),
# #            ylab= "sd vst BR") #fewer features above 1.5 sd
# # dev.off()
# # 
# 
# ##exploratory PCA
# COAD_trans <- preProcess(as.data.frame(rc_vst_BR),
#                        method=c(#"BoxCox",
#                          "center","scale", 
#                          "pca"),
#                        thresh = list(thresh = 0.60))
# COAD_PC <- predict(COAD_trans, rc_vst_BR)
# 
# library(ggsci)
# pdf(paste0(plotDir, "COAD_miR_rc_vst_BR_TSS_PCA.pdf"),
#     onefile = T)
# ggplot(COAD_PC,aes(x=PC1,y=PC2,
#                  colour= tssAnnot$TSS[match(rownames(COAD_PC), rownames(tssAnnot))]))+ #clinical.coad$race[match(rownames(tssAnnot),
#                                #                          clinical.coad$submitter_id)]))+#sub("[a,b,c]$","",clinical.coad$tumor_stage[match(rownames(tssAnnot),
#                                #clinical.coad$submitter_id)])))+
#   geom_point(na.rm = F) +
#   theme_minimal() +
#    scale_color_manual(name="TSS",values = wes_palette("Darjeeling1", 25,
#                                           type = c( "continuous")) )
#   #scale_color_npg(name="Ethnicity")#lancet() #jco # npg # aaas
# dev.off()
# 
# #for plot annotation after BR
# tssAnnot$meanBR <- rowMeans(rc_vst_limBatchRem[, grep("hsa",
#                                             colnames(rc_vst_limBatchRem))])
# tssAnnot$sumsBR <- rowSums(rc_vst_limBatchRem[, grep("hsa",
#                                            colnames(rc_vst_limBatchRem))])
# tssAnnot$SDsBR <- rowSds(as.matrix(rc_vst_limBatchRem[, grep("hsa",
#                                                    colnames(rc_vst_limBatchRem))]))
# ## for plot annotation before BR
# tssAnnot$meanVST <- rowMeans(miR_COAD_vst[, grep("hsa",
#                                                   colnames(miR_COAD_vst))])
# tssAnnot$sumsVST <- rowSums(miR_COAD_vst[, grep("hsa",
#                                                  colnames(miR_COAD_vst))])
# tssAnnot$SDsVST <- rowSds(as.matrix(miR_COAD_vst[, grep("hsa",
#                                                          colnames(miR_COAD_vst))]))
# ## for plot before BR with quantile normalization
# tssAnnot$meanQN <- rowMeans(miR_COAD_qn[, grep("hsa",
#                                                 colnames(miR_COAD_qn))])
# tssAnnot$sumsQN <- rowSums(miR_COAD_qn[, grep("hsa",
#                                                colnames(miR_COAD_qn))])
# tssAnnot$SDsQN <- rowSds(as.matrix(miR_COAD_qn[, grep("hsa",
#                                                        colnames(miR_COAD_qn))]))
# 
# ## visualize batch effects based on boxplots of samples or batches:
# library(ggrepel)
# library(wesanderson)
# ggplot(tssAnnot, aes(x=TSS, y=SDsVST, colour=TSS)) +
#   geom_boxplot(outlier.size=0.8) +
#   geom_jitter(shape=16, position=position_jitter(0.2), size=0.6) +
#   ggtitle("SDs of rc_vst pre BR") +
#   theme_classic() +
#   #geom_text_repel(aes(label = outls), na.rm = TRUE) +
#   scale_color_manual(values = wes_palette("Darjeeling1", 25,
#                                           type = c( "continuous")) ) ## used this for the TSS
# ggsave(paste0(plotDir, "rc_vst_preBR_TSS_boxplotSDs_allSamples.pdf"),
#        width=6, height=6)
# 
# ### remove outliers and uninformative CMS
# rc_vst_BR <- rc_vst_limBatchRem[grep("CMS",rc_vst_limBatchRem$CMS),]
# rc_vst_BR <- rc_vst_BR[-which(rownames(rc_vst_BR) %in% outlierS), ]
# dim(rc_vst_BR)
# # ### find threshold for most discriminative features
# # hist(colSums(as.matrix(rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))])), 500)
# # densityplot(colVars(as.matrix(rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))])))
# # #biphasic peak is the most pronounced in variances, eyeballing threshold for feature preselection.
# feature_names <- colnames(rc_vst_BR[,which(colVars(as.matrix(
#   rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))]))>0.5)])
# rc_vst_BR <- rc_vst_BR[,c(feature_names, "CMS")]
# 
# ##identifying low info features, plot their correlation
# nzv <- nearZeroVar(rc_vst_BR, saveMetrics=TRUE)
# which(nzv$nzv) #no features to remove
# descrCor <-  cor(rc_vst_BR[, setdiff(grep("hsa", colnames(rc_vst_BR)), 
#                                      which(nzv$nzv))]) #exclude zero var
# library(corrplot)
# highlyCorDescr <- unique(c(colnames(descrCor[,findCorrelation(descrCor, cutoff = .75, exact=T)]), 
#                            colnames(rc_vst_BR[,nzv$nzv])))
# ### plot correlating features
# corrplot(descrCor[highlyCorDescr, 
#                   highlyCorDescr], type="upper", order="hclust",
#          col=brewer.pal(n=8, name="RdYlBu"), addgrid.col = NA, tl.cex = .5,
#          tl.col = "#999999")
# ### exclude higly correlating(redundant) features
# rc_vst_BR <- rc_vst_BR[,setdiff(colnames(rc_vst_BR),highlyCorDescr)]
# dim(rc_vst_BR)
# 
# write.table(rc_vst_BR, "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/rc_vst_BR-outl_mostVar-highlyCorr.txt")

rc_vst_BR <- read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/rc_vst_BR-outl_mostVar-highlyCorr.txt")


###### check ICGC for overlap or extra samples: ######
icgc.coad.tx <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/rnaseq.extended.metadata.aliquot_id.V4_colo.csv",
          sep="\t", header = T)
icgc.coad.tx$TCGA.patID <- sub("-01A-.*","",icgc.coad.tx$sample_id)
icgc.coad.mir <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/DCC_PCAWG_transcriptome_miRNA_sample_info_COAD.txt",
                          sep="\t", header = F)
icgc.coad.tx$TCGA.patID[which(!icgc.coad.tx$TCGA.patID %in% clinical$submitter_id & #new to TCGA set
                                          icgc.coad.tx$submitted_donor_id %in% icgc.coad.mir$V1)] # has also miR
#->none


###### get purity estimates for TCGA samples: 
purity.ABS <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                          sep="\t", header = T, dec=".",
                        colClasses = c("character","character", "factor","numeric",
                                       "numeric", "numeric", "numeric", "numeric",
                                       "numeric", "factor"))
purity.ABS$TCGA.patID <- sub("-01$","",purity.ABS$array)

purity.ABS.coad <- purity.ABS[match(rc_vst_BR$Sample.ID, purity.ABS$TCGA.patID), ] 




############# VALIDATION SETS ############

### get VUdata, it has one row by p3/p5, summarize by summing up
VUdata <- read.csv( "data/raw/Main_merged_rounded_okt_2014.txt",
                  sep = "\t", colClasses = c(rep("character", 2), rep("numeric", 221)))
VUdata_miR <- as.data.frame(VUdata[,] %>% dplyr::group_by(miRNA) %>% 
  dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                      .funs = sum))
VUdata_miR$miRNA <- sub("-1_.+$","",VUdata_miR$miRNA) #crop names end from summarized isomirs
VUdata_miR <- VUdata_miR[-grep("chr",VUdata_miR$miRNA),] #kick VU-defined miRs, they are not in TCGA anyways

# ## some exploration on feature overlap...
# VUdata_miR$miRNA[-which(VUdata_miR$miRNA %in% rownames(miR_COAD_comb))]
# rownames(miR_COAD_comb)[-which(rownames(miR_COAD_comb) %in% VUdata_miR$miRNA)]
# commonMir <- intersect(VUdata_miR$miRNA, rownames(miR_COAD_comb))
# ##compare distributions of total counts: look pretty similar
# hist(log(VUdata_miR[commonMir, "total"]+1), 30)
# hist(log(miR_COAD_comb[commonMir, "total"]+1), 30)

## adapt miR names
rownames(VUdata_miR) <- make.names(VUdata_miR$miRNA, unique=T)
VUdata_miR$miRNA <- sub(".mir.",".mir.",VUdata_miR$miRNA)
VUdata_miR$miRNA <- sub(".let.",".let.",VUdata_miR$miRNA)
VUdata_miR$miRcomb <- sub("-.+","",VUdata_miR$miRNA) 

###summarize the different isoforms
#View(VUdata_miR[!isUnique(VUdata_miR$miRcomb),]) #isomirs are often similar but if there's also the main gene it' s much higher so try sum
VUdata_miR_comb <- as.data.frame(VUdata_miR %>% dplyr::group_by(miRcomb) %>% 
                              dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                  .funs = sum)) #affects around 180miRNA
row.names(VUdata_miR_comb) <- VUdata_miR_comb$miRcomb

### apply vst 
VU_rc_vst <- VUdata_miR_comb[grep("hsa", rownames(VUdata_miR_comb)), -c(1,2)]
#VU_rc_vst <- as.data.frame( t( quantile_normalisation(VU_rc_vst)))
VU_rc_vst <- as.data.frame( t( varianceStabilizingTransformation( as.matrix( 
  round( VU_rc_vst, 0 ) ) ) ) ) #in right direction for RF

###get clinical data
clinVU <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/Main_merged_rounded_okt_2014 lokatie van biopt en tumor percentage.csv",
                        strip.white = T, sep=";", blank.lines.skip = T, nrows = 220)[,1:3]
#all patients for miR-data
colnames(clinVU) <- c("sampleID", "locSample","Tu%")
clinVU <- clinVU[match(rownames(VU_rc_vst), clinVU$sampleID), ] #sort like miR-data

clinVU.Surv <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/clinDataVU_manuallyCombined.csv",
                    sep="\t") #extensive clin data
clinVU.MSS <- read.csv2("/Users/ronjaadam/projects/miRNA_mCRC/Identifier-OSknown_DPTB.csv",
                        sep="\t")

#match sample names to patient ID in extended clin info table
clinVU.Surv$sampleID <- rownames(VU_rc_vst)[
  match(sub(" .*","",clinVU.Surv$ID),
        sub("\\..*","",rownames(VU_rc_vst)) )]
clinVU <- left_join(clinVU, clinVU.Surv, by=c("sampleID"))
clinVU <- left_join(clinVU, clinVU.MSS, by=c("sampleID"))
clinVU$Date <- gsub("/",".", clinVU$Date.x)
# clinVU$DateYY <- sub(".{6}","", clinVU$Date)
# clinVU$DateMM <- sub("\\..*","", sub(".{3}","", clinVU$Date))
# clinVU$DateDD <- sub("\\..*","", clinVU$Date)
# clinVU$Date <- as.factor(paste0(clinVU$DateYY, clinVU$DateMM, clinVU$DateDD))
# clinVU$Date <- factor(clinVU$Date, levels= levels(clinVU$Date),
#                       labels = c("NA",levels(clinVU$Date)[2:27],"NA"),
#                       ordered=T) #for checking batch effect
clinVU$StageT <- sub("N.*","",clinVU$Stage)
clinVU$StageM <- sub(".*M","M",clinVU$Stage)

## simplify location data:
clinVU$sampleType <- factor( clinVU$locSample,
                             labels=c("Met_perit","Met_perit","Met_NOS",
                                      "Met_liver","normal_liver","recur_CRC",    
                                      "Met_lung","normal_lung","Met_LN",  
                                      "Met_LN","Met_LN","Met_LN", 
                                      "Met_LN", "Met_NOS", "Met_perit",      
                                      "Met_perit", "normal_colon", "Met_perit",              
                                      "Met_perit", "Met_ovary", "Met_ovary", 
                                      "normal_ovary" ,"primary_CRC", "normal_colon",       
                                      "Met_NOS", "Met_perit"))
clinVU[,c("sampleType", "SampleOrigin", "locSample")]
clinVU$sampleType <- factor(clinVU$sampleType, levels=c("primary_CRC","recur_CRC", 
                                                        "Met_perit", "Met_liver", "Met_LN",
                                                        "Met_ovary", "Met_lung", "Met_NOS",
                                                       "normal_colon" , "normal_liver",   
                                                        "normal_lung","normal_ovary"),
                            ordered = T)

##extract patient of sample for samples with paired samples
clinVU$patient <- sub("_.*","", sub("\\..*","", sub("^[P,M]", "", clinVU$sampleID)))
#add the survival data for mets from the matching primaries
clinVU[grep("Met",clinVU$sampleType), "OS"] <- clinVU[grep("CRC",clinVU$sampleType), 
                                                      "OS"][match(clinVU[grep("Met",clinVU$sampleType),"patient"],
                                                    clinVU[grep("CRC",clinVU$sampleType),"patient"])]

clinVU[grep("Met",clinVU$sampleType), "Syn_metachroom"] <- clinVU[grep("CRC",clinVU$sampleType), 
                                                      "Syn_metachroom"][match(clinVU[grep("Met",clinVU$sampleType),"patient"],
                                                                  clinVU[grep("CRC",clinVU$sampleType),"patient"])]
clinVU[grep("Met",clinVU$sampleType), "Event_OS"] <- clinVU[grep("CRC",clinVU$sampleType), 
                                                                  "Event_OS"][match(clinVU[grep("Met",clinVU$sampleType),"patient"],
                                                                                          clinVU[grep("CRC",clinVU$sampleType),"patient"])]

clinVU$perc <- factor(clinVU$`Tu%`, ordered=T, #make ordered
          levels=c("0" ,  "35" , "40",  "45",  "50",  "55" , "60",  "65" ,"<70" , "70",">70" ))
#simplified location in colon 
clinVU$LRcolon <- factor(clinVU$Location, labels=c("R","R","L","R", "L",
                                                   "Met.omentum","ND","L", "Rectum","L",
                                                   "T" ))
clinVU$LRcolon <- factor(clinVU$LRcolon, ordered = T,
                         levels=c("R","T", "L","Rectum",
                                  "Met.omentum","ND" ),
                         labels=c("right.colon", "transv.colon","left.colon", "rectum",
                                  "Met_perit", "ND"))
summary(clinVU$sampleType)

# ###convert to rpm table
# VU_rpm <- apply( VUdata_miR_comb[1:1714, -c(1,2)], 1,
#                 function(x){x/(VUdata_miR_comb["Total",-c(1,2)]/1000000)})
# VU_rpm_df <- data.frame(Reduce(rbind,VU_rpm))
# row.names(VU_rpm_df) <- row.names(VUdata_miR_comb[1:1714,])
# VU_rpm_df$total <- rowSums(VU_rpm_df, na.rm = T)
# hist(log(rpm_t_comb[commonMir, "total"]),30)
# hist(log(VU_rpm_df[commonMir, "total"]),30)
# t_VU_rpm_df <- data.frame(t(VU_rpm_df[,]))
# write_csv2(t_VU_rpm_df, "data/cleaned/RA_VU_combMir_rpm.csv") 
# rm(VUdata_miR_comb)
# rm(VUdata)
# rm(VUdata_miR)

### what is the origin 

##exploratory PCA
VU_trans <- preProcess(VU_rc_vst[grep("prim", clinVU$sampleType),1:1714],#[which(sourceS=="P"),], 
  method=c(#"BoxCox",
    "center","scale", 
    "pca"),
  thresh = list(thresh = 0.60))
VU_PC <- predict(VU_trans, VU_rc_vst[,1:1714])

#### plot kmeans class labels on PCA #####
VU_kM <- kmeans(VU_rc_vst[grep("prim", clinVU$sampleType),1:1714], 4, nstart = 20)
VU_kM$cluster <- as.factor(VU_kM$cluster)
summary(VU_kM$cluster)

## to check for batch effects
colorBy <- clinVU$Date[match(rownames(VU_rc_vst[grep("prim", clinVU$sampleType),]), 
                             clinVU$sampleID)]

pdf("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/plots/VUdata_vst_prim_kM_PCA.pdf",
    onefile = T)
ggplot(VU_PC,aes(x=PC1,y=PC2,
              colour= VU_kM$cluster))+#clinVU$sampleType)) +
  geom_point(na.rm = F) +
  theme_minimal() +
  # scale_color_manual(values = wes_palette("Darjeeling1", 24,
  #                                        type = c( "continuous")) )
  scale_color_npg()#lancet() #jco # npg # aaas
dev.off()

##exploratory tSNE
set.seed(9)  
library(Rtsne)
tsne_model_VU <- Rtsne( as.matrix( VU_rc_vst[grep("prim", 
                                                  clinVU$sampleType),1:1714] ), 
                      check_duplicates=FALSE, 
                      pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_VU <- as.data.frame(tsne_model_VU$Y) 
ggplot(d_tsne_VU, aes(x=V1, y=V2, colour=colorBy))+
  geom_point() +
  xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme_minimal()+
  #scale_color_npg()
  scale_color_manual(values = c(wes_palette("Darjeeling1", 24, 
                                                    type = c( "continuous")) ))
scale_colour_manual(values = paletteCMS)
ggsave("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/plots/VUdata_vst_tSNE_prim_date.pdf")


############
## fetch TCGA data READ:
# library(TCGAbiolinks)
clinical.read <- GDCquery_clinic(project = c("TCGA-READ"), type = "clinical")
# CMS_samples$read <- clinical$tissue_or_organ_of_origin[match(CMS_samples$SampleId, 
#                                                              as.character(clinical.read$submitter_id))]
# 
rownames(clinical.read) <- clinical.read$submitter_id
# query.exp <- GDCquery(project = c("TCGA-READ"), 
#                       legacy = T,
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       file.type = "results", # or normalized_results
#                       experimental.strategy = "RNA-Seq", 
#                       platform = "Illumina GA"#OR: GA (for READ there seem to be duplicates with Illumina HiSeq & Illumina GA)
# )
# GDCdownload(query.exp)
# TCGARnaseqSE <- as.data.frame(GDCprepare(query= query.exp, 
#                                          summarizedExperiment = F, save = F))
# TCGARnaseqRSEM.GA <- TCGARnaseqSE[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_
# query.exp <- GDCquery(project = c("TCGA-READ"), 
#                       legacy = T,
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       file.type = "results", # or normalized_results
#                       experimental.strategy = "RNA-Seq", 
#                       platform = "Illumina HiSeq"#OR: GA (for READ there seem to be duplicates with Illumina HiSeq & Illumina GA)
# )
# GDCdownload(query.exp)
# TCGARnaseqSE <- as.data.frame(GDCprepare(query= query.exp, 
#                                          summarizedExperiment = F, save = F))
# TCGARnaseqRSEM.HiSeq <- TCGARnaseqSE[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_
# 
# unique(c(colnames(TCGARnaseqRSEM.HiSeq), colnames(TCGARnaseqRSEM.GA)))
# TCGARnaseqRaw <- cbind( TCGARnaseqRSEM.HiSeq, TCGARnaseqRSEM.GA )
# colnames( TCGARnaseqRaw ) <- c(paste0("HS.",colnames(TCGARnaseqRSEM.HiSeq)), 
#                                paste0("GA.",colnames(TCGARnaseqRSEM.GA)))
# TCGARnaseqRaw <- TCGARnaseqRaw[,grep("-01A-.*", colnames(TCGARnaseqRaw))] #keep only first primary tu sample
# remove(TCGARnaseqSE) #163 unique primary READ samples when combining GA and HiSeq
# colnames(TCGARnaseqRaw) <- sub("scaled_estimate_","",colnames(TCGARnaseqRaw))
# colnames(TCGARnaseqRaw) <- sub("-01A-.*", "",colnames(TCGARnaseqRaw))
# 
# batchInfo <- data.frame("sample"=colnames(TCGARnaseqRaw),
#                         "batch"=sub("\\..*","",colnames(TCGARnaseqRaw)))
# TCGARnaseqDF.READ <- limma::removeBatchEffect(as.matrix(log2(TCGARnaseqRaw+0.000001)), #the entire dataset
#                                               batch=batchInfo$batch)
# colnames(TCGARnaseqDF.READ) <- sub("^...","",colnames(TCGARnaseqDF.READ)) #GA./HS. being the batch identifier in the colname, keeps HiSeq source in dupliciates
# TCGARnaseqDF.READ <- TCGARnaseqDF.READ[,!duplicated(colnames(TCGARnaseqDF.READ))]
# rownames(TCGARnaseqDF.READ) <- oldRowNames
#
# write.csv(data.frame("LV.table"=CMS_samples$CMS_final_netw_RF[match(rownames(CMSread.rf$predictedCMS),
#                                                                CMS_samples$SampleId)],
#                 CMSread.rf$predictedCMS,
#                 CMSread.ss$predictedCMS), row.names = T,
#             "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-READ-allRNA-01Aonly_CMS-labels.csv")
# write.csv(TCGARnaseqDF.READ,
#             "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-READ-mRNA-GA-HiSeq_scaled_estimate_01Aonly.csv")

TCGARnaseqDF.READ <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-READ-mRNA-GA-HiSeq_scaled_estimate_01Aonly.csv",
                              row.names = 1, header = T)
colnames(TCGARnaseqDF.READ) <- gsub("\\.","-",colnames(TCGARnaseqDF.READ))
## get the classification from CMSclassifier package
# library(CMSclassifier)
# rownames(TCGARnaseqDF.READ) <- sub(".*\\|", "", rownames(TCGARnaseqDF.READ))#keep only the ENTREZ for classififer
# CMSread.rf <-classifyCMS(TCGARnaseqDF.READ, method="RF")
# CMSread.ss <- classifyCMS(TCGARnaseqDF.READ, method="SSP")
# CMSsymbols <- mapIds(org.Hs.eg.db, keys=listModelGenes(method = "RF"), 
#                      'SYMBOL', 'ENTREZID')
# View(data.frame("LV.table"=CMS_samples$CMS_final_netw_RF[match(rownames(CMSread.rf$predictedCMS), 
#                                                                CMS_samples$SampleId)], 
#                 CMSread.rf$predictedCMS, 
#                 CMSread.ss$predictedCMS))
predictedCMS.READ <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-READ-allRNA-01Aonly_CMS-labels.csv",
                              row.names = 1)

###################
##get READ miRNA 
## fetch TCGA data miRNAs:
# library(TCGAbiolinks)
# query.exp <- GDCquery(project = c("TCGA-READ"), 
#                       #legacy = F,
#                       data.category = "Transcriptome Profiling",
#                       data.type = "miRNA Expression Quantification",
#                       experimental.strategy = "miRNA-Seq")
# GDCdownload(query.exp)
# TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
# TCGARnaseqDF <- as.data.frame(TCGARnaseqSE[,grep("read_count_", colnames(TCGARnaseqSE)) ]) #or normalized_count_ for normalized_results
# colnames(TCGARnaseqDF) <- sub("-01A-.*","",sub("read_count_","",colnames(TCGARnaseqDF)))
# rownames(TCGARnaseqDF) <- TCGARnaseqSE$miRNA_ID
# remove(TCGARnaseqSE)
# TCGARnaseqDF$miRcomb <- sub("-[1-9]$","",row.names(TCGARnaseqDF))
# miR_READ_comb <- as.data.frame(TCGARnaseqDF %>% dplyr::group_by(miRcomb) %>% 
#                                  dplyr::summarise_if(.predicate = function(x) is.numeric(x),
#                                                      .funs = mean)) #affects around 180miRNA
# row.names(miR_READ_comb) <- miR_READ_comb$miRcomb #afterwards remove column 1
# 
# ##normalize, in right direction for RF
# miR_READ_qn <- as.data.frame( t( quantile_normalisation(miR_READ_comb[,-1])))
# miR_READ_vst <- as.data.frame( t( varianceStabilizingTransformation( as.matrix( 
#   round( miR_READ_comb[,-1], 0 ) ) ) ) )
# colnames(miR_READ_vst) <- gsub("-", ".", colnames(miR_READ_vst))
# #get the CMS labels from CMSclassifier after vst
# miR_READ_vst <- miR_READ_vst[-grep("-.+-.+-.+-.+", rownames(miR_READ_vst)),] #get rid of normal and extra samples
# 
# rownames(miR_READ_vst) <- sub("-01A-.*","",rownames(miR_READ_vst))
# miR_READ_vst$CMS.cl.rf<- predictedCMS.READ[match(rownames(miR_READ_vst), 
#                         rownames(predictedCMS.READ)),"RF"]
# # miR_READ_vst$CMS.cl.rf <- CMSread.rf$nearestCMS[match(rownames(miR_READ_vst), 
# #                                                       rownames(CMSread.rf$nearestCMS)),"RF"]
# miR_READ_CMS.info <- as.character(miR_READ_vst$CMS.cl.rf)
# miR_READ_CMS.info[which(is.na(miR_READ_CMS.info))] <- "NA"
# 
# ##batch effect removal
# batchInfo <- data.frame("sample"=rownames(miR_READ_vst),
#                         "batch"=factor(sub("-.*","", 
#                                            sub("TCGA-","", 
#                                                rownames(miR_READ_vst)))))#[
#                                                  #grep("CMS", miR_READ_vst$CMS.cl.rf),]
# designR <- model.matrix( ~ 0 + miR_READ_CMS.info)#[
#   #grep("CMS", miR_READ_vst$CMS.cl.rf)] ) ## do i need an intercept? otherwise it compares all to CMS1 or CMS2 if releveled? 
# miR_READ_vst_BR <- as.data.frame(t(limma::removeBatchEffect(as.matrix(t(
#   miR_READ_vst[#grep("CMS", miR_READ_vst$CMS.cl.rf),
#                grep("hsa", colnames(miR_READ_vst))])), #the entire dataset
#                                               design=designR,
#                                               batch=batchInfo$batch)))
# purity.ABS.read <- purity.ABS[match(row.names(miR_READ_vst_BR), purity.ABS$TCGA.patID), ] 
# 
# 
# ##get the CMS labels from Louis table
# miR_READ_vst_BR$CMS.lv <- CMS_samples$CMS[match(rownames(miR_READ_vst_BR), 
#                                                 CMS_samples$SampleId)]#get the CMS labels from Louis table
# miR_READ_vst_BR$CMS.cl.rf <- miR_READ_vst$CMS.cl.rf#[grep("CMS", miR_READ_vst$CMS.cl.rf)]
# 
# write.csv(miR_READ_vst,
#            "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-miR_READ_vst_CMS-labels.csv")
# write.csv(miR_READ_vst_BR,
#            "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-miR_READ_vst_BR_CMS-labels.csv")
miR_READ_vst_BR <- read.csv(file=
          "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-miR_READ_vst_BR_CMS-labels.csv",
          row.names = 1)


##check batch effects
tsne_model_READ <- Rtsne( as.matrix( miR_READ_vst_BR[,]), 
                        check_duplicates=FALSE, 
                        pca=TRUE, perplexity=20, theta=0.5, dims=2)
colREAD <- sub("-.*","",sub("TCGA-","",
                            clinical.read$submitter_id[match(rownames(miR_READ_vst), 
                                 clinical.read$submitter_id)]))
ggplot(as.data.frame(tsne_model_READ$Y), 
       aes(x=V1, y=V2, color =batchInfo$batch))+
  geom_point() +
  xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme_minimal()+
  scale_color_manual(values = c(wes_palette("Darjeeling1", 24, 
                                            type = c( "continuous")) ))
ggsave(paste0(plotDir,"/READ-miR_vst_tSNE_TSS_after-BR.pdf" ))
##-> removed the AG difference while retaining mRNA-based predicted CMS differences




################################################################

####### independent data ########


###### CPTAC-2: has mir-Seq and CMS from mRNA-Seq, plus clinical data, mutations etc. 
###### but  has no survival data

query.exp <- GDCquery(project = c("CPTAC-2"), 
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification", 
                      experimental.strategy = "miRNA-Seq"
)
GDCdownload(query.exp)
getResults(query.exp )
TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
# write.csv(TCGARnaseqSE, "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/CPTAC2_GDCquery_miRNA-Seq_SE.csv",
#           row.names = F)
TCGAmiR.Cptac2 <- TCGARnaseqSE[,grep("read_count_", colnames(TCGARnaseqSE)) ] #or reads_per_million_miRNA_mapped_ for normalized_results
colnames(TCGAmiR.Cptac2) <- sub("read_count_","",colnames(TCGAmiR.Cptac2))
rownames(TCGAmiR.Cptac2) <- TCGARnaseqSE$miRNA_ID
remove(TCGARnaseqSE)
##TCGA data has apparently similar counts for isomirs, summarize by mean up
TCGAmiR.Cptac2$miRcomb <- sub("-[1-9]$","",row.names(TCGAmiR.Cptac2))
miR_Cptac_comb <- as.data.frame(TCGAmiR.Cptac2 %>% dplyr::group_by(miRcomb) %>% 
                                 dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                     .funs = mean)) #affects around 180miRNA
row.names(miR_Cptac_comb) <- miR_Cptac_comb$miRcomb #afterwards remove column 1
miR_Cptac_comb$miRcomb <- NULL
clinical.cptac <- GDCquery_clinic(project = c("CPTAC-2"), type = "clinical")
##identify how submitter id and sample id map
cptac.ids <- data.frame("coln"=colnames(miR_Cptac_comb),
                        "subm"=NA)
for(i in 1:length(colnames(miR_Cptac_comb))){
cptac.ids$subm[i] <- as.character(clinical.cptac$submitter_id[ grep(colnames(miR_Cptac_comb)[i], 
     as.character(clinical.cptac$submitter_sample_ids) ) ])
}
colnames(miR_Cptac_comb)==cptac.ids$coln

clinSupp.cptac <- read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/Vasaikar_CPTAC2_STable1_mmc1.txt", 
                             sep="\t", header = T)
##normalize, in right direction for RF
miR_Cptac_vst <- as.data.frame( t( varianceStabilizingTransformation( 
  as.matrix( round( miR_Cptac_comb[, which( cptac.ids$subm %in% clinSupp.cptac$SampleID ) ] ##only the subset from the publication suppl clin data
, 0 ) ) ) ) )
colnames(miR_Cptac_vst) <- gsub("-", ".", colnames(miR_Cptac_vst)) # format mir names for RF

miR_Cptac_qn <- as.data.frame( t( quantile_normalisation( miR_Cptac_comb[, 
                                   which( cptac.ids$subm %in% clinSupp.cptac$SampleID ) ] ##only the subset from the publication suppl clin data
                    ) ) )
colnames(miR_Cptac_qn) <- gsub("-", ".", colnames(miR_Cptac_qn)) # format mir names for RF

clinSupp.cptac$sample.label <- cptac.ids$coln[match(clinSupp.cptac$SampleID,
                                                    cptac.ids$subm)]
rownames(clinSupp.cptac) <- clinSupp.cptac$SampleID

write.csv(miR_Cptac_vst,
          "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-miR_CPTAC_vst.csv")
miR_Cptac_vst <- read.csv(file="/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-miR_CPTAC_vst.csv",
                            row.names = 1)


###GSE29623: mRNA and miRNA with survival data 
mRNA21 <- read.table(sep = "\t", header=T,
  "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE29623_FFPE_array_mRNA_miRNA/GSE29621_series_matrix_expression.txt")
rownames(mRNA21) <- mRNA21$ID_REF
mRNA21 <- mRNA21[,-1]
##clinical data mRNA
mRNA21.clin <- read.csv2(sep = "\t", header = F, as.is = T,
                         "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE29623_FFPE_array_mRNA_miRNA/GSE29621_series_matrix_clinical.txt")
row.names(mRNA21.clin) <- mRNA21.clin$V1
colnames(mRNA21.clin) <- mRNA21.clin["Sample_title",]
mRNA21.clin <- data.frame(t(mRNA21.clin[,-1]))
row.names(mRNA21.clin) <- mRNA21.clin$Sample_title

library("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c("affy_hg_u133_plus_2","entrezgene_id", "hgnc_symbol"),
  filter = "affy_hg_u133_plus_2",
  values = rownames(mRNA21),
  uniqueRows=TRUE)
mRNA21$entrez <- annotLookup$entrezgene_id[match(rownames(mRNA21), 
                                                 annotLookup$affy_hg_u133_plus_2)]
mRNA21 <- mRNA21[!is.na(mRNA21$entrez),] # won't be needed for classififer anyways
mRNA21_comb <- as.data.frame(mRNA21 %>% dplyr::group_by(entrez) %>% 
                              dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                  .funs = mean)) 
row.names(mRNA21_comb) <- mRNA21_comb$entrez #afterwards remove column 1 
mRNA21_comb <- mRNA21_comb[,match(rownames(miR22_qn), colnames(mRNA21_comb))]
##TODO normalizations
CMS.gse23.ss <- classifyCMS(mRNA21_comb[,-1], method="RF")
summary(factor(CMS.gse23.ss$predictedCMS$RF))
CMS.gse23.ss <- classifyCMS(mRNA21_comb[,-1], method="SSP")
summary(factor(CMS.gse23.ss$predictedCMS$SSP))
mRNA21.clin$CMSmRNA <- factor(CMS.gse23.ss$predictedCMS$SSP)

###mir data
miR22 <- read.table(sep = "\t", header=T,
                     "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE29623_FFPE_array_mRNA_miRNA/GSE29622_series_matrix_expression.txt")
miR22.B <- read.table(sep = "\t", header=T,
                                 "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE29623_FFPE_array_mRNA_miRNA/GSE29622_non-normalized_B.txt")
rownames(miR22) <- miR22$ID_REF
miR22 <- miR22[,-1]
View(miR22[-which(isUnique(sub("-[0-9]+$","",sub("\\*","",miR22$ID_REF)))), 1:10])

miR22.clin <- read.csv2(sep = "\t", header = F, as.is = T,
                         "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE29623_FFPE_array_mRNA_miRNA/GSE29622_series_matrix_clinical.txt")
row.names(miR22.clin) <- miR22.clin$V1
miR22.clin["Series_sample_id",] <- strsplit(miR22.clin["Series_sample_id",1], " ")[[1]]
colnames(miR22.clin) <- miR22.clin["Sample_geo_accession",]
miR22.clin <- data.frame(t(miR22.clin[,-1]))
row.names(miR22.clin) <- miR22.clin$Sample_title
summary(miR22.clin$Sample_characteristics_ch1_M)
rownames(miR22.clin)==rownames(mRNA21.clin) ##same order for both datasets
as.character(miR22.clin$Sample_geo_accession)==colnames(miR22) ##same order for both datasets
as.character(mRNA21.clin$Sample_geo_accession)==colnames(mRNA21) ##same order for both datasets

### combine isomirs to match TCGA (only one row per mir gene)
miR22 <- miR22[-which(rowMeans(miR22) <= 1 ),] #-310
miR22$miRcomb <- sub("\\*?-[0-9]+$","",row.names(miR22))
miR22$miRcomb[grep("-.*-.*-.*", miR22$miRcomb)] <- sub("-[1-9]p?$","",miR22$miRcomb[grep("-.*-.*-.*", miR22$miRcomb)])
miR22$miRcomb[grep("-.*-.*-.*", miR22$miRcomb)] <- sub("-[1-9]p?$","",miR22$miRcomb[grep("-.*-.*-.*", miR22$miRcomb)])
miR22$miRcomb <- gsub("-",".",miR22$miRcomb)
miR22$miRcomb <- gsub("miR","mir",miR22$miRcomb)

miR22_comb <- as.data.frame(miR22 %>% dplyr::group_by(miRcomb) %>% 
                                 dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                     .funs = mean)) #affects around 210miRNA
row.names(miR22_comb) <- miR22_comb$miRcomb #afterwards remove column 1 
##normalize, in right direction for RF
miR22_qn <- as.data.frame( t( quantile_normalisation(miR22_comb[,-1])))



### GSE18392: microarray miRNA data with dMMR/pMMR and stage information






#### GSE121842: miRNA from 3 patients, 3 CRC, 3 Normal colon and theres also mRNA for each
miR42 <- read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE121842_mRNA-miRNA/GSE121842_miRNA.All.Counts.exp.txt",
                    colClasses = c(rep("character", 3), rep("integer", 6)), skip = 1)
colnames(miR42) <- as.character(read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE121842_mRNA-miRNA/GSE121842_miRNA.All.Counts.exp.txt",
                                          as.is = T )[1,])
miR42$miRNA <- sub("-[3,5]p","",miR42$miRNAName) #ignoring the strands
miR42$miRNA[-grep("miR-.$",miR42$miRNA)] <- sub("-.$","",miR42$miRNA[-grep("miR-.$",miR42$miRNA)])
miRC42 <- as.data.frame(miR42 %>% dplyr::group_by(miRNA) %>% 
                              dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                  .funs = sum)) #summarizing isoforms and strands
rownames(miRC42) <- sub("miR","mir",gsub("-",".",miRC42$miRNA)) #make name format comparable

miRC42_vst <- as.data.frame(t(DESeq2::varianceStabilizingTransformation(
  as.matrix(miRC42[,grep("C", colnames(miRC42))])))) #vst (only cancer samples) and transpose
View(miRC42_vst)

library(devtools)
install_github("Sage-Bionetworks/CMSclassifier")
library(CMSclassifier)
mRNA42 <- read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE121842_mRNA-miRNA/GSE121842_all.counts.exp.txt",
                   skip = 1)
colnames(mRNA42) <- as.character(read.table("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/GSE121842_mRNA-miRNA/GSE121842_all.counts.exp.txt",
                                           as.is = T )[1,])
rownames(mRNA42) <- mRNA42$GeneID
mRNA42_vst <- as.data.frame(DESeq2::varianceStabilizingTransformation(
  as.matrix(mRNA42[,grep("C", colnames(mRNA42))])))
library('org.Hs.eg.db')
entrezIds <- mapIds(org.Hs.eg.db, keys=rownames(mRNA42), 'ENTREZID', 'SYMBOL')
mRNA42_vst <- mRNA42_vst[-which(is.na(entrezIds)),]

rownames(mRNA42_vst) <- entrezIds[-which(is.na(entrezIds))]
CMS42 <- CMSclassifier::classifyCMS(mRNA42_vst,method="RF")
CMS42ss <- CMSclassifier::classifyCMS(mRNA42_vst,method="SSP")
