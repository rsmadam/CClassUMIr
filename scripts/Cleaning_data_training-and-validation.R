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

##function for quantile normalisation 
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

##function for outlier detection
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
} #Tukey's mild outlier definition



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


############## training data COAD #############

############# get COAD raw mRNA from TCGA ################
clinical.coad <- GDCquery_clinic(project = c("TCGA-COAD"), type = "clinical")
oldRowNames <- rownames(read.table(paste0(projDir,"/Data/TCGA-COAD-allRNA-01Aonly-RSEM-ComBat.csv"),
                                   header = T)) ##genes
#after packages update, rownames are gone. using the old rownames from a saved table because the number of rows is the same

rownames(clinical.coad) <- clinical.coad$submitter_id
query.exp <- GDCquery(project = c("TCGA-COAD"),
                      legacy = T,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      file.type = "results", # or normalized_results
                      experimental.strategy = "RNA-Seq",
                      platform = "Illumina GA"
)
GDCdownload(query.exp)
TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
TCGARnaseqRSEM.GA <- as.data.frame(TCGARnaseqSE)[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_
query.exp <- GDCquery(project = c("TCGA-COAD"),
                      legacy = T,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      file.type = "results", # or normalized_results
                      experimental.strategy = "RNA-Seq",
                      platform = "Illumina HiSeq"#OR: GA (for READ there seem to be duplicates with Illumina HiSeq & Illumina GA)
)
GDCdownload(query.exp)
TCGARnaseqSE <- as.data.frame(GDCprepare(query= query.exp, summarizedExperiment = F, save = F))
TCGARnaseqRSEM.HiSeq <- as.data.frame(TCGARnaseqSE)[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_

##for the CMSclassifier use RSEM normalization, log transform and remove batch effect GA/HiSeq
unique(c(colnames(TCGARnaseqRSEM.HiSeq), colnames(TCGARnaseqRSEM.GA)))
TCGARnaseqRaw <- cbind( TCGARnaseqRSEM.HiSeq, TCGARnaseqRSEM.GA )
colnames( TCGARnaseqRaw ) <- c(paste0("HS.",colnames(TCGARnaseqRSEM.HiSeq)),
                               paste0("GA.",colnames(TCGARnaseqRSEM.GA)))
TCGARnaseqRaw <- TCGARnaseqRaw[,grep("-01A-.*", colnames(TCGARnaseqRaw))] #keep only first primary tu sample
remove(TCGARnaseqSE) #455 unique COAD patients when combining GA and HiSeq
colnames(TCGARnaseqRaw) <- sub("scaled_estimate_","",colnames(TCGARnaseqRaw))
colnames(TCGARnaseqRaw) <- sub("-01A-.*", "",colnames(TCGARnaseqRaw))
batchInfo <- data.frame("sample"=colnames(TCGARnaseqRaw),
                        "batch"=sub("\\..*","",colnames(TCGARnaseqRaw)))
TCGARnaseqDF.COAD <- limma::removeBatchEffect(as.matrix(log2(TCGARnaseqRaw+0.000001)), #the entire dataset
                                              batch=batchInfo$batch)
colnames(TCGARnaseqDF.COAD) <- sub("^...","",colnames(TCGARnaseqDF.COAD)) #GA./HS. being the batch identifier in the colname, keeps HiSeq source in dupliciates
TCGARnaseqDF.COAD <- TCGARnaseqDF.COAD[,!duplicated(colnames(TCGARnaseqDF.COAD))]
rownames(TCGARnaseqDF.COAD) <- oldRowNames

## get the classification from CMSclassifier package
rownames(TCGARnaseqDF.COAD) <- sub(".*\\|", "", rownames(TCGARnaseqDF.COAD))#keep only the ENTREZ for classififer
CMScoad.rf <-classifyCMS(TCGARnaseqDF.COAD, method="RF")
CMScoad.ss <- classifyCMS(TCGARnaseqDF.COAD, method="SSP")
CMSsymbols <- mapIds(org.Hs.eg.db, keys=listModelGenes(method = "RF"),
                     'SYMBOL', 'ENTREZID')
write.csv(data.frame("LV.table"=CMS_samples$CMS_final_netw_RF[match(rownames(CMScoad.rf$predictedCMS),
                                                               CMS_samples$SampleId)],
                CMScoad.rf$predictedCMS,
                CMScoad.ss$predictedCMS), row.names = T,
            paste0(projHome, "miRNA_classifier/Data/TCGA-COAD-allRNA-01Aonly_CMS-labels.csv"))
write.csv(TCGARnaseqDF.COAD,
            paste0(projHome, "miRNA_classifier/Data/TCGA-COAD-mRNA-GA-HiSeq_scaled_estimate_01Aonly_BR.csv"))
TCGARnaseqDF.COAD <- read.csv(paste0(projDir,"/Data/TCGA-COAD-mRNA-GA-HiSeq_scaled_estimate_01Aonly_BR.csv"),
                                   row.names = 1)
predictedCMS.COAD <- read.csv(paste0(projDir, "/Data/TCGA-COAD-allRNA-01Aonly_CMS-labels.csv"),
                    row.names = 1)


################### get COAD miRNA from TCGA ################

## fetch TCGA COAD miRNAs:
query.exp <- GDCquery(project = c("TCGA-COAD"),
                      legacy = F,
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification",
                      experimental.strategy = "miRNA-Seq"
)
GDCdownload(query.exp)
TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
TCGAmiR.COAD <- TCGARnaseqSE[,grep("read_count_", colnames(TCGARnaseqSE)) ] #or normalized_count_ for normalized_results
colnames(TCGAmiR.COAD) <- sub("read_count_","",colnames(TCGAmiR.COAD))
rownames(TCGAmiR.COAD) <- TCGARnaseqSE$miRNA_ID
remove(TCGARnaseqSE)
##TCGA data has similar counts for isomirs, summarize by mean up
TCGAmiR.COAD$miRcomb <- sub("-[1-9]$","",row.names(TCGAmiR.COAD))
miR_COAD_comb <- as.data.frame(TCGAmiR.COAD %>% dplyr::group_by(miRcomb) %>%
                                 dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                     .funs = mean)) #affects around 180miRNA
row.names(miR_COAD_comb) <- miR_COAD_comb$miRcomb #afterwards remove column 1
miR_COAD_comb$miRcomb <- NULL
##keep only one sample per patient
colnames(miR_COAD_comb) <- sub("-01A-.*","",colnames(miR_COAD_comb))
miR_COAD_comb <- miR_COAD_comb[,-grep("-.+-.+-.+-.+", colnames(miR_COAD_comb) )] #get rid of 20 normal and extra samples
miR_COAD_comb <- miR_COAD_comb[,!duplicated(colnames(miR_COAD_comb))] #no more duplicated patients
##normalize, in right direction for RF
miR_COAD_vst <- as.data.frame( t( varianceStabilizingTransformation(
  as.matrix( round( miR_COAD_comb, 0 ) ) ) ) )
colnames(miR_COAD_vst) <- gsub("-", ".", colnames(miR_COAD_vst))
##also created quantile normalization, just to check
miR_COAD_qn <-as.data.frame( t(quantile_normalisation(miR_COAD_comb))) # classifier not improved
colnames(miR_COAD_qn) <- gsub("-", ".", colnames(miR_COAD_qn))



########## CMS labels ############
####combine miR data with labels
##get the CMS labels from CMSclassifier
rownames(miR_COAD_vst) <- sub("//.*","",rownames(miR_COAD_vst))
miR_COAD_vst$CMS.cl.rf <- predictedCMS.COAD[match(rownames(miR_COAD_vst),
                                                        rownames(predictedCMS.COAD)),"RF"]
miR_COAD_vst$CMS.cl.rf <- factor(miR_COAD_vst$CMS.cl.rf)

##get the CMS labels Guinney et al. publication
miR_COAD_vst$CMS.lv <- CMS_samples$CMS[match(rownames(miR_COAD_vst),
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
rc_vst_limBatchRem <- t(limma::removeBatchEffect(t(miR_COAD_vst[,grep("hsa",
                                                                colnames(miR_COAD_vst))]), #the entire dataset, 339 samples
                                                 design = designM,
                                                 batch=tssAnnot$TSS) )
#function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects
rc_vst_limBatchRem <- as.data.frame(rc_vst_limBatchRem)
rc_vst_limBatchRem$CMS <- miR_COAD_vst$CMS


######### further quality control of miR data ########
##statistics outlier detection
tssAnnot$mean <- rowMeans(rc_vst_limBatchRem[,grep("hsa", colnames(rc_vst_limBatchRem))])
tssAnnot$sums <- rowSums(rc_vst_limBatchRem[,grep("hsa", colnames(rc_vst_limBatchRem))])
tssAnnot$sDs <- rowSds(as.matrix(rc_vst_limBatchRem[,grep("hsa", colnames(rc_vst_limBatchRem))]))
tssAnnot$outls <- is_outlier(tssAnnot$sDs)

#### meanSD Plot (should show horizontal line)
pdf(paste0(plotDir, "meanRankSDPlot_COAD.pdf"),
    onefile = T)
#before:
meanSdPlot(t(as.matrix(log2(t(miR_COAD_comb[grep("hsa",
                                                rownames(miR_COAD_comb)),]+1)))),
           ylab= "sd log2 raw")
#vst:
meanSdPlot(t(as.matrix(miR_COAD_vst[!tssAnnot$outls,grep("hsa",
                                          colnames(miR_COAD_vst))])),
           ylab= "sd vst") #peak sd is pushed under 1
#after vst & BR
meanSdPlot(t(as.matrix(rc_vst_limBatchRem[!tssAnnot$outls,grep("hsa",
                                                    colnames(rc_vst_limBatchRem))])),
           ylab= "sd vst BR") #fewer features above 1.5 sd
dev.off()

#### exploratory PCA
COAD_trans <- preProcess(as.data.frame(rc_vst_BR),
                       method=c(#"BoxCox",
                         "center","scale",
                         "pca"),
                       thresh = list(thresh = 0.60))
COAD_PC <- predict(COAD_trans, rc_vst_BR)

pdf(paste0(plotDir, "COAD_miR_rc_vst_BR_TSS_PCA.pdf"),
    onefile = T)
ggplot(COAD_PC,aes(x=PC1,y=PC2,
                 colour= tssAnnot$TSS[match(rownames(COAD_PC), rownames(tssAnnot))]))+ #clinical.coad$race[match(rownames(tssAnnot),
                               #                          clinical.coad$submitter_id)]))+#sub("[a,b,c]$","",clinical.coad$tumor_stage[match(rownames(tssAnnot),
                               #clinical.coad$submitter_id)])))+
  geom_point(na.rm = F) +
  theme_minimal() +
   scale_color_manual(name="TSS",values = wes_palette("Darjeeling1", 25,
                                          type = c( "continuous")) )
  #scale_color_npg(name="Ethnicity")#lancet() #jco # npg # aaas
dev.off()

## for plot annotation after BR
tssAnnot$meanBR <- rowMeans(rc_vst_limBatchRem[, grep("hsa",
                                            colnames(rc_vst_limBatchRem))])
tssAnnot$sumsBR <- rowSums(rc_vst_limBatchRem[, grep("hsa",
                                           colnames(rc_vst_limBatchRem))])
tssAnnot$SDsBR <- rowSds(as.matrix(rc_vst_limBatchRem[, grep("hsa",
                                                   colnames(rc_vst_limBatchRem))]))
## for plot annotation before BR
tssAnnot$meanVST <- rowMeans(miR_COAD_vst[, grep("hsa",
                                                  colnames(miR_COAD_vst))])
tssAnnot$sumsVST <- rowSums(miR_COAD_vst[, grep("hsa",
                                                 colnames(miR_COAD_vst))])
tssAnnot$SDsVST <- rowSds(as.matrix(miR_COAD_vst[, grep("hsa",
                                                         colnames(miR_COAD_vst))]))

#### visualize batch effects based on boxplots of samples or batches:
ggplot(tssAnnot, aes(x=TSS, y=SDsVST, colour=TSS)) +
  geom_boxplot(outlier.size=0.8) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.6) +
  ggtitle("SDs of rc_vst pre BR") +
  theme_classic() +
  #geom_text_repel(aes(label = outls), na.rm = TRUE) +
  scale_color_manual(values = wes_palette("Darjeeling1", 25,
                                          type = c( "continuous")) ) ## used this for the TSS
ggsave(paste0(plotDir, "rc_vst_preBR_TSS_boxplotSDs_allSamples.pdf"),
       width=6, height=6)



######### cleaning of miR data ########

## remove outliers and uninformative CMS
outlierS <- rownames(tssAnnot[is_outlier(tssAnnot$sDs),])#outliers with low mean AND high SD (across all samples)
outlierS
rc_vst_BR <- rc_vst_limBatchRem[grep("CMS",rc_vst_limBatchRem$CMS),]
rc_vst_BR <- rc_vst_BR[-which(rownames(rc_vst_BR) %in% outlierS), ]
dim(rc_vst_BR)

## find threshold for most discriminative features
hist(colSums(as.matrix(rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))])), 500)
densityplot(colVars(as.matrix(rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))])))
#biphasic peak is the most pronounced in variances, eyeballing threshold for feature preselection.
feature_names <- colnames(rc_vst_BR[,which(colVars(as.matrix(
  rc_vst_BR[,grep("hsa", colnames(rc_vst_BR))]))>0.5)])
rc_vst_BR <- rc_vst_BR[,c(feature_names, "CMS")]

##identifying low info features, plot their correlation
nzv <- nearZeroVar(rc_vst_BR, saveMetrics=TRUE)
which(nzv$nzv) #no features to remove
descrCor <-  cor(rc_vst_BR[, setdiff(grep("hsa", colnames(rc_vst_BR)),
                                     which(nzv$nzv))]) #exclude zero var
highlyCorDescr <- unique(c(colnames(descrCor[,findCorrelation(descrCor, cutoff = .75, exact=T)]),
                           colnames(rc_vst_BR[,nzv$nzv])))
## plot correlating features
corrplot(descrCor[highlyCorDescr,
                  highlyCorDescr], type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"), addgrid.col = NA, tl.cex = .5,
         tl.col = "#999999")

## exclude higly correlating(redundant) features
rc_vst_BR <- rc_vst_BR[,setdiff(colnames(rc_vst_BR),highlyCorDescr)]
dim(rc_vst_BR)

## this is the final table 
write.table(rc_vst_BR, paste0(projHome, "miRNA classifier/Data/rc_vst_BR-outl_mostVar-highlyCorr.txt"))

rc_vst_BR <- read.table(paste0(projDir,"/Data/rc_vst_BR-outl_mostVar-highlyCorr.txt"))


###### get purity estimates for TCGA samples: 
purity.ABS <- read.csv2(paste0(projHome, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"),
                          sep="\t", header = T, dec=".",
                        colClasses = c("character","character", "factor","numeric",
                                       "numeric", "numeric", "numeric", "numeric",
                                       "numeric", "factor"))
purity.ABS$TCGA.patID <- sub("-01$","",purity.ABS$array)
purity.ABS.coad <- purity.ABS[match(rc_vst_BR$Sample.ID, purity.ABS$TCGA.patID), ] 
purity.ABS.read <- purity.ABS[match(rownames(miR_READ_vst_BR), purity.ABS$TCGA.patID), ] 
remove(purity.ABS)





############# VALIDATION SETS ############



#############  EGAS1127 ############
### get EGAS1127 miR data, here called VUdata, it has one row by p3/p5, summarize by summing up
VUdata <- read.csv( "data/raw/Main_merged_rounded_okt_2014.txt",
                  sep = "\t", colClasses = c(rep("character", 2), rep("numeric", 221)))
VUdata_miR <- as.data.frame(VUdata[,] %>% dplyr::group_by(miRNA) %>% 
  dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                      .funs = sum))
VUdata_miR$miRNA <- sub("-1_.+$","",VUdata_miR$miRNA) #crop names end from summarized isomirs
VUdata_miR <- VUdata_miR[-grep("chr",VUdata_miR$miRNA),] #kick VU-defined miRs, they are not in TCGA anyways
remove(VUdata)

## adapt miR names
rownames(VUdata_miR) <- make.names(VUdata_miR$miRNA, unique=T)
VUdata_miR$miRNA <- sub(".mir.",".mir.",VUdata_miR$miRNA)
VUdata_miR$miRNA <- sub(".let.",".let.",VUdata_miR$miRNA)
VUdata_miR$miRcomb <- sub("-.+","",VUdata_miR$miRNA) 

## summarize the different isoforms
#View(VUdata_miR[!isUnique(VUdata_miR$miRcomb),]) #isomirs are often similar but if there's also the main gene it' s much higher so try sum
VUdata_miR_comb <- as.data.frame(VUdata_miR %>% dplyr::group_by(miRcomb) %>% 
                              dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                  .funs = sum)) #affects around 180miRNA
row.names(VUdata_miR_comb) <- VUdata_miR_comb$miRcomb

### apply vst 
VU_rc_vst <- VUdata_miR_comb[grep("hsa", rownames(VUdata_miR_comb)), -c(1,2)]
VU_rc_vst <- as.data.frame( t( varianceStabilizingTransformation( as.matrix( 
  round( VU_rc_vst, 0 ) ) ) ) ) #in right direction for RF


######## EGAS1127 (VU) clinical data #######
### get clinical data
clinVU.Surv <- read.csv2(paste0(projHome, "clinDataVU_manuallyCombined.csv"),
                    sep="\t") #extensive clin data
clinVU.MSS <- read.csv2(paste0(projHome, "Identifier-OSknown_DPTB.csv"),
                        sep="\t")
clinVU.age1 <- read.csv2(paste0(projHome, "Identifier-OSknown age and gender.csv"),
                         sep="\t")
clinVU.age2 <- read.csv2(paste0(projHome, "Identifier-OSunknown age and gender.csv"),
                         sep="\t")
clinVU <- read.csv2(paste0(projHome, "Main_merged_rounded_okt_2014 lokatie van biopt en tumor percentage.csv"),
                    strip.white = T, sep=";", blank.lines.skip = T, nrows = 220)[,1:3]
#all patients for miR-data
colnames(clinVU) <- c("sampleID", "locSample","Tu%")
clinVU <- clinVU[match(rownames(VU_rc_vst), clinVU$sampleID), ] #sort like miR-data

#match sample names to patient ID in extended clin info table
clinVU.Surv$sampleID <- rownames(VU_rc_vst)[
  match(sub(" .*","",clinVU.Surv$ID),
        sub("\\..*","",rownames(VU_rc_vst)) )]
clinVU <- left_join(clinVU, clinVU.Surv, by=c("sampleID"))
clinVU <- left_join(clinVU, clinVU.MSS, by=c("sampleID"))
clinVU <- left_join(clinVU, clinVU.age1, by=c("sampleID"= "SampleID"))
clinVU[match(clinVU.age2$sampleID, clinVU$sampleID),
       c("Gender", "Age", "MSI.MSS", "Stage")] <- clinVU.age2[,c("Gender", "Age", "MSI.MSS", "Stage")]
clinVU$StageT <- sub("N.*","",clinVU$TNM)
clinVU$StageM <- sub(".*M","M",clinVU$TNM)

clinVU$Date <- gsub("/",".", clinVU$Date.x)
clinVU$DateYY <- sub(".{6}","", clinVU$Date)
clinVU$DateMM <- sub("\\..*","", sub(".{3}","", clinVU$Date))
clinVU$DateDD <- sub("\\..*","", clinVU$Date)
clinVU$Date <- as.factor(paste0(clinVU$DateYY, clinVU$DateMM, clinVU$DateDD))
clinVU$Date <- factor(clinVU$Date, levels= levels(clinVU$Date),
                      labels = c("NA",levels(clinVU$Date)[2:27],"NA"),
                      ordered=T) #for checking batch effect

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
clinVU$SampleOrigin <- NULL #is most incomplete

##extract patient number from sample for samples with paired samples
clinVU$patient <- sub("_.*","", sub("\\..*","", sub("^[P,M]", "", clinVU$sampleID)))
##add the survival data for mets from the matching primaries
clinVU[grep("Met",clinVU$sampleType), "OS"] <- clinVU[grep("CRC",clinVU$sampleType), 
                                                      "OS"][match(clinVU[grep("Met",clinVU$sampleType),"patient"],
                                                    clinVU[grep("CRC",clinVU$sampleType),"patient"])]

clinVU[grep("Met",clinVU$sampleType), "Syn_metachroom"] <- clinVU[grep("CRC",clinVU$sampleType), 
                                                      "Syn_metachroom"][match(clinVU[grep("Met",clinVU$sampleType),"patient"],
                                                                  clinVU[grep("CRC",clinVU$sampleType),"patient"])]
clinVU[grep("Met",clinVU$sampleType), "Event_OS"] <- clinVU[grep("CRC",clinVU$sampleType), 
                                                                  "Event_OS"][match(clinVU[grep("Met",clinVU$sampleType),"patient"],
                                                                                          clinVU[grep("CRC",clinVU$sampleType),"patient"])]
##simplify purity estimates 
clinVU$perc <- factor(clinVU$`Tu%`, ordered=T, #make ordered
          levels=c("0" ,  "35" , "40",  "45",  "50",  "55" , "60",  "65" ,"<70" , "70",">70" ))
##simplify location in colon 
clinVU$LRcolon <- factor(clinVU$Location, labels=c("R","R","L","R", "L",
                                                   "Met.omentum","ND","L", "Rectum","L",
                                                   "T" ))
clinVU$LRcolon <- factor(clinVU$LRcolon, ordered = T,
                         levels=c("R","T", "L","Rectum",
                                  "Met.omentum","xND" ),
                         labels=c("right.colon", "transv.colon","left.colon", "rectum",
                                  "Met_perit", "ND"))
summary(clinVU$sampleType)



################ quality control EGAS1127 (VU) ################
##exploratory PCA
VU_trans <- preProcess(VU_rc_vst[grep("prim", clinVU$sampleType),
                                 which(apply(VU_rc_vst[,1:1437], 2, var) != 0 )],#[which(sourceS=="P"),], 
  method=c(#"BoxCox",
    "center","scale", 
    "pca"),
  thresh = list(thresh = 0.60))
VU_PC <- predict(VU_trans, VU_rc_vst[,1:1437])

#### plot kmeans class labels on PCA #####
VU_kM <- kmeans(VU_rc_vst[grep("prim", clinVU$sampleType),1:1437], 4, nstart = 20)
VU_kM$cluster <- as.factor(VU_kM$cluster)
summary(VU_kM$cluster)

pdf(paste0(projHome, "miRNA classifier/analyses/plots/VUdata_vst_prim_kM_PCA.pdf"),
    onefile = T)
ggplot(VU_PC,aes(x=PC1,y=PC2,
              colour= VU_kM$cluster))+#clinVU$sampleType)) +
  geom_point(na.rm = F) +
  theme_minimal() +
  # scale_color_manual(values = wes_palette("Darjeeling1", 24,
  #                                        type = c( "continuous")) )
  scale_color_npg()#lancet() #jco # npg # aaas
dev.off()

## to check for batch effects
colorBy <- clinVU$Date[match(rownames(VU_rc_vst[grep("prim", clinVU$sampleType),]), 
                             clinVU$sampleID)]

##exploratory tSNE
set.seed(9)  
tsne_model_VU <- Rtsne( as.matrix( VU_rc_vst[grep("prim", 
                                                  clinVU$sampleType),1:1437] ), 
                      check_duplicates=FALSE, 
                      pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_VU <- as.data.frame(tsne_model_VU$Y) 
ggplot(d_tsne_VU, aes(x=V1, y=V2, colour=pred_VU_RF$CMS_20[grep("prim", 
                                                             clinVU$sampleType)]))+
  geom_point() +
  xlab("") + ylab("") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme_minimal()+
  scale_color_manual(values = pal_npg("nrc")(5)[c(1,3,4,5)])
  scale_color_manual(values = c(wes_palette("Darjeeling1", 24, 
                                                    type = c( "continuous")) ))
scale_colour_manual(values = paletteCMS)
ggsave(paste0(projHome, "miRNA_classifier/analyses/plots/VUdata_vst_tSNE_prim_date.pdf"))






############ READ dataset #################

## fetch TCGA data READ:
clinical.read <- GDCquery_clinic(project = c("TCGA-READ"), type = "clinical")
CMS_samples$read <- clinical$tissue_or_organ_of_origin[match(CMS_samples$SampleId,
                                                             as.character(clinical.read$submitter_id))]

rownames(clinical.read) <- clinical.read$submitter_id
query.exp <- GDCquery(project = c("TCGA-READ"),
                      legacy = T,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      file.type = "results", # or normalized_results
                      experimental.strategy = "RNA-Seq",
                      platform = "Illumina GA"#OR: GA (for READ there seem to be duplicates with Illumina HiSeq & Illumina GA)
)
GDCdownload(query.exp)
TCGARnaseqSE <- as.data.frame(GDCprepare(query= query.exp,
                                         summarizedExperiment = F, save = F))
TCGARnaseqRSEM.GA <- TCGARnaseqSE[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_
query.exp <- GDCquery(project = c("TCGA-READ"),
                      legacy = T,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      file.type = "results", # or normalized_results
                      experimental.strategy = "RNA-Seq",
                      platform = "Illumina HiSeq"#OR: GA (for READ there seem to be duplicates with Illumina HiSeq & Illumina GA)
)
GDCdownload(query.exp)
TCGARnaseqSE <- as.data.frame(GDCprepare(query= query.exp,
                                         summarizedExperiment = F, save = F))
TCGARnaseqRSEM.HiSeq <- TCGARnaseqSE[,grep("scaled_estimate", colnames(TCGARnaseqSE)) ] #or raw_count_

unique(c(colnames(TCGARnaseqRSEM.HiSeq), colnames(TCGARnaseqRSEM.GA)))
TCGARnaseqRaw <- cbind( TCGARnaseqRSEM.HiSeq, TCGARnaseqRSEM.GA )
colnames( TCGARnaseqRaw ) <- c(paste0("HS.",colnames(TCGARnaseqRSEM.HiSeq)),
                               paste0("GA.",colnames(TCGARnaseqRSEM.GA)))
TCGARnaseqRaw <- TCGARnaseqRaw[,grep("-01A-.*", colnames(TCGARnaseqRaw))] #keep only first primary tu sample
remove(TCGARnaseqSE) #163 unique primary READ samples when combining GA and HiSeq
colnames(TCGARnaseqRaw) <- sub("scaled_estimate_","",colnames(TCGARnaseqRaw))
colnames(TCGARnaseqRaw) <- sub("-01A-.*", "",colnames(TCGARnaseqRaw))

batchInfo <- data.frame("sample"=colnames(TCGARnaseqRaw),
                        "batch"=sub("\\..*","",colnames(TCGARnaseqRaw)))
TCGARnaseqDF.READ <- limma::removeBatchEffect(as.matrix(log2(TCGARnaseqRaw+0.000001)), #the entire dataset
                                              batch=batchInfo$batch)
colnames(TCGARnaseqDF.READ) <- sub("^...","",colnames(TCGARnaseqDF.READ)) #GA./HS. being the batch identifier in the colname, keeps HiSeq source in dupliciates
TCGARnaseqDF.READ <- TCGARnaseqDF.READ[,!duplicated(colnames(TCGARnaseqDF.READ))]
rownames(TCGARnaseqDF.READ) <- oldRowNames

write.csv(data.frame("LV.table"=CMS_samples$CMS_final_netw_RF[match(rownames(CMSread.rf$predictedCMS),
                                                               CMS_samples$SampleId)],
                CMSread.rf$predictedCMS,
                CMSread.ss$predictedCMS), row.names = T,
            paste0(projHome, "miRNA_classifier/Data/TCGA-READ-allRNA-01Aonly_CMS-labels.csv"))
write.csv(TCGARnaseqDF.READ,
            paste0(projHome, "miRNA_classifier/Data/TCGA-READ-mRNA-GA-HiSeq_scaled_estimate_01Aonly.csv"))

TCGARnaseqDF.READ <- read.csv(paste0(projDir,"/Data/TCGA-READ-mRNA-GA-HiSeq_scaled_estimate_01Aonly.csv"),
                              row.names = 1, header = T)
colnames(TCGARnaseqDF.READ) <- gsub("\\.","-",colnames(TCGARnaseqDF.READ))

## get the classification from CMSclassifier package
rownames(TCGARnaseqDF.READ) <- sub(".*\\|", "", rownames(TCGARnaseqDF.READ))#keep only the ENTREZ for classififer
CMSread.rf <-classifyCMS(TCGARnaseqDF.READ, method="RF")
CMSread.ss <- classifyCMS(TCGARnaseqDF.READ, method="SSP")
CMSsymbols <- mapIds(org.Hs.eg.db, keys=listModelGenes(method = "RF"),
                     'SYMBOL', 'ENTREZID')
View(data.frame("LV.table"=CMS_samples$CMS_final_netw_RF[match(rownames(CMSread.rf$predictedCMS),
                                                               CMS_samples$SampleId)],
                CMSread.rf$predictedCMS,
                CMSread.ss$predictedCMS))
predictedCMS.READ <- read.csv(paste0(projDir,"/Data/TCGA-READ-allRNA-01Aonly_CMS-labels.csv"),
                              row.names = 1)

############ READ miRNA data ###########
## fetch TCGA READ data miRNAs:
query.exp <- GDCquery(project = c("TCGA-READ"),
                      #legacy = F,
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification",
                      experimental.strategy = "miRNA-Seq")
GDCdownload(query.exp)
TCGARnaseqSE <- GDCprepare(query= query.exp, summarizedExperiment = F, save = F)
TCGARnaseqDF <- as.data.frame(TCGARnaseqSE[,grep("read_count_", colnames(TCGARnaseqSE)) ]) #or normalized_count_ for normalized_results
colnames(TCGARnaseqDF) <- sub("-01A-.*","",sub("read_count_","",colnames(TCGARnaseqDF)))
rownames(TCGARnaseqDF) <- TCGARnaseqSE$miRNA_ID
remove(TCGARnaseqSE)
TCGARnaseqDF$miRcomb <- sub("-[1-9]$","",row.names(TCGARnaseqDF))
miR_READ_comb <- as.data.frame(TCGARnaseqDF %>% dplyr::group_by(miRcomb) %>%
                                 dplyr::summarise_if(.predicate = function(x) is.numeric(x),
                                                     .funs = mean)) #affects around 180miRNA
row.names(miR_READ_comb) <- miR_READ_comb$miRcomb #afterwards remove column 1

##normalize, in right direction for RF
miR_READ_vst <- as.data.frame( t( varianceStabilizingTransformation( as.matrix(
  round( miR_READ_comb[,-1], 0 ) ) ) ) )
colnames(miR_READ_vst) <- gsub("-", ".", colnames(miR_READ_vst))
#get the CMS labels from CMSclassifier after vst
miR_READ_vst <- miR_READ_vst[-grep("-.+-.+-.+-.+", rownames(miR_READ_vst)),] #get rid of normal and extra samples

rownames(miR_READ_vst) <- sub("-01A-.*","",rownames(miR_READ_vst))
miR_READ_vst$CMS.cl.rf<- predictedCMS.READ[match(rownames(miR_READ_vst),
                        rownames(predictedCMS.READ)),"RF"]
miR_READ_CMS.info <- as.character(miR_READ_vst$CMS.cl.rf)
miR_READ_CMS.info[which(is.na(miR_READ_CMS.info))] <- "NA"

##batch effect removal
batchInfo <- data.frame("sample"=rownames(miR_READ_vst),
                        "batch"=factor(sub("-.*","",
                                           sub("TCGA-","",
                                               rownames(miR_READ_vst)))))#[
                                                 #grep("CMS", miR_READ_vst$CMS.cl.rf),]
designR <- model.matrix( ~ 0 + miR_READ_CMS.info)#[
  #grep("CMS", miR_READ_vst$CMS.cl.rf)] ) ## do i need an intercept? otherwise it compares all to CMS1 or CMS2 if releveled?
miR_READ_vst_BR <- as.data.frame(t(limma::removeBatchEffect(as.matrix(t(
  miR_READ_vst[#grep("CMS", miR_READ_vst$CMS.cl.rf),
               grep("hsa", colnames(miR_READ_vst))])), #the entire dataset
                                              design=designR,
                                              batch=batchInfo$batch)))
purity.ABS.read <- purity.ABS[match(row.names(miR_READ_vst_BR), purity.ABS$TCGA.patID), ]


##get the CMS labels from Louis table
miR_READ_vst_BR$CMS.lv <- CMS_samples$CMS[match(rownames(miR_READ_vst_BR),
                                                CMS_samples$SampleId)]#get the CMS labels from Louis table
miR_READ_vst_BR$CMS.cl.rf <- miR_READ_vst$CMS.cl.rf#[grep("CMS", miR_READ_vst$CMS.cl.rf)]

write.csv(miR_READ_vst,
           paste0(projHome, "miRNA_classifier/Data/TCGA-miR_READ_vst_CMS-labels.csv"))
write.csv(miR_READ_vst_BR,
           paste0(projHome, "miRNA_classifier/Data/TCGA-miR_READ_vst_BR_CMS-labels.csv"))
miR_READ_vst_BR <- read.csv(file=paste0(projDir,"/Data/TCGA-miR_READ_vst_BR_CMS-labels.csv"),
          row.names = 1)


##check batch effects
set.seed(5678)
tsne_model_READ <- Rtsne( as.matrix( miR_READ_vst_BR[,1:156]), 
                        check_duplicates=FALSE, 
                        pca=TRUE, perplexity=20, theta=0.5, dims=2)
colREAD <- sub("-.*","",sub("TCGA-","",
                            clinical.read$submitter_id[match(rownames(miR_READ_vst), 
                                 clinical.read$submitter_id)]))
ggplot(as.data.frame(tsne_model_READ$Y), 
       aes(x=V1, y=V2, color = miR_READ_vst_BR$CMS.lv))+
  geom_point() +
  xlab("") + ylab("") +
  theme_minimal()+
  # scale_colour_manual(values = paletteCMS, name = "CMS")
  # 
  scale_color_manual(values = c(wes_palette("Darjeeling1", 24, 
                                            type = c( "continuous")) ))
ggsave(paste0(plotDir,"/READ-tSNE_miR_vst_BR_CMS.lv.pdf" ))
##-> removed the AG difference while retaining mRNA-based predicted CMS differences


batchInfo <- data.frame("sample"=rownames(miR_READ_vst_BR),
                        "batch"=factor(sub("-.*","",
                                           sub("TCGA-","",
                                               rownames(miR_READ_vst_BR)))))#[
                                                 #grep("CMS", miR_READ_vst$CMS.cl.rf),]

