######## RTN network ########
library(RTN)
library(DESeq2)
library(RedeR)

#### top 20 miR from RF (50x, same for 80x) ####
tfs <- c("hsa-mir-625", "hsa-mir-592", "hsa-mir-552", "hsa-mir-218" ,
         "hsa-mir-31", "hsa-mir-375", "hsa-mir-143",  "hsa-mir-615",  
         "hsa-mir-335",  "hsa-mir-146b", "hsa-mir-99a", "hsa-mir-141", 
         "hsa-mir-92b" ,  "hsa-mir-942", "hsa-mir-3170",  "hsa-mir-362",  
         "hsa-mir-30a",   "hsa-mir-155", "hsa-mir-582", "hsa-mir-92a" )
#### just in case WNT genes need anotation ####
markerGenesFlierWNTup_AdCa <- c("Cemip","(Kiaa1199)", "Foxq1", "Cldn2", "Etv4", "Phlda1", "Hes6", "Cd44", "Slc7a5", "Sord", "Rnf43", "(Flj20315)", "Ctps", "Enc1", "Tead4", "Met", "Myc", "Rrp12", "(Kiaa0690)", "C19orf48","(Mgc13170)", "Tspan5","(Tm4sf9)", "Cdc25a", "Cdkn3", "Bysl", "Scd", "Sox9","Pdcd2l", "(Mgc13096)", "Sox4", "Wdr77", "Trmt1", "Nob1p", "Psf1", "Nle1", "Pmscl1", "Eif5a2", "Tgif", "Ppil1", "Ets2", "Slc39a10", "Rrp9", "(Rnu3ip2)", "Utp4","(Cirh1a)", "Znf511", "Ddx21", "Dcun1d5", "Hilpda", "Znf593", "Mrto4", "Dkc1", "Odc1", "Cdk4","Psmg4", "(Loc389362)", "Ppan", "Mac30","Cytor", "(Mgc4677)", "Nop16", "Aurkb", "C12orf45","(Mgc40397)", "Cd320", "Znf703", "Wdr12", "Ns", "Cdt1", "Pno1","(Loc56902)", "Sdccag16", "Noc3l", "Peo1", "Cad", "Xpot", "Lyar","(Flj20425)", "Polr1c", "Nat10","(Flj10774)", "Cdca7", "Wdr74", "Yap1", "Mtvr1", "Srm", "Nxt1", "Snhg15","(Loc285958)", "Ece2","(Mgc2408)", "Mrpl36", "Ddx10", "Mrps12", "Rabepk", "(Rab9p40)", "Rheb", "Wdr3", "Rfp", "Polr1b", "Gpatch4", "(Flj20249)", "Ddx20", "Nifk","(Mki67ip)", "Farsa", "Crtap", "Sf3a2", "Gemin5", "Polr1d", "Mettl1", "Rrp1b","(Kiaa0179)", "Znf259", "Rrp41", "Hnrpdl", "Nprl5","(C16orf35)", "Ccdc86", "(Mgc2574)", "Ddx56", "Trex2", "Nola1","Rrp7a", "(Cgi-96)", "Imp4", "Mrps30", "Yrdc", "(Flj23476)", "Heatr1", "Slc29a2", "Tere1", "Mrps26", "Slc19a1", "Ydjc","(Loc150223)", "Tp53rk", "(C20orf64)", "Cog8","Fam60a", "(C12orf14)", "Mybbp1a", "Pdcd11", "Nop2", "Nol5a", "Slc3a2", "Nudcd3", "(Kiaa1068)")
markerGenesFlierWNTup_onlyCa <- c("CXCL5", "ZIC2", "CBX2",  "C20orf82", "GLS2", "BAG2", "ZNF503", "SNHG17","(LOC388796)", "LEF1", "FLJ12683", "", "NELF", "HOXA11", "CKLFSF7", "LIMS1", "SLC19A2", "SCLY", "HIRA", "SLC52A2","(GPR172A)", "DOCK4", "AMPD2", "MAT2A", "C1orf109","(FLJ20508)", "TPST2", "GASK1B", "(DKFZp434L142)", "PPIF", "DUS3L","(LOC56931)", "ISG20L2","(FLJ12671)", "SLC25A19", "SLC1A5", "TIMM10", "POLRMT", "ITGB4BP", "ZNF275", "SFRS7", "RBM12", "UMPK")
markerGenesFlierWNTup_onlyAd <- c("Ascl2", "Axin2", "Gpr49", "Tdgf1", "Tbx3", "Ephb3", "Ephb2", "Znrf3", "Ucc1", "Il20ra", "Rhobtb3", "Slc5a1", "STAMBPL1","(Amsh-Lp)", "Gemin4", "Gpx2", "Id3", "Skb1", "Kcne3", "CEP76","(C18orf9)","RRP15", "(Cgi-115)", "ZNHIT2", "(C11orf5)", "Nse1", "Fut4", "Dhx33", "KIZ", "(C20orf19)", "Gemin6", "Jag1", "Pigw", "Myb", "Apex1", "SPIN4", "(Loc139886)", "NIP7", "(Cgi-37)", "Hsa9761", "SYBU","(Flj20366)", "Ptma", "Mgc16824", "Hspc242", "Znf278",  "PAAF1", "(Flj11848)", "Srprb", "Znrd1", "Tfb2m", "Rpc62", "Arid5b", "Hnrpa0", "Mettl2", "Jtv1", "Loc20525", "Mrpl15", "Znf339", "Mrpl12")
wntGenes <- toupper( c(markerGenesFlierWNTup_AdCa, markerGenesFlierWNTup_onlyCa, markerGenesFlierWNTup_onlyAd))

### get the data ####
TCGA.COAD.mRna.scale <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-mRna-scale.csv", 
                                 row.names=1)

TCGA.COAD.miR.scale <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/Data/TCGA-COAD-miR-scale.csv", 
                                row.names=1)
rownames(TCGA.COAD.miR.scale) <- gsub("-",".", rownames(TCGA.COAD.miR.scale))


DEgenes <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/DE-results/resSigTable-COAD-GA-mRNA-LV-CMS-CMS1.csv",
                    as.is = T)

DEgenes <- DEgenes[which( DEgenes$GeneID %in% rownames(TCGA.COAD.mRna.scale)), ]


DEmiRS <- read.csv("/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/DE-results/resSigTable-p0.05fc0.71-COAD-GA-miRs-LV-CMS-CMS3.csv",
                   as.is = T)

DEmiRS$GeneID <- gsub("-",".", DEmiRS$GeneID)
DEmiRS <- DEmiRS[which( DEmiRS$GeneID %in% rownames(TCGA.COAD.miR.scale)), ]

##identified 200 top DE genes for each CMS and |log2FC|>0.85, FDR<0.001
##from TCGA-COAD-GA
#### merge data according to RTN vignette ####
expData <- rbind( TCGA.COAD.miR.scale[ as.character( DEmiRS$GeneID ) ,
                                       colnames(TCGA.COAD.miR.scale)],
                  TCGA.COAD.mRna.scale[ as.character( DEgenes$GeneID ) ,
                                        colnames(TCGA.COAD.miR.scale)])
rowAnnotation <- data.frame("PROBEID"=rownames(expData), 
                            "SYMBOL"=sub("\\|.*","", rownames(expData)),
                            "WNT"=as.numeric( sub("\\|.*","", rownames(expData)) %in% wntGenes ),
                            "L2FC"=round( c( DEmiRS$log2FoldChange, DEgenes$log2FoldChange), 1),
                            "NRF"=round( importances_df[match(rownames(expData), 
                                                              rownames( importances_df)),"Overall"]) )
rowAnnotation$NRF[is.na(rowAnnotation$NRF)] <- 0
colAnnot <- cbind( data.frame("ID"=rownames(cmsAnnot[colnames(TCGA.COAD.miR.scale),]) ),
                   cmsAnnot[colnames(TCGA.COAD.miR.scale),])
colAnnot$LV.table[is.na(colAnnot$LV.table)] <- "NOLBL"


#### RTN analysis ####
rtni <- tni.constructor(expData = as.matrix(expData), 
                        regulatoryElements = DEmiRS$GeneID[which(DEmiRS$pval<0.001 & #this is according to Fessler paper
                                                                   DEmiRS$log2FoldChange > 0)], #check the most DE miRs
                        ## also plot upregulated miRs but separately from downregualted miRs  
                        rowAnnotation = rowAnnotation, 
                        colAnnotation = colAnnot)
rtni <- tni.permutation(rtni)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni, eps = NA) 

#regulons <- tni.get(rtni, what = "regulons.and.mode", idkey = "SYMBOL")
# write.csv(data.frame(names(regulons$hsa.mir.141), regulons$hsa.mir.141 ), "/Users/ronjaadam/projects/miRNA_mCRC/miRNA classifier/analyses/Network-RTN/TCGA-COAD-GA-CMS-DEgenes-all-CMS4-miR141-regulon.txt",
#           row.names = F, quote = F)

#### RTN graph via RedeR ####
g <- tni.graph(rtni) #default parameters got CMS4 down close to Fessler Figure 
g <- att.setv(g, from = "L2FC", to = "nodeColor", breaks=seq(-5,5,0.8), pal=2)
g <- att.setv(g, from = "NRF", to = "nodeSize")
rdp <- RedPort()
calld(rdp, checkcalls=TRUE)
resetd(rdp)
addGraph(rdp, g, layout=NULL)
addLegend.color(rdp, g, type="edge")
addLegend.color(rdp, g, type="node")
addLegend.size(rdp, g, type="node")
addLegend.shape(rdp, g)
relax(rdp, ps = TRUE)
