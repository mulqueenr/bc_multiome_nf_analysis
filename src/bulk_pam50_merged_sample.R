
"""

######################
For Bulk
######################

"""

```R genefu
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(genefu)
library(dplyr)
library(org.Hs.eg.db)
args = commandArgs(trailingOnly=TRUE)

dat<-readRDS(args[1])

#Using genefu per pseudobulked sample
#data: Matrix of annotations with at least one column named "EntrezGene.ID"
#'   (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust
#'   model), dimnames being properly defined.
#do.mapping TRUE if the mapping through Entrez Gene ids must be performed
#'   (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.

#CDCA1 KNTC2 ORC6L use different names in our data
#NUF2, NDC80, ORC6 resp.
pam50_genes<-c('ACTR3B', 'ANLN', 'BAG1', 'BCL2', 'BIRC5', 'BLVRA', 
  'CCNB1', 'CCNE1', 'CDC20', 'CDC6', 'NUF2', 'CDH3', 'CENPF', 'CEP55', 
  'CXXC5', 'EGFR', 'ERBB2', 'ESR1', 'EXO1', 'FGFR4', 'FOXA1', 'FOXC1', 
  'GPR160', 'GRB7', 'KIF2C', 'NDC80', 'KRT14', 'KRT17', 'KRT5', 'MAPT', 
  'MDM2', 'MELK', 'MIA', 'MKI67', 'MLPH', 'MMP11', 'MYBL2', 'MYC', 'NAT1', 
  'ORC6', 'PGR', 'PHGDH', 'PTTG1', 'RRM2', 'SFRP1', 'SLC39A6', 'TMEM45B', 
  'TYMS', 'UBE2C', 'UBE2T')

sample_names<-paste(unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(1))),
  unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(2))),sep="_") #this might not be needed anymore

counts<-as.data.frame(t(dat[["RNA"]]@counts)) 
counts<-cbind(counts,sample_names)
counts<-as.data.frame(counts %>% group_by(sample_names) %>% summarize_all(funs(sum)))
row.names(counts)<-counts$sample_name
counts<-counts[,2:ncol(counts)]
counts<-counts[,colSums(counts)>0]
dat_in<-counts

dat_in<-dat_in[!(row.names(dat_in) %in% c("RM_4","sample_15","sample_19")),] #exclude NAT samples
dat_in<-NormalizeData(dat_in,normalization.method="CLR")
dannot<-as.data.frame(cbind(Gene.Symbol=colnames(dat_in),EntrezGene.ID=mapIds(org.Hs.eg.db, colnames(dat_in), 'ENTREZID', 'SYMBOL'),probe=colnames(dat_in)))
pam50_out<-molecular.subtyping(sbt.model="pam50",data=dat_in,annot=dannot,do.mapping=TRUE,verbose=T)

#try this as well
#pam50_out_model<-intrinsic.cluster(data=dat_in,annot=dannot,do.mapping=TRUE,std="robust",intrinsicg=pam50$centroids.map[,c("probe","EntrezGene.ID")],verbose=T,mins=0)#,mapping=dannot)
#pam50_out<-intrinsic.cluster.predict(sbt.model=pam50_out_model$model, data=dat_in, annot=dannot, do.mapping=TRUE,do.prediction.strength=TRUE,verbose=TRUE)
#saveRDS(pam50_out,file="pseudobulk_pam50.rds")

pam50_meta<-setNames(nm=row.names(dat@meta.data),pam50_out$subtype[match(dat$sample, names(pam50_out$subtype))])
dat<-AddMetaData(dat,pam50_meta,col.name="pseudobulk_genefu_pam50")
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

#tried just epithelial, tried both old method (intrinsic cluster) and updated method (molecular subtyping). maybe play around with normalizing first?
#limit to epithelial? or maybe read up on proper normalization? our HER2+ isn't being labelled as such

```

```R ssbc
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(genefu)
library(dplyr)
library(org.Hs.eg.db)
library(sspbc)

args = commandArgs(trailingOnly=TRUE)

dat<-readRDS(args[1])

sample_names<-paste(unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(1))),
  unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(2))),sep="_")
counts<-as.data.frame(t(dat[["RNA"]]@counts)) 
counts<-cbind(counts,sample_names)
counts<-as.data.frame(counts %>% group_by(sample_names) %>% summarize_all(funs(sum)))
row.names(counts)<-counts$sample_name
counts<-counts[,2:ncol(counts)]
counts<-counts[,colSums(counts)>0]
dat_in<-counts
dat_in<-dat_in[!(row.names(dat_in) %in% c("RM_4","sample_15","sample_19")),] #exclude NAT samples
dat_in<-as.data.frame(t(dat_in))

#set up matrix by unique entrez gene names
dat_in<-dat_in[!duplicated(mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')),]
dat_in<-dat_in[!isNA(mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')),]
row.names(dat_in)<-mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')
myresults <- applySSP(gex=as.matrix(dat_in), id=row.names(dat_in), ssp.name="ssp.pam50",id.type="EntrezGene",report=TRUE)

dat_pam50<-setNames(nm=row.names(dat@meta.data),myresults[match(dat@meta.data$sample,row.names(myresults)),1])
dat<-AddMetaData(dat,dat_pam50,col.name="pseudobulk_sspbc_PAM50")
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```


```R
library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(genefu)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
dat<-readRDS("phase2.QC.filt.SeuratObject.rds")

#cell lineage
  load("/home/groups/CEDAR/mulqueen/ref/embo/Human-PosSigGenes.RData")
  ls()
  #[1] "Basal" "LP"    "ML"    "Str" #lineage types
  lineage_in=list(EMBO_Basal=Basal,EMBO_LP=LP,EMBO_ML=ML,EMBO_Str=Str)

#Immune markers
  immune_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/ImmuneMarkers2.txt",header=T,sep="\t")
  immune_in<-lapply(split(immune_in,immune_in$CellType),function(x) x$Signatures)#split up data frame to a named list of genes per cell type
  names(immune_in)<-paste0("EMBO_",names(immune_in))#rename the list just so we can track the source

#using both the given PAM50 short list and the Swarbrick supplied more extensive gene list below
  PAM50_in<-read.table("/home/groups/CEDAR/mulqueen/ref/embo/PAM50.txt",header=T,sep="\t")
  PAM50_in<-lapply(split(PAM50_in,PAM50_in$Subtype),function(x) x$Gene)#split up data frame to a named list of genes per cell type
  names(PAM50_in)<-paste0("PAM50_",names(PAM50_in))
  features_in=c(immune_in,PAM50_in)   

molecular.subtyping(sbt.model="pam50")
#SCSubtype Features determined by Swarbrick manuscript (Supp Table 4)
  module_feats<-list()
  module_feats[["Basal_SC"]]=c('EMP1', 'TAGLN', 'TTYH1', 'RTN4', 'TK1', 'BUB3', 'IGLV3.25', 'FAM3C', 'TMEM123', 'KDM5B', 'KRT14', 'ALG3', 'KLK6', 'EEF2', 'NSMCE4A', 'LYST', 'DEDD', 'HLA.DRA', 'PAPOLA', 'SOX4', 'ACTR3B', 'EIF3D', 'CACYBP', 'RARRES1', 'STRA13', 'MFGE8', 'FRZB', 'SDHD', 'UCHL1', 'TMEM176A', 'CAV2', 'MARCO', 'P4HB', 'CHI3L2', 'APOE', 'ATP1B1', 'C6orf15', 'KRT6B', 'TAF1D', 'ACTA2', 'LY6D', 'SAA2', 'CYP27A1', 'DLK1', 'IGKV1.5', 'CENPW', 'RAB18', 'TNFRSF11B', 'VPS28', 'HULC', 'KRT16', 'CDKN2A', 'AHNAK2', 'SEC22B', 'CDC42EP1', 'HMGA1', 'CAV1', 'BAMBI', 'TOMM22', 'ATP6V0E2', 'MTCH2', 'PRSS21', 'HDAC2', 'ZG16B', 'GAL', 'SCGB1D2', 'S100A2', 'GSPT1', 'ARPC1B', 'NIT1', 'NEAT1', 'DSC2', 'RP1.60O19.1', 'MAL2', 'TMEM176B', 'CYP1B1', 'EIF3L', 'FKBP4', 'WFDC2', 'SAA1', 'CXCL17', 'PFDN2', 'UCP2', 'RAB11B', 'FDCSP', 'HLA.DPB1', 'PCSK1N', 'C4orf48', 'CTSC')
  module_feats[["Her2E_SC"]]=c('PSMA2', 'PPP1R1B', 'SYNGR2', 'CNPY2', 'LGALS7B', 'CYBA', 'FTH1', 'MSL1', 'IGKV3.15', 'STARD3', 'HPD', 'HMGCS2', 'ID3', 'NDUFB8', 'COTL1', 'AIM1', 'MED24', 'CEACAM6', 'FABP7', 'CRABP2', 'NR4A2', 'COX14', 'ACADM', 'PKM', 'ECH1', 'C17orf89', 'NGRN', 'ATG5', 'SNHG25', 'ETFB', 'EGLN3', 'CSNK2B', 'RHOC', 'PSENEN', 'CDK12', 'ATP5I', 'ENTHD2', 'QRSL1', 'S100A7', 'TPM1', 'ATP5C1', 'HIST1H1E', 'LGALS1', 'GRB7', 'AQP3', 'ALDH2', 'EIF3E', 'ERBB2', 'LCN2', 'SLC38A10', 'TXN', 'DBI', 'RP11.206M11.7', 'TUBB', 'CRYAB', 'CD9', 'PDSS2', 'XIST', 'MED1', 'C6orf203', 'PSMD3', 'TMC5', 'UQCRQ', 'EFHD1', 'BCAM', 'GPX1', 'EPHX1', 'AREG', 'CDK2AP2', 'SPINK8', 'PGAP3', 'NFIC', 'THRSP', 'LDHB', 'MT1X', 'HIST1H4C', 'LRRC26', 'SLC16A3', 'BACE2', 'MIEN1', 'AR', 'CRIP2', 'NME1', 'DEGS2', 'CASC3', 'FOLR1', 'SIVA1', 'SLC25A39', 'IGHG1', 'ORMDL3', 'KRT81', 'SCGB2B2', 'LINC01285', 'CXCL8', 'KRT15', 'RSU1', 'ZFP36L2', 'DKK1', 'TMED10', 'IRX3', 'S100A9', 'YWHAZ')
  module_feats[["LumA_SC"]]=c('SH3BGRL', 'HSPB1', 'PHGR1', 'SOX9', 'CEBPD', 'CITED2', 'TM4SF1', 'S100P', 'KCNK6', 'AGR3', 'MPC2', 'CXCL13', 'RNASET2', 'DDIT4', 'SCUBE2', 'KRT8', 'MZT2B', 'IFI6', 'RPS26', 'TAGLN2', 'SPTSSA', 'ZFP36L1', 'MGP', 'KDELR2', 'PPDPF', 'AZGP1', 'AP000769.1', 'MYBPC1', 'S100A1', 'TFPI2', 'JUN', 'SLC25A6', 'HSP90AB1', 'ARF5', 'PMAIP1', 'TNFRSF12A', 'FXYD3', 'RASD1', 'PYCARD', 'PYDC1', 'PHLDA2', 'BZW2', 'HOXA9', 'XBP1', 'AGR2', 'HSP90AA1') 
  module_feats[["LumB_SC"]]=c('UGCG', 'ARMT1', 'ISOC1', 'GDF15', 'ZFP36', 'PSMC5', 'DDX5', 'TMEM150C', 'NBEAL1', 'CLEC3A', 'GADD45G', 'MARCKS', 'FHL2', 'CCDC117', 'LY6E', 'GJA1', 'PSAP', 'TAF7', 'PIP', 'HSPA2', 'DSCAM.AS1', 'PSMB7', 'STARD10', 'ATF3', 'WBP11', 'MALAT1', 'C6orf48', 'HLA.DRB1', 'HIST1H2BD', 'CCND1', 'STC2', 'NR4A1', 'NPY1R', 'FOS', 'ZFAND2A', 'CFL1', 'RHOB', 'LMNA', 'SLC40A1', 'CYB5A', 'SRSF5', 'SEC61G', 'CTSD', 'DNAJC12', 'IFITM1', 'MAGED2', 'RBP1', 'TFF1', 'APLP2', 'TFF3', 'TRH', 'NUPR1', 'EMC3', 'TXNIP', 'ARPC4', 'KCNE4', 'ANPEP', 'MGST1', 'TOB1', 'ADIRF', 'TUBA1B', 'MYEOV2', 'MLLT4', 'DHRS2', 'IFITM2')
  module_feats[["proliferation_score"]]<-c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")

#Swarbrick Gene Module Classification (Supp Table 5)
gene_module<-list()
  gene_module[["gene_module_1"]]<-c('ATF3', 'JUN', 'NR4A1', 'IER2', 'DUSP1', 'ZFP36', 'JUNB', 'FOS', 'FOSB', 'PPP1R15A', 'KLF6', 'DNAJB1', 'EGR1', 'BTG2', 'HSPA1B', 'HSPA1A', 'RHOB', 'CLDN4', 'MAFF', 'GADD45B', 'IRF1', 'EFNA1', 'SERTAD1', 'TSC22D1', 'CEBPD', 'CCNL1', 'TRIB1', 'MYC', 'ELF3', 'LMNA', 'NFKBIA', 'TOB1', 'HSPB1', 'BRD2', 'MCL1', 'PNRC1', 'IER3', 'KLF4', 'ZFP36L2', 'SAT1', 'ZFP36L1', 'DNAJB4', 'PHLDA2', 'NEAT1', 'MAP3K8', 'GPRC5A', 'RASD1', 'NFKBIZ', 'CTD-3252C9.4', 'BAMBI', 'RND1', 'HES1', 'PIM3', 'SQSTM1', 'HSPH1', 'ZFAND5', 'AREG', 'CD55', 'CDKN1A', 'UBC', 'CLDN3', 'DDIT3', 'BHLHE40', 'BTG1', 'ANKRD37', 'SOCS3', 'NAMPT', 'SOX4', 'LDLR', 'TIPARP', 'TM4SF1', 'CSRNP1', 'GDF15', 'ZFAND2A', 'NR4A2', 'ERRFI1', 'RAB11FIP1', 'TRAF4', 'MYADM', 'ZC3H12A', 'HERPUD1', 'CKS2', 'BAG3', 'TGIF1', 'ID3', 'JUND', 'PMAIP1', 'TACSTD2', 'ETS2', 'DNAJA1', 'PDLIM3', 'KLF10', 'CYR61', 'MXD1', 'TNFAIP3', 'NCOA7', 'OVOL1', 'TSC22D3', 'HSP90AA1', 'HSPA6', 'C15orf48', 'RHOV', 'DUSP4', 'B4GALT1', 'SDC4', 'C8orf4', 'DNAJB6', 'ICAM1', 'DNAJA4', 'MRPL18', 'GRB7', 'HNRNPA0', 'BCL3', 'DUSP10', 'EDN1', 'FHL2', 'CXCL2', 'TNFRSF12A', 'S100P', 'HSPB8', 'INSIG1', 'PLK3', 'EZR', 'IGFBP5', 'SLC38A2', 'DNAJB9', 'H3F3B', 'TPM4', 'TNFSF10', 'RSRP1', 'ARL5B', 'ATP1B1', 'HSPA8', 'IER5', 'SCGB2A1', 'YPEL2', 'TMC5', 'FBXO32', 'MAP1LC3B', 'MIDN', 'GADD45G', 'VMP1', 'HSPA5', 'SCGB2A2', 'TUBA1A', 'WEE1', 'PDK4', 'STAT3', 'PERP', 'RBBP6', 'KCNQ1OT1', 'OSER1', 'SERP1', 'UBE2B', 'HSPE1', 'SOX9', 'MLF1', 'UBB', 'MDK', 'YPEL5', 'HMGCS1', 'PTP4A1', 'WSB1', 'CEBPB', 'EIF4A2', 'S100A10', 'ELMSAN1', 'ISG15', 'CCNI', 'CLU', 'TIMP3', 'ARL4A', 'SERPINH1', 'SCGB1D2', 'UGDH', 'FUS', 'BAG1', 'IFRD1', 'TFF1', 'SERTAD3', 'IGFBP4', 'TPM1', 'PKIB', 'MALAT1', 'XBP1', 'HEBP2', 'GEM', 'EGR2', 'ID2', 'EGR3', 'HSPD1', 'GLUL', 'DDIT4', 'CDC42EP1', 'RBM39', 'MT-ND5', 'CSNK1A1', 'SLC25A25', 'PEG10', 'DEDD2')

gene_module[["gene_module_2"]]<-c('AZGP1', 'ATP5C1', 'ATP5F1', 'NHP2', 'MGP', 'RPN2', 'C14orf2', 'NQO1', 'REEP5', 'SSR2', 'NDUFA8', 'ATP5E', 'SH3BGRL', 'PIP', 'PRDX2', 'RAB25', 'EIF3L', 'PRDX1', 'USMG5', 'DAD1', 'SEC61G', 'CCT3', 'NDUFA4', 'APOD', 'CHCHD10', 'DDIT4', 'MRPL24', 'NME1', 'DCXR', 'NDUFAB1', 'ATP5A1', 'ATP5B', 'ATOX1', 'SLC50A1', 'POLR2I', 'TIMM8B', 'VPS29', 'TIMP1', 'AHCY', 'PRDX3', 'RBM3', 'GSTM3', 'ABRACL', 'RBX1', 'PAFAH1B3', 'AP1S1', 'RPL34', 'ATPIF1', 'PGD', 'CANX', 'SELENBP1', 'ATP5J', 'PSME2', 'PSME1', 'SDHC', 'AKR1A1', 'GSTP1', 'RARRES3', 'ISCU', 'NPM1', 'SPDEF', 'BLVRB', 'NDUFB3', 'RPL36A', 'MDH1', 'MYEOV2', 'MAGED2', 'CRIP2', 'SEC11C', 'CD151', 'COPE', 'PFN2', 'ALDH2', 'SNRPD2', 'TSTD1', 'RPL13A', 'HIGD2A', 'NDUFC1', 'PYCARD', 'FIS1', 'ITM2B', 'PSMB3', 'G6PD', 'CST3', 'SH3BGRL3', 'TAGLN2', 'NDUFA1', 'TMEM183A', 'S100A10', 'NGFRAP1', 'DEGS2', 'ARPC5', 'TM7SF2', 'RPS10', 'LAMTOR5', 'TMEM256', 'UQCRB', 'TMEM141', 'KRTCAP2', 'HM13', 'NDUFS6', 'PARK7', 'PSMD4', 'NDUFB11', 'TOMM7', 'EIF6', 'UQCRHL', 'ADI1', 'VDAC1', 'C9orf16', 'ETFA', 'LSM3', 'UQCRH', 'CYB5A', 'SNRPE', 'BSG', 'SSR3', 'DPM3', 'LAMTOR4', 'RPS11', 'FAM195A', 'TMEM261', 'ATP5I', 'EIF5A', 'PIN4', 'ATXN10', 'ATP5G3', 'ARPC3', 'UBA52', 'BEX4', 'ROMO1', 'SLC25A6', 'SDCBP', 'EIF4EBP1', 'PFDN6', 'PSMA3', 'RNF7', 'SPCS2', 'CYSTM1', 'CAPG', 'CD9', 'GRHPR', 'SEPP1', 'ESF1', 'TFF3', 'ARPC1B', 'ANXA5', 'WDR83OS', 'LYPLA1', 'COMT', 'MDH2', 'DNPH1', 'RAB13', 'EIF3K', 'PTGR1', 'LGALS3', 'TPI1', 'COPZ1', 'LDHA', 'PSMD8', 'EIF2S3', 'NME3', 'EIF3E', 'MRPL13', 'ZFAND6', 'FAM162A', 'ATP6V0E1', 'TMED10', 'HNRNPA3', 'PPA1', 'SNX17', 'APOA1BP', 'TUFM', 'ECHS1', 'GLTSCR2', 'RPS27L', 'NDUFB1', 'SSBP1', 'PRDX6', 'ENO1', 'PPP4C', 'COA3', 'TCEAL4', 'MRPL54', 'LAMTOR2', 'PAIP2', 'DAP', 'RPL22L1', 'C6orf203', 'TECR', 'PEBP1', 'TMED9', 'ATP6V1F', 'ESD', 'EIF3I', 'SCO2', 'ATP5D', 'UAP1', 'TMEM258', 'COX17')

gene_module[["gene_module_3"]]<-c('HLA-B', 'HLA-A', 'VIM', 'CD74', 'SRGN', 'HLA-C', 'IFI27', 'HLA-E', 'IFITM1', 'PSMB9', 'RGCC', 'S100A4', 'HLA-DRA', 'ISG15', 'IL32', 'SPARC', 'TAGLN', 'IFITM3', 'IFITM2', 'IGFBP7', 'CALD1', 'HLA-DPB1', 'HLA-DPA1', 'B2M', 'TIMP1', 'RGS1', 'FN1', 'ACTA2', 'HLA-DRB1', 'SERPING1', 'ANXA1', 'TPM2', 'TMSB4X', 'CD69', 'CCL4', 'LAPTM5', 'GSN', 'APOE', 'STAT1', 'SPARCL1', 'IFI6', 'DUSP1', 'CXCR4', 'CCL5', 'UBE2L6', 'MYL9', 'SLC2A3', 'BST2', 'CAV1', 'CD52', 'ZFP36L2', 'HLA-DQB1', 'PDLIM1', 'TNFAIP3', 'CORO1A', 'RARRES3', 'TYMP', 'C1S', 'PTRF', 'PSME2', 'CYTIP', 'COL1A1', 'PSMB8', 'NNMT', 'HLA-DQA1', 'DUSP2', 'COL1A2', 'ARHGDIB', 'COL6A2', 'FOS', 'CCL2', 'BGN', 'ID3', 'TUBA1A', 'RAC2', 'LBH', 'HLA-DRB5', 'FCER1G', 'GBP1', 'C1QA', 'COTL1', 'LUM', 'MYL6', 'GBP2', 'BTG1', 'CD37', 'HCST', 'LIMD2', 'IFIT3', 'IL7R', 'PTPRC', 'NKG7', 'FYB', 'TAP1', 'LTB', 'S100A6', 'COL3A1', 'EMP3', 'A2M', 'JUNB', 'TPM1', 'FABP4', 'TXNIP', 'SAT1', 'FXYD5', 'CD3E', 'HLA-DMA', 'CTSC', 'TSC22D3', 'MYL12A', 'CST3', 'CNN2', 'PHLDA1', 'LYZ', 'IFI44L', 'MARCKS', 'ID1', 'DCN', 'TGFBI', 'BIRC3', 'THY1', 'LGALS1', 'GPX1', 'C1QB', 'CD2', 'CST7', 'COL6A3', 'ACAP1', 'IFI16', 'ITM2B', 'POSTN', 'LDHB', 'FLNA', 'FILIP1L', 'CDKN1A', 'IRF1', 'LGALS3', 'SERPINH1', 'EFEMP1', 'PSME1', 'SH3BGRL3', 'IL2RG', 'CD3D', 'SFRP2', 'TIMP3', 'ALOX5AP', 'GMFG', 'CYBA', 'TAGLN2', 'LAP3', 'RGS2', 'CLEC2B', 'TRBC2', 'NR4A2', 'S100A8', 'PSMB10', 'OPTN', 'CTSB', 'FTL', 'KRT17', 'AREG', 'MYH9', 'MMP7', 'COL6A1', 'GZMA', 'RNASE1', 'PCOLCE', 'PTN', 'PYCARD', 'ARPC2', 'SGK1', 'COL18A1', 'GSTP1', 'NPC2', 'SOD3', 'MFGE8', 'COL4A1', 'ADIRF', 'HLA-F', 'CD7', 'APOC1', 'TYROBP', 'C1QC', 'TAPBP', 'STK4', 'RHOH', 'RNF213', 'SOD2', 'TPM4', 'CALM1', 'CTGF', 'PNRC1', 'CD27', 'CD3G', 'PRKCDBP', 'PARP14', 'IGKC', 'IGFBP5', 'IFIT1', 'LY6E')

gene_module[["gene_module_4"]]<-c('STMN1', 'H2AFZ', 'UBE2C', 'TUBA1B', 'BIRC5', 'HMGB2', 'ZWINT', 'TUBB', 'HMGB1', 'DEK', 'CDK1', 'HMGN2', 'UBE2T', 'TK1', 'RRM2', 'RANBP1', 'TYMS', 'CENPW', 'MAD2L1', 'CKS2', 'CKS1B', 'NUSAP1', 'TUBA1C', 'PTTG1', 'KPNA2', 'PCNA', 'CENPF', 'HIST1H4C', 'CDKN3', 'UBE2S', 'CCNB1', 'HMGA1', 'DTYMK', 'SNRPB', 'CDC20', 'NASP', 'MCM7', 'PLP2', 'TUBB4B', 'PLK1', 'CCNB2', 'MKI67', 'TOP2A', 'TPX2', 'PKMYT1', 'PRC1', 'SMC4', 'CENPU', 'RAN', 'DUT', 'PA2G4', 'BUB3', 'RAD21', 'SPC25', 'HN1', 'CDCA3', 'H2AFV', 'HNRNPA2B1', 'CCNA2', 'PBK', 'LSM5', 'DNAJC9', 'RPA3', 'TMPO', 'SNRPD1', 'CENPA', 'KIF20B', 'USP1', 'H2AFX', 'PPM1G', 'NUF2', 'SNRPG', 'KIF22', 'KIAA0101', 'DEPDC1', 'RNASEH2A', 'MT2A', 'STRA13', 'ANLN', 'CACYBP', 'NCL', 'NUDT1', 'ECT2', 'LSM4', 'ASF1B', 'CENPN', 'TMEM106C', 'CCT5', 'HSPA8', 'HMMR', 'SRSF3', 'AURKB', 'GGH', 'AURKA', 'TRIP13', 'CDCA8', 'HMGB3', 'HNRNPAB', 'FAM83D', 'CDC25B', 'GGCT', 'KNSTRN', 'CCT6A', 'PTGES3', 'ANP32E', 'CENPK', 'MCM3', 'DDX21', 'HSPD1', 'SKA2', 'CALM2', 'UHRF1', 'HINT1', 'ORC6', 'MZT1', 'MIS18BP1', 'WDR34', 'NAP1L1', 'TEX30', 'SFN', 'HSPE1', 'CENPM', 'TROAP', 'CDCA5', 'RACGAP1', 'SLC25A5', 'ATAD2', 'DBF4', 'KIF23', 'CEP55', 'SIVA1', 'SAC3D1', 'PSIP1', 'CLSPN', 'CCT2', 'DLGAP5', 'PSMA4', 'SMC2', 'AP2S1', 'RAD51AP1', 'MND1', 'ILF2', 'DNMT1', 'NUCKS1', 'LMNB1', 'RFC4', 'EIF5A', 'NPM3', 'ARL6IP1', 'ASPM', 'GTSE1', 'TOMM40', 'HNRNPA1', 'GMNN', 'FEN1', 'CDCA7', 'SLBP', 'TNFRSF12A', 'TM4SF1', 'CKAP2', 'CENPE', 'SRP9', 'DDX39A', 'COMMD4', 'RBM8A', 'CALM3', 'RRM1', 'ENO1', 'ANP32B', 'SRSF7', 'FAM96A', 'TPRKB', 'FABP5', 'PPIF', 'SERPINE1', 'TACC3', 'RBBP7', 'NEK2', 'CALM1', 'GMPS', 'EMP2', 'HMG20B', 'SMC3', 'HSPA9', 'NAA20', 'NUDC', 'RPL39L', 'PRKDC', 'CDCA4', 'HIST1H1A', 'HES6', 'SUPT16H', 'PTMS', 'VDAC3', 'PSMC3', 'ATP5G1', 'PSMA3', 'PGP', 'KIF2C', 'CARHSP1')

gene_module[["gene_module_5"]]<-c('GJA1', 'SCGB2A2', 'ARMT1', 'MAGED2', 'PIP', 'SCGB1D2', 'CLTC', 'MYBPC1', 'PDZK1', 'MGP', 'SLC39A6', 'CCND1', 'SLC9A3R1', 'NAT1', 'SUB1', 'CYP4X1', 'STC2', 'CROT', 'CTSD', 'FASN', 'PBX1', 'SLC4A7', 'FOXA1', 'MCCC2', 'IDH1', 'H2AFJ', 'CYP4Z1', 'IFI27', 'TBC1D9', 'ANPEP', 'DHRS2', 'TFF3', 'LGALS3BP', 'GATA3', 'LTF', 'IFITM2', 'IFITM1', 'AHNAK', 'SEPP1', 'ACADSB', 'PDCD4', 'MUCL1', 'CERS6', 'LRRC26', 'ASS1', 'SEMA3C', 'APLP2', 'AMFR', 'CDV3', 'VTCN1', 'PREX1', 'TP53INP1', 'LRIG1', 'ANK3', 'ACLY', 'CLSTN1', 'GNB1', 'C1orf64', 'STARD10', 'CA12', 'SCGB2A1', 'MGST1', 'PSAP', 'GNAS', 'MRPS30', 'MSMB', 'DDIT4', 'TTC36', 'S100A1', 'FAM208B', 'STT3B', 'SLC38A1', 'DMKN', 'SEC14L2', 'FMO5', 'DCAF10', 'WFDC2', 'GFRA1', 'LDLRAD4', 'TXNIP', 'SCGB3A1', 'APOD', 'N4BP2L2', 'TNC', 'ADIRF', 'NPY1R', 'NBPF1', 'TMEM176A', 'GLUL', 'BMP2K', 'SLC44A1', 'GFPT1', 'PSD3', 'CCNG2', 'CGNL1', 'TMED7', 'NOVA1', 'ARCN1', 'NEK10', 'GPC6', 'SCGB1B2P', 'IGHG4', 'SYT1', 'SYNGR2', 'HSPA1A', 'ATP6AP1', 'TSPAN13', 'MT-ND2', 'NIFK', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'EVL', 'GRN', 'ERH', 'CD81', 'NUPR1', 'SELENBP1', 'C1orf56', 'LMO3', 'PLK2', 'HACD3', 'RBBP8', 'CANX', 'ENAH', 'SCD', 'CREB3L2', 'SYNCRIP', 'TBL1XR1', 'DDR1', 'ERBB3', 'CHPT1', 'BANF1', 'UGDH', 'SCUBE2', 'UQCR10', 'COX6C', 'ATP5G1', 'PRSS23', 'MYEOV2', 'PITX1', 'MT-ND4L', 'TPM1', 'HMGCS2', 'ADIPOR2', 'UGCG', 'FAM129B', 'TNIP1', 'IFI6', 'CA2', 'ESR1', 'TMBIM4', 'NFIX', 'PDCD6IP', 'CRIM1', 'ARHGEF12', 'ENTPD5', 'PATZ1', 'ZBTB41', 'UCP1', 'ANO1', 'RP11-356O9.1', 'MYB', 'ZBTB44', 'SCPEP1', 'HIPK2', 'CDK2AP1', 'CYHR1', 'SPINK8', 'FKBP10', 'ISOC1', 'CD59', 'RAMP1', 'AFF3', 'MT-CYB', 'PPP1CB', 'PKM', 'ALDH2', 'PRSS8', 'NPW', 'SPR', 'PRDX3', 'SCOC', 'TMED10', 'KIAA0196', 'NDP', 'ZSWIM7', 'AP2A1', 'PLAT', 'SUSD3', 'CRABP2', 'DNAJC12', 'DHCR24', 'PPT1', 'FAM234B', 'DDX17', 'LRP2', 'ABCD3', 'CDH1', 'NFIA') 

gene_module[["gene_module_6"]]<-c('AGR2', 'TFF3', 'SELM', 'CD63', 'CTSD', 'MDK', 'CD74', 'S100A13', 'IFITM3', 'HLA-B', 'AZGP1', 'FXYD3', 'IFITM2', 'RABAC1', 'S100A14', 'CRABP2', 'LTF', 'RARRES1', 'HLA-A', 'PPIB', 'HLA-C', 'S100A10', 'S100A9', 'TIMP1', 'DDIT4', 'S100A16', 'LGALS1', 'LAPTM4A', 'SSR4', 'S100A6', 'CD59', 'BST2', 'PDIA3', 'KRT19', 'CD9', 'FXYD5', 'SCGB2A2', 'NUCB2', 'TMED3', 'LY6E', 'CFD', 'ITM2B', 'PDZK1IP1', 'LGALS3', 'NUPR1', 'SLPI', 'CLU', 'TMED9', 'HLA-DRA', 'SPTSSB', 'TMEM59', 'KRT8', 'CALR', 'HLA-DRB1', 'IFI6', 'NNMT', 'CALML5', 'S100P', 'TFF1', 'ATP1B1', 'SPINT2', 'PDIA6', 'S100A8', 'HSP90B1', 'LMAN1', 'RARRES3', 'SELENBP1', 'CEACAM6', 'TMEM176A', 'EPCAM', 'MAGED2', 'SNCG', 'DUSP4', 'CD24', 'PERP', 'WFDC2', 'HM13', 'TMBIM6', 'C12orf57', 'DKK1', 'MAGED1', 'PYCARD', 'RAMP1', 'C11orf31', 'STOM', 'TNFSF10', 'BSG', 'TMED10', 'ASS1', 'PDLIM1', 'CST3', 'PDIA4', 'NDUFA4', 'GSTP1', 'TYMP', 'SH3BGRL3', 'PRSS23', 'P4HA1', 'MUC5B', 'S100A1', 'PSAP', 'TAGLN2', 'MGST3', 'PRDX5', 'SMIM22', 'NPC2', 'MESP1', 'MYDGF', 'ASAH1', 'APP', 'NGFRAP1', 'TMEM176B', 'C8orf4', 'KRT81', 'VIMP', 'CXCL17', 'MUC1', 'COMMD6', 'TSPAN13', 'TFPI', 'C15orf48', 'CD151', 'TACSTD2', 'PSME2', 'CLDN7', 'ATP6AP2', 'CUTA', 'MT2A', 'CYB5A', 'CD164', 'TM4SF1', 'SCGB1D2', 'GSTM3', 'EGLN3', 'LMAN2', 'IFI27', 'PPP1R1B', 'B2M', 'ANXA2', 'SARAF', 'MUCL1', 'CSRP1', 'NPW', 'SLC3A2', 'PYDC1', 'QSOX1', 'TSPAN1', 'GPX1', 'TMSB4X', 'FGG', 'GUK1', 'IL32', 'ATP6V0E1', 'BCAP31', 'CHCHD10', 'TSPO', 'TNFRSF12A', 'MT1X', 'PDE4B', 'HSPA5', 'SCD', 'SERINC2', 'PSCA', 'VAMP8', 'ELF3', 'TSC22D3', 'S100A7', 'GLUL', 'ZG16B', 'TMEM45A', 'APMAP', 'RPS26', 'CALU', 'OSTC', 'NCCRP1', 'SQLE', 'RPS28', 'SSR2', 'SOX4', 'CLEC3A', 'TMEM9', 'RPL10', 'MUC5AC', 'HLA-DPA1', 'ZNHIT1', 'AQP5', 'CAPG', 'SPINT1', 'NDFIP1', 'FKBP2', 'C1S', 'LDHA', 'NEAT1', 'RPL36A', 'S100A11', 'LCN2', 'TUBA1A', 'GSTK1', 'SEPW1', 'P4HB') 

gene_module[["gene_module_7"]]<-c('KCNQ1OT1', 'AKAP9', 'RHOB', 'SOX4', 'VEGFA', 'CCNL1', 'RSRP1', 'RRBP1', 'ELF3', 'H1FX', 'FUS', 'NEAT1', 'N4BP2L2', 'SLC38A2', 'BRD2', 'PNISR', 'CLDN4', 'MALAT1', 'SOX9', 'DDIT3', 'TAF1D', 'FOSB', 'ZNF83', 'ARGLU1', 'DSC2', 'MACF1', 'GTF2I', 'SEPP1', 'ANKRD30A', 'PRLR', 'MAFB', 'NFIA', 'ZFAS1', 'MTRNR2L12', 'RNMT', 'NUPR1', 'MT-ND6', 'RBM39', 'HSPA1A', 'HSPA1B', 'RGS16', 'SUCO', 'XIST', 'PDIA6', 'VMP1', 'SUGP2', 'LPIN1', 'NDRG1', 'PRRC2C', 'CELF1', 'HSP90B1', 'JUND', 'ACADVL', 'PTPRF', 'LMAN1', 'HEBP2', 'ATF3', 'BTG1', 'GNAS', 'TSPYL2', 'ZFP36L2', 'RHOBTB3', 'TFAP2A', 'RAB6A', 'KMT2C', 'POLR2J3', 'CTNND1', 'PRRC2B', 'RNF43', 'CAV1', 'RSPO3', 'IMPA2', 'FAM84A', 'FOS', 'IGFBP5', 'NCOA3', 'WSB1', 'MBNL2', 'MMP24-AS1', 'DDX5', 'AP000769.1', 'MIA3', 'ID2', 'HNRNPH1', 'FKBP2', 'SEL1L', 'PSAT1', 'ASNS', 'SLC3A2', 'EIF4EBP1', 'HSPH1', 'SNHG19', 'RNF19A', 'GRHL1', 'WBP1', 'SRRM2', 'RUNX1', 'ASH1L', 'HIST1H4C', 'RBM25', 'ZNF292', 'RNF213', 'PRPF38B', 'DSP', 'EPC1', 'FNBP4', 'ETV6', 'SPAG9', 'SIAH2', 'RBM33', 'CAND1', 'CEBPB', 'CD44', 'NOC2L', 'LY6E', 'ANGPTL4', 'GABPB1-AS1', 'MTSS1', 'DDX42', 'PIK3C2G', 'IAH1', 'ATL2', 'ADAM17', 'PHIP', 'MPZ', 'CYP27A1', 'IER2', 'ACTR3B', 'PDCD4', 'COLCA1', 'KIAA1324', 'TFAP2C', 'CTSC', 'MYC', 'MT1X', 'VIMP', 'SERHL2', 'YPEL3', 'MKNK2', 'ZNF552', 'CDH1', 'LUC7L3', 'DDIT4', 'HNRNPR', 'IFRD1', 'RASSF7', 'SNHG8', 'EPB41L4A-AS1', 'ZC3H11A', 'SNHG15', 'CREB3L2', 'ERBB3', 'THUMPD3-AS1', 'RBBP6', 'GPBP1', 'NARF', 'SNRNP70', 'RP11-290D2.6', 'SAT1', 'GRB7', 'H1F0', 'EDEM3', 'KIAA0907', 'ATF4', 'DNAJC3', 'DKK1', 'SF1', 'NAMPT', 'SETD5', 'DYNC1H1', 'GOLGB1', 'C4orf48', 'CLIC3', 'TECR', 'HOOK3', 'WDR60', 'TMEM101', 'SYCP2', 'C6orf62', 'METTL12', 'HIST1H2BG', 'PCMTD1', 'PWWP2A', 'HIST1H3H', 'NCK1', 'CRACR2B', 'NPW', 'RAB3GAP1', 'TMEM63A', 'MGP', 'ANKRD17', 'CALD1', 'PRKAR1A', 'PBX1', 'ATXN2L', 'FAM120A', 'SAT2', 'TAF10', 'SFRP1', 'CITED2') 

#maybe add Ecotypes as well, no given gene lists from the publication??



sample_cell_signature_transfer<-function(){
  dat_epi<-subset(dat,EMBO_predicted.id=="epithelial")

  #embo lineage
  dat_epi<-AddModuleScore(dat_epi,features=lineage_in,names=names(lineage_in),assay="SoupXRNA",seed=123,search=TRUE)

  colnames(dat_epi@meta.data)[which(colnames(dat_epi@meta.data) %in% c("Cluster1","Cluster2","Cluster3","Cluster4"))]<-c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str") #Rename them

  #Immune cell features and PAM50 canonical short list of genes
  for(i in 1:length(features_in)){
    features_in[[i]]<-features_in[[i]][features_in[[i]] %in% row.names(dat_epi[["SoupXRNA"]])] #make sure gene names match
    dat_epi<-MetaFeature(dat_epi,features=c(features_in[[i]]),meta.name=names(features_in)[i],assay="SoupXRNA")}

  #SCSubype List of genes
  #run only on epithelial cells
  module_scores<-AddModuleScore(dat_epi,features=module_feats,assay="SoupXRNA",search=TRUE,name=names(module_feats)) #use add module function to add cell scores
  module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-(length(module_feats)-1),ncol(module_scores@meta.data))]
  colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like
  dat_epi<-AddMetaData(dat,metadata=module_scores)

  #Swarbrick Gene Modules
  #run only on epithelial cells
  gene_module_out<-AddModuleScore(dat_epi,features=gene_module,assay="SoupXRNA",search=TRUE,name=names(gene_module)) #use add module function to add cell scores
  gene_module_out<-gene_module_out@meta.data[seq(ncol(gene_module_out@meta.data)-(length(gene_module)-1),ncol(gene_module_out@meta.data))]#get the 7 added gene modules
  colnames(gene_module_out)<-names(gene_module) 
  dat_epi<-AddMetaData(dat_epi,metadata=gene_module_out)
  out<-dat_epi@meta.data[c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str",names(module_feats),names(gene_module))]
  return(out)
}

single_sample_PAM50_assignment<-function(){
  met<-dat@meta.data
  met<-met[met$EMBO_predicted.id %in% c("epithelial"),]
  pam50_list<-  c("EMBO_Basal","EMBO_LP","EMBO_ML","EMBO_Str" )
  max_pam50<-lapply(1:nrow(met),function(i) pam50_list[which(met[i,pam50_list]==max(met[i,pam50_list],na.rm=T))])
  max_pam50<-unlist(lapply(1:length(max_pam50),function(i) do.call("paste",as.list(max_pam50[[i]]))))
  max_pam50<-unlist(lapply(max_pam50,function(i) gsub("EMBO_","",i)))
  names(max_pam50)<-row.names(met)
  return(max_pam50)
}

single_sample_SCtype_assignment<-function(){
  met<-dat@meta.data
  met<-met[met$EMBO_predicted.id %in% c("epithelial"),]
  scsubtype_list<-  c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")
  max_scsubtype<-lapply(1:nrow(met),function(i) scsubtype_list[which(met[i,scsubtype_list]==max(met[i,scsubtype_list],na.rm=T))])
  max_scsubtype<-unlist(lapply(1:length(max_scsubtype),function(i) do.call("paste",as.list(max_scsubtype[[i]]))))
  names(max_scsubtype)<-row.names(met)
  return(max_scsubtype)
}

#Generate scores per epithelial cell
epithelial_metadata<-sample_cell_signature_transfer()
dat<-AddMetaData(dat,metadata=epithelial_metadata) #add to master data frame metadata

#assign top PAM50 designation by epithelial cells
max_pam50<-single_sample_SCtype_assignment()
dat<-AddMetaData(dat,max_pam50,col.name="PAM50_epi_designation")

#assign top scsubtype by epithelial cells
max_scsubtype<-single_sample_PAM50_assignment()
dat<-AddMetaData(dat,max_scsubtype,col.name="SCSubtype_epi_designation")

saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")

```