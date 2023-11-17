library(Seurat)
library(ggplot2)
library(Signac)
library(plyr)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)


setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
dat_merged<-readRDS("phase2.QC.filt.SeuratObject.rds")

#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT")
hbca_snmarkers[["basal"]]=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2")
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
features<-llply(hbca_snmarkers, unlist)
dat<-dat_merged

#Idents(dat)<-factor(dat$EMBO_predicted.id,levels=rev(c("epithelial","cycling.epithelial","CAFs","Endothelial","Pericytes","Myeloid","TAMs","TAMs_2","Plasma.cells","B.cells","T.cells","NA")))
Idents(dat)<-factor(dat$HBCA_predicted.id,levels=c("luminal epithelial cell of mammary gland",
  "basal cell","fibroblast","endothelial cell of lymphatic vessel",
  "endothelial cell of vascular tree","pericyte","myeloid cell","T cell","mast cell","adipocyte of breast" ))

plt1<-DotPlot(dat,assay = "RNA",features, dot.scale = 10,cluster.idents = FALSE)+ RotatedAxis()+ scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))

ggsave(plt1,file="hbca_markers.pdf",height=10,width=30,limitsize = FALSE)
system("slack -F hbca_markers.pdf ryan_todo")
#I dont see a AC044810.3 gene, so I'm assuming its a typo and its AC044810.2??

#plt2<-DotPlot(dat,assay = "GeneActivity",features,cols=c("grey","blue"), dot.scale = 5,cluster.idents = FALSE,)+ RotatedAxis()

#scRNA TF markers https://www.nature.com/articles/s41586-023-06252-9/figures/7
hbca_sctfs=list()
hbca_sctfs[["lumhr"]]=c("BATF","AR","FOXA1","MYB","PGR","TBX3","SPDEF")
hbca_sctfs[["lumsec"]]=c("MESP1","FOXI1","KLF5")
hbca_sctfs[["basal"]]=c("KLF16","TP63","IRX2","SIX4","IRX4","LMX1B","MAFF")
hbca_sctfs[["fibro"]]=c("TWIST2","NFATC2","WISP2","CREB3L1","BHLHE22","ALX4","TWIST1")
hbca_sctfs[["lymphatic"]]=c("TBX1","NFATC1","HOXD10","TAL1","TFF3","KLF2","HOXD8")
hbca_sctfs[["vascular"]]=c("SOX17","NPDC1","MECOM","HOXB9","SOX7","HOXD1","HOXB6")
hbca_sctfs[["perivasc"]]=c("SOX5","SOX13","EMX2","NFIL3","GLIS2","FOXL1","BCL6")
hbca_sctfs[["myeloid"]]=c("SP1","TFEC","NFKB1","RXRA","IRF5","MXD1","MAFB")
hbca_sctfs[["tcells"]]=c("SNAI3","TBX21","EOMES","ZNF831","LEF1","PBX4","ETS1")
hbca_sctfs[["bcells"]]=c("BCL11A","CREB3","POU2F2","CREB3L2","PAX5","IRF4","SPIB")
features<-llply(hbca_sctfs, unlist)
features_2<-ConvertMotifID(object=dat,name=unlist(features))
Idents(dat)<-factor(dat$EMBO_predicted.id,levels=rev(c("epithelial","cycling.epithelial","CAFs","Endothelial","Pericytes","Myeloid","TAMs","TAMs_2","Plasma.cells","B.cells","T.cells","NA")))
plt1<-DotPlot(dat,assay = "chromvar",features_2,cols=c("grey","red"), dot.scale = 5,cluster.idents = FALSE,)+ RotatedAxis()
plt2<-DotPlot(dat,assay = "GeneActivity",features,cols=c("grey","blue"), dot.scale = 5,cluster.idents = FALSE,)+ RotatedAxis()
ggsave(plt1/plt2,file="hbca_markers_sctfs.pdf",height=20,width=30,limitsize = FALSE)
system("slack -F hbca_markers_sctfs.pdf ryan_todo")

hbca_sntfs=list()
hbca_sntfs[["lumhr"]]=c("FOXA1","XBP1","MYB","SPDEF","GATA3","RUNX1","ESR1","KLF5")
hbca_sntfs[["lumsec"]]=c("BCL11A","SMARCA4","SOX9","GRHL1","EHF","RARG")
hbca_sntfs[["basal"]]=c("VDR","BACH2","HLF","ZNF561","TEAD3","RARB","MXI1")
hbca_sntfs[["fibro"]]=c("GLI2","HMGA2","WISP2","PLAGL1","ZEB1","FOXP2","SMAD3")
hbca_sntfs[["lymphatic"]]=c("NFATC1","ZMAT4","HVEP2","ELK3","NFATC2","SMAD1","CREB5")
hbca_sntfs[["vascular"]]=c("NR5A2","PML","ELV4","HVEP1","ERG","SMAD9","ETS2")
hbca_sntfs[["perivasc"]]=c("NR2F2","NFATC4","EBF2","PRDM16","PRRX1","NR1H3","RXRA")
hbca_sntfs[["myeloid"]]=c("TFEC","IRF8","MAFB","SPI1","CEBPB","MAF","JDP2")
hbca_sntfs[["tcells"]]=c("SP4","CTCF","TCF7","RUNX3","PRDM1","MYBL1","RUNX2")
hbca_sntfs[["mast"]]=c("GATA2","NFKB1","BATF","NFKB2","FOSB","GABPB1")
hbca_sntfs[["adipo"]]=c("STAT5A","TEAD1","BCL6","FOXO1","MLXIPL","SOX6","STAT5B")
features<-llply(hbca_sntfs, unlist)
features_2<-llply(features,function(x)ConvertMotifID(object=dat,name=x))
features_2
Idents(dat)<-factor(dat$HBCA_predicted.id,levels=c("luminal epithelial cell of mammary gland",
  "basal cell","fibroblast","endothelial cell of lymphatic vessel",
  "endothelial cell of vascular tree","pericyte","myeloid cell","T cell","mast cell","adipocyte of breast" ))
plt1<-DotPlot(dat,assay = "chromvar",features_2, dot.scale = 10,cluster.idents = TRUE)+ RotatedAxis() + scale_color_gradient2(low="#7f3b08",mid="white",high="#542788")
plt2<-DotPlot(dat,assay = "GeneActivity",features, dot.scale = 10,cluster.idents = TRUE)+ RotatedAxis()+ scale_color_gradient2(low="#7f3b08",mid="white",high="#542788")
ggsave(plt1/plt2,file="hbca_markers_sntfs.pdf",height=20,width=30,limitsize = FALSE)
system("slack -F hbca_markers_sntfs.pdf ryan_todo")






#Running TFs deteremined through our data set
chromvar_markers<-FindAllMarkers(dat,assay="chromvar", min.pct = 0.1)
chromvar_markers_listed=list()
for(x in unique(chromvar_markers$cluster)){
  tmp<-head(chromvar_markers[chromvar_markers$cluster==x,],n=3)
  chromvar_markers_listed[[x]]<-tmp$gene
  }
features_2<-llply(chromvar_markers_listed,function(x)ConvertMotifID(object=dat,id=x))

plt1<-DotPlot(dat,assay = "chromvar",unique(chromvar_markers_listed), dot.scale = 10,cluster.idents = TRUE)+ RotatedAxis() + scale_color_gradient2(low="#7f3b08",mid="white",high="#542788")+scale_x_discrete(labels=c(features_2))
ggsave(plt1,file="hbca_markers_sntfs.pdf",height=20,width=30,limitsize = FALSE)
system("slack -F hbca_markers_sntfs.pdf ryan_todo")


#Running TFs deteremined through our data set
rna_markers<-FindAllMarkers(dat,assay="SoupXRNA", min.pct = 0.1)
rna_markers_listed=list()
for(x in unique(rna_markers$cluster)){
  tmp<-head(rna_markers[rna_markers$cluster==x,],n=3)
  rna_markers_listed[[x]]<-tmp$gene
  }
plt1<-DotPlot(dat,assay = "SoupXRNA",rna_markers_listed, dot.scale = 10,cluster.idents = TRUE)+ RotatedAxis() + scale_color_gradient2(low="#7f3b08",mid="white",high="#542788")+ scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))
ggsave(plt1,file="hbca_markers_denovo_rna.pdf",height=20,width=30,limitsize = FALSE)
system("slack -F hbca_markers_denovo_rna.pdf ryan_todo")




