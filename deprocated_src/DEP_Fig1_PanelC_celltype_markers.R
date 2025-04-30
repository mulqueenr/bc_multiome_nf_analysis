module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggalluvial)
library(reshape2)
library(optparse)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
set.seed(123)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.passqc.SeuratObject.rds"

# #snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT")
hbca_snmarkers[["basal"]]=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
#hbca_snmarkers[["macro_mono"]]=c("CXCL5","SERPINB2","EREG","VCAN","C3","FCGBP","SDS","APOE","RNASE1","HMOX1","LYVE1","SELENOP")
hbca_snmarkers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
hbca_snmarkers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
hbca_snmarkers[["pdc"]]=c("IGKC","PTGDS","IRF8","DNASE1L3","LGALS2","C1orf54","CLIC3")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
#hbca_snmarkers[["treg"]]=c("TIGIT","TBC1D4","BATF","CARD16")
#hbca_snmarkers[["cd4_8_tcell"]]=c("HSPA1A","HSPA1B","DNAJB1","HSPA6","CCL5","ZNF683","TRAF3IP3")
#hbca_snmarkers[["nkt"]]=c("GNLY","NKG7","FGFBP2","CCL3","GZMH","CCL4")
features<-llply(hbca_snmarkers, unlist)

Idents(dat)<-factor(dat$assigned_celltype,levels=c("luminal_epithelial","basal_epithelial",
"adipocyte","endothelial_vascular","endothelial_lymphatic","pericyte","fibroblast",
"mast","myeloid","bcell","plasma","pDC",
"tcell"
))
DefaultAssay(dat)<-"RNA"
p10<-DotPlot(subset(dat,cells=names(Idents(dat))[!is.na(Idents(dat))]),features=features,cluster.idents=FALSE,dot.scale=8)+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))+
  theme(axis.text.x = element_text(angle=90))


ggsave(p10,file="assigned_celltypes.features.pdf",height=10,width=40,limitsize=F)

