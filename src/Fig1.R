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
opt$object_input="merged.geneactivity.SeuratObject.rds"

###########PANEL B UPDATE ################
dat<-readRDS(opt$object_input)
dat[["RNA"]] <- as(dat[["RNA"]], Class = "Assay5")

#clustering function for processing stuff
multimodal_cluster<-function(dat=dat,prefix="allcells",res=0.5,dotsize=5){
  # Perform standard analysis of each modality independently 
  #RNA analysis
  DefaultAssay(dat) <- 'RNA'
  dat<-NormalizeData(dat) %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  dat <- RunUMAP(dat, 
    reduction="pca", 
    dims = 1:30, 
    reduction.name = paste(prefix,"umap","rna",sep="."),
    reduction.key = "rnaUMAP_")

  #ATAC analysis
  DefaultAssay(dat)<- 'peaks'
  dat<-RunTFIDF(dat) %>%  FindTopFeatures() %>% RunSVD()
  dat <- RunUMAP(dat, 
    reduction = "lsi", 
    dims = 2:30, 
    reduction.name=paste(prefix,"umap","atac",sep="."), 
    reduction.key = "atacUMAP_")

  # build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(object = dat,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:40),
    modality.weight.name = "RNA.weight",
    weighted.nn.name=paste(prefix,"weighted.nn",sep="."),
    snn.graph.name=paste(prefix,"wsnn",sep="."),
    verbose = TRUE)

  dat <- RunUMAP(dat, 
    nn.name = paste(prefix,"weighted.nn",sep="."), 
    reduction.name = paste(prefix,"wnn.umap",sep="."), 
    reduction.key = "wnnUMAP_")

  dat <- FindClusters(dat, 
    graph.name = paste(prefix,"wsnn",sep="."), 
    algorithm = 3, 
    resolution = res, 
    verbose = FALSE)
  p1<-DimPlot(dat, pt.size=dotsize,group.by='HBCA_predicted.id',label = TRUE, repel = TRUE, reduction = paste(prefix,"umap.rna",sep="."),raster=T) + NoLegend()
  p2<-DimPlot(dat, pt.size=dotsize,group.by='HBCA_predicted.id',label = TRUE, repel = TRUE, reduction = paste(prefix,"umap.atac",sep="."),raster=T) + NoLegend()
  p3<-DimPlot(dat, pt.size=dotsize,group.by='HBCA_predicted.id',label = TRUE, repel = TRUE, reduction = paste(prefix,"wnn.umap",sep="."),raster=T) + NoLegend()

  p4<-DimPlot(dat, pt.size=dotsize,group.by = 'seurat_clusters',label = TRUE, repel = TRUE, reduction = paste(prefix,"umap.rna",sep="."),raster=T) + NoLegend()
  p5<-DimPlot(dat, pt.size=dotsize,group.by = 'seurat_clusters',label = TRUE, repel = TRUE, reduction = paste(prefix,"umap.atac",sep="."),raster=T) + NoLegend()
  p6<-DimPlot(dat, pt.size=dotsize,group.by = 'seurat_clusters',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=T) + NoLegend()
  return(c(dat,list(p1,p2,p3,p4,p5,p6)))
  }

#initial run before qc filtering
#set filters for all cells scrublet score <0.35, ncount_RNA>1000, ncount_ATAC>1000
prefix="allcells"
out<-multimodal_cluster(dat,prefix=prefix,res=0.5)
dat<-out[[1]]
plt_out<-patchwork::wrap_plots(out[2:length(out)], nrow = 2, ncol = 3)
p7<-FeaturePlot(dat,features=c("scrublet_Scores","nCount_RNA","nCount_peaks"),reduction = "allcells.wnn.umap",ncol=3)+ NoLegend()
dat$scrublet_DropletType<-ifelse(as.numeric(dat$scrublet_Scores)<0.35,"singlet","doublet")
p8<-DimPlot(dat,group.by="scrublet_DropletType",reduction = "allcells.wnn.umap")+ NoLegend()
dat$passqc<-ifelse(as.numeric(dat$scrublet_Scores)<0.35 & as.numeric(dat$nCount_RNA)>1000 & as.numeric(dat$nCount_peaks)>1000,"PASS","FAIL")
p9<-DimPlot(dat,group.by="passqc",reduction = "allcells.wnn.umap")+ NoLegend()
plt<-plt_out/(p7 )/(p8+p9)
ggsave(plt,file="allcells.umap.pdf",height=40,width=30)
saveRDS(dat,file="merged.geneactivity.SeuratObject.rds")


scrublet_filter_check<-function(scrublet_score){
  dat$passqc<-ifelse(as.numeric(dat$scrublet_Scores)<scrublet_score & as.numeric(dat$nCount_RNA)>1000 & as.numeric(dat$nCount_peaks)>1000,"PASS","FAIL")
  dat2<-subset(dat,passqc=="PASS")
  cell_count=nrow(dat2@meta.data)
  out<-multimodal_cluster(dat2,prefix=prefix,res=0.5)
  dat2<-out[[1]]
  plt<-DimPlot(dat2,group.by="HBCA_predicted.id",reduction = "allcells.wnn.umap")+ NoLegend()+ggtitle(paste("Scrublet:",scrublet_score,"Cells:",cell_count))
  return(plt)

}

plt_list<-lapply(c(0.35,0.2,0.1),scrublet_filter_check)
plt_out<-patchwork::wrap_plots(plt_list, nrow = 3, ncol = 1)
ggsave(plt_out,file="allcells.umap.passqc.pdf",height=30,width=10)

#subset to cells passing qc.
prefix="allcells_passqc"
dat<-subset(dat,passqc=="PASS")
out<-multimodal_cluster(dat,prefix="allcells_passqc")
dat<-out[[1]]
plt_out<-patchwork::wrap_plots(out[2:length(out)], nrow = 2, ncol = 3)
dat$log10_nCount_RNA<-log10(dat$nCount_RNA)
dat$log10_nCount_peaks<-log10(dat$nCount_peaks)

p7<-FeaturePlot(dat,features=c("scrublet_Scores","log10_nCount_RNA","log10_nCount_peaks"),reduction = paste(prefix,"wnn.umap",sep="."),ncol=3,raster=T)+ NoLegend()
dat$high_cell_count_sample<-ifelse(dat@meta.data$sample %in% names(which(table(dat$sample)>10000)),"TRUE","FALSE")
dat$sample_cell_count<-table(dat$sample)[dat@meta.data$sample]

p8<-DimPlot(dat,group.by="high_cell_count_sample",reduction = "allcells.wnn.umap")+ggtitle("High Cell Count Sample")
p9<-FeaturePlot(dat,features=c("sample_cell_count"),reduction = "allcells.wnn.umap")+ggtitle("Sample Cell Count")

plt<-(p7/p8 /p9)
ggsave(plt,file="allcells.umap.passqc_cellcountsamples.pdf",height=40,width=30)

plt<-ggplot(dat@meta.data,aes(x=high_cell_count_sample,y=scrublet_Scores))+geom_jitter(aes(color=as.numeric(dat$sample_cell_count)))+geom_violin()
ggsave(plt,file="cellcountsamples_scrublet_scores.pdf",height=40,width=30)

#Based on this using scrublet cutoff of 0.2
#initial run before qc filtering
#set filters for all cells scrublet score <0.35, ncount_RNA>1000, ncount_ATAC>1000
prefix="allcells"
dat$passqc<-ifelse(as.numeric(dat$scrublet_Scores)<0.2 & as.numeric(dat$nCount_RNA)>1000 & as.numeric(dat$nCount_peaks)>1000,"PASS","FAIL")
dat<-subset(dat,cell=row.names(dat@meta.data[dat$passqc=="PASS",]))

out<-multimodal_cluster(dat,prefix=prefix,res=0.5)
dat<-out[[1]]
plt_out<-patchwork::wrap_plots(out[2:length(out)], nrow = 2, ncol = 3)
p7<-FeaturePlot(dat,features=c("scrublet_Scores","log10_nCount_RNA","log10_nCount_peaks"),reduction = "allcells.wnn.umap",ncol=3)+ NoLegend()
dat$scrublet_DropletType<-ifelse(as.numeric(dat$scrublet_Scores)<0.2,"singlet","doublet")
p8<-DimPlot(dat,group.by="scrublet_DropletType",reduction = "allcells.wnn.umap")+ NoLegend()
p9<-DimPlot(dat,group.by="passqc",reduction = "allcells.wnn.umap")+ NoLegend()
plt<-plt_out/(p7 )/(p8+p9)
ggsave(plt,file="allcells.umap.pdf",height=40,width=30)
saveRDS(dat,file="merged.geneactivity.passqc.SeuratObject.rds")

#Make stacked barplot on identities per cluster
DF<-as.data.frame(dat@meta.data %>% group_by(seurat_clusters,HBCA_predicted.id) %>% tally())
clus_by_celltype<-reshape2::dcast(DF,seurat_clusters~HBCA_predicted.id,value.var="n",fill=0)
row.names(clus_by_celltype)<-paste0("cluster_",clus_by_celltype[,1])
clus_by_celltype<-clus_by_celltype[,2:ncol(clus_by_celltype)]
pdf("allcells.cluster_by_celltype.pdf")
Heatmap(scale(t(as.matrix(clus_by_celltype))))
dev.off()

DF<-as.data.frame(dat@meta.data %>% group_by(seurat_clusters,sample) %>% tally())
clus_by_sample<-reshape2::dcast(DF,seurat_clusters~sample,value.var="n",fill=0)
row.names(clus_by_sample)<-paste0("cluster_",clus_by_sample[,1])
clus_by_sample<-clus_by_sample[,2:ncol(clus_by_sample)]
pdf("allcells.cluster_by_sample.pdf")
Heatmap(scale(t(as.matrix(clus_by_sample))))
dev.off()

DF<-as.data.frame(dat@meta.data %>% group_by(seurat_clusters,sample,HBCA_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=seurat_clusters,fill=HBCA_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()
plt2<-ggplot(DF,aes(x=seurat_clusters,fill=sample,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()
plt3<-ggplot(dat@meta.data,aes(x=seurat_clusters,y=scrublet_Scores))+geom_violin()+theme_minimal()
ggsave(plt1/plt2/plt3,file="allcells.cluster_barplots.pdf",width=50,limitsize=F)

scrublet_scores<-dat@meta.data %>% group_by(seurat_clusters) %>% summarize(mean=mean(scrublet_Scores))
#42 43 35 5 doublet clusters (>0.15 score mean for cluster)
prefix="allcells"
dat$passqc<-ifelse(dat$seurat_clusters %in% c("42","43","35","5"),"FAIL","PASS")
dat<-subset(dat,cell=row.names(dat@meta.data[dat$passqc=="PASS",]))

out<-multimodal_cluster(dat,prefix=prefix,res=0.5)
dat<-out[[1]]
plt_out<-patchwork::wrap_plots(out[2:length(out)], nrow = 2, ncol = 3)
p7<-FeaturePlot(dat,features=c("scrublet_Scores","log10_nCount_RNA","log10_nCount_peaks"),reduction = "allcells.wnn.umap",ncol=3)+ NoLegend()
dat$scrublet_DropletType<-ifelse(as.numeric(dat$scrublet_Scores)<0.2,"singlet","doublet")
p8<-DimPlot(dat,group.by="scrublet_DropletType",reduction = "allcells.wnn.umap")+ NoLegend()
p9<-DimPlot(dat,group.by="passqc",reduction = "allcells.wnn.umap")+ NoLegend()
plt<-plt_out/(p7 )/(p8+p9)
ggsave(plt,file="allcells.umap.pdf",height=40,width=30)
saveRDS(dat,file="merged.geneactivity.passqc.SeuratObject.rds")


#Make stacked barplot on identities per cluster


DF<-as.data.frame(dat@meta.data %>% group_by(seurat_clusters,sample,HBCA_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=seurat_clusters,fill=HBCA_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()
plt2<-ggplot(DF,aes(x=seurat_clusters,fill=sample,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()
plt3<-ggplot(dat@meta.data,aes(x=seurat_clusters,y=scrublet_Scores))+geom_violin()+theme_minimal()
ggsave(plt1/plt2/plt3,file="allcells.cluster_barplots.pdf",width=50,limitsize=F)

#assign clusters by predicted ID to start
dat$reclust<-"luminal_epithelial"
dat@meta.data[dat$seurat_clusters=="35",]$reclust<-"pericyte"
dat@meta.data[dat$seurat_clusters %in% c("13","16","19","41"),]$reclust<-"fibroblast"
dat@meta.data[dat$seurat_clusters %in% c("5","18","31"),]$reclust<-"myeloid"
dat@meta.data[dat$seurat_clusters %in% c("17","21","34"),]$reclust<-"basal_epithelial"
dat@meta.data[dat$seurat_clusters %in% c("8"),]$reclust<-"tcell"
dat@meta.data[dat$seurat_clusters %in% c("12"),]$reclust<-"endothelial_vascular"
dat@meta.data[dat$seurat_clusters %in% c("40"),]$reclust<-"adipocyte"
dat@meta.data[dat$seurat_clusters %in% c("36"),]$reclust<-"endothelial_lymphatic"

p1<-DimPlot(dat,group.by="seurat_clusters",reduction = "allcells.wnn.umap")
p2<-DimPlot(dat,group.by="HBCA_predicted.id",reduction = "allcells.wnn.umap")
p3<-DimPlot(dat,group.by="reclust",reduction = "allcells.wnn.umap")
ggsave(p1/p2/p3,file="umap_reclust.pdf",width=10,height=30)

saveRDS(dat,file="merged.geneactivity.passqc.SeuratObject.rds")

dat_nonepi<-subset(dat,cell=row.names(dat@meta.data)[!(dat$reclust %in% c("luminal_epithelial","basal_epithelial"))])
prefix="nonepi_passqc"
out<-multimodal_cluster(dat_nonepi,prefix="nonepi_passqc",res=0.8)
dat_nonepi<-out[[1]]
p9<-DimPlot(dat_nonepi, group.by = 'reclust',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) + NoLegend()
p10<-DimPlot(dat_nonepi, group.by = 'sample',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) + NoLegend()
p11<-DimPlot(dat_nonepi, group.by = 'seurat_clusters',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) + NoLegend()
plt<-(p9+p10+p11)
ggsave(plt,file="nonepi.umap.passqc.pdf",height=10,width=30)
Idents(dat_nonepi)<-dat_nonepi$seurat_clusters
dat_nonepi<-JoinLayers(dat_nonepi,assay="RNA")
marker_sets<-FindAllMarkers(dat_nonepi,assay="RNA",only.pos=T,logfc.threashold=1,method="roc")
write.table(marker_sets,file="nonepi_subcluster_passqc.markers.tsv",sep="\t",col.names=T)
saveRDS(dat_nonepi,file="merged.geneactivity.passqc.nonepi.SeuratObject.rds")

dat<-readRDS(file="merged.geneactivity.passqc.SeuratObject.rds")
dat_nonepi<-readRDS(file="merged.geneactivity.passqc.nonepi.SeuratObject.rds")
#assign clusters by predicted ID to start
dat_nonepi$assigned_celltype<-"myeloid"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("10","17"),]$assigned_celltype<-"plasma"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("13"),]$assigned_celltype<-"low_quality"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("18"),]$assigned_celltype<-"bcell"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("27"),]$assigned_celltype<-"pDC"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("28"),]$assigned_celltype<-"mast"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("16","1"),]$assigned_celltype<-"tcell"#"cd4_8_tcell"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("19","9"),]$assigned_celltype<-"tcell"#"treg"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("9"),]$assigned_celltype<-"tcell"#"natural_killer_t"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("23"),]$assigned_celltype<-"endothelial_lymphatic"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("5","8","20"),]$assigned_celltype<-"endothelial_vascular"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("25"),]$assigned_celltype<-"adipocyte"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("21"),]$assigned_celltype<-"pericyte"
dat_nonepi@meta.data[dat_nonepi$seurat_clusters %in% c("24","11","15","4","12","24","0","26","2"),]$assigned_celltype<-"fibroblast"

dat_nonepi$Diag_MolDiag<-paste(dat_nonepi$Diagnosis,dat_nonepi$Mol_Diagnosis)
p9<-DimPlot(dat_nonepi, group.by = 'reclust',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) 
p10<-DimPlot(dat_nonepi, group.by = 'sample',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) 
p11<-DimPlot(dat_nonepi, group.by = 'seurat_clusters',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) 
p12<-DimPlot(dat_nonepi, group.by = 'assigned_celltype',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) 
p13<-DimPlot(dat_nonepi, group.by = 'Diagnosis',label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) 
p14<-DimPlot(dat_nonepi, group.by = "Diag_MolDiag",label = TRUE, repel = TRUE,reduction = paste(prefix,"wnn.umap",sep="."),raster=F) 

plt<-(p9+p10+p11+p12+p13+p14)
ggsave(plt,file="nonepi.umap.passqc.pdf",height=30,width=50,limitsize=F)

dat$assigned_celltype<-dat$reclust
dat<-AddMetaData(dat,col.name="assigned_celltype",dat_nonepi$assigned_celltype)

dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)

p9<-DimPlot(dat, group.by = 'reclust',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
p10<-DimPlot(dat, group.by = 'sample',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
p11<-DimPlot(dat, group.by = 'seurat_clusters',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
p12<-DimPlot(dat, group.by = 'assigned_celltype',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
p13<-DimPlot(dat, group.by = 'Diagnosis',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
p14<-DimPlot(dat, group.by = 'Diag_MolDiag',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
plt<-(p9+p10+p11+p12+p13+p14)
ggsave(plt,file="allcells.umap.passqc.pdf",height=30,width=50,limitsize=F)

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


saveRDS(dat,file="merged.geneactivity.passqc.SeuratObject.rds")
saveRDS(dat_nonepi,file="merged.geneactivity.passqc.nonepi.SeuratObject.rds")


Idents(dat)<-dat$sample
DefaultAssay(dat)<-"RNA"
p10<-DotPlot(subset(dat,cells=names(dat$assigned_celltype=="basal_epithelial")),features=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2"),cluster.idents=FALSE,dot.scale=8)+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))+
  theme(axis.text.x = element_text(angle=90))
p11<-FeaturePlot(dat,features=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2"),reduction = "allcells.wnn.umap",raster=T)
ggsave(p10/p11,file="basal.features_by_sample.pdf",limitsize=F,height=30,width=30)


