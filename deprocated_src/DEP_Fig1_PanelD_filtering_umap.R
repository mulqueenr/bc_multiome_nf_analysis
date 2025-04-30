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
    dims.list = list(1:30, 2:30),
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

saveRDS(dat,file="merged.geneactivity.passqc.SeuratObject.rds")
saveRDS(dat_nonepi,file="merged.geneactivity.passqc.nonepi.SeuratObject.rds")


Idents(dat)<-dat$sample
DefaultAssay(dat)<-"RNA"
p10<-DotPlot(subset(dat,cells=names(dat$assigned_celltype=="basal_epithelial")),features=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2"),cluster.idents=FALSE,dot.scale=8)+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))+
  theme(axis.text.x = element_text(angle=90))
p11<-FeaturePlot(dat,features=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2"),reduction = "allcells.wnn.umap",raster=T)
ggsave(p10/p11,file="basal.features_by_sample.pdf",limitsize=F,height=30,width=30)


