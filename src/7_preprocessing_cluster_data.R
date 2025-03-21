#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(plyr)
library(optparse)


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects/merged.geneactivity.SeuratObject.rds", 
              help="Input seurat object", metavar="character"),
    make_option(c("-r", "--ref_dir"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref", 
              help="Reference directory containing genome information. default: %default]", metavar="character"),
    make_option(c("-o","--plot_output_directory"), type="character", default=NULL,
              help="Directory to publish output plots to.", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
ref_dir=opt$ref_dir
met<-read.csv(opt$metadata,header=T,sep=",")
outdir<-opt$plot_output_directory
dat<-readRDS(file=opt$object_input)

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

#define scrublet score based on cell counts per sample (doublets less likely with lower cell counts)


saveRDS(dat,file="merged.geneactivity.SeuratObject.rds")

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