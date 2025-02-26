#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(reshape2)
library(optparse)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
set.seed(123)

option_list = list(
  make_option(c("-s", "--sample_dir"), type="character", default="NAT_1", 
              help="Sample directory from cellranger output.", metavar="character"),
    make_option(c("-p", "--peaks_bed"), type="character", default="NULL", 
              help="Bed file formated peaks.", metavar="character"),
        make_option(c("-o", "--output_directory"), type="character", default=NULL, 
              help="Output directory, defined in nextflow parameters.", metavar="character")

); 

#######testing#########
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
#seurat_obj_list="IDC_12.SeuratObject.rds IDC_4.SeuratObject.rds DCIS_2.SeuratObject.rds IDC_7.SeuratObject.rds IDC_2.SeuratObject.rds IDC_5.SeuratObject.rds IDC_9.SeuratObject.rds IDC_6.SeuratObject.rds NAT_4.SeuratObject.rds NAT_11.SeuratObject.rds IDC_3.SeuratObject.rds ILC_1.SeuratObject.rds IDC_10.SeuratObject.rds IDC_1.SeuratObject.rds DCIS_1.SeuratObject.rds IDC_11.SeuratObject.rds DCIS_3.SeuratObject.rds NAT_14.SeuratObject.rds IDC_8.SeuratObject.rds"
#ref_dir="/home/groups/CEDAR/mulqueen/bc_multiome/ref" 
#met=read.csv("sample_metadata.csv",header=T,sep=",")
#outdir<-"/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/plots"
###################

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample_in=opt$sample_dir
#outname<-"NAT_11"
outname<-sample_in
nf_dir=getwd()
wd=paste0(nf_dir,"/",outname,"/","outs")

#peaks=read.csv(file="merged.nf.bed",sep="\t",col.names=c("chr","start","end"))
peaks=read.csv(file=opt$peaks_bed,sep="\t",col.names=c("chr","start","end"))
peaks<-peaks[peaks$chr %in% c(paste0("chr",1:22),"chrX"),]
peaks<-peaks[peaks$start>0,]
peaks<-makeGRangesFromDataFrame(peaks)

#outdir="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis"
outdir=paste0(opt$output_directory,"/plots")
system(paste0("mkdir -p ",outdir))

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

counts <- Read10X_h5(paste0(wd,"/filtered_feature_bc_matrix.h5")) #count data
fragpath <- paste0(wd,"/atac_fragments.tsv.gz") #atac fragments
metadata_cellranger<-read.csv(paste0(wd,"/per_barcode_metrics.csv")) #metadata
row.names(metadata_cellranger)<-metadata_cellranger$barcode
soupx_output<-readRDS(paste0(wd,"/soupx_corrected_counts.rds")) #load SoupX contamination corrected output
#scrublet is broken due to an annoy versioning error interaction with docker/singularity. I'm just skipping for now
scrublet_output<-read.table(paste0(wd,"/","scrublet_results.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
scrublet_output$Barcode<-substr(scrublet_output$Barcode,2,nchar(scrublet_output$Barcode))
row.names(scrublet_output)<-scrublet_output$Barcode

# create a Seurat object containing the RNA data
dat <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
dat[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

#Create corrected RNA data and add to object
dat[["SoupXRNA"]]<-CreateAssayObject(
  counts=soupx_output)

#QC cells
DefaultAssay(dat) <- "ATAC"
dat<-AddMetaData(dat,metadata=metadata_cellranger)
dat<-AddMetaData(dat,metadata=scrublet_output)

# quantify counts in each peak (using merged peak set)
#i actually used macs3
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object

peaks_assay <-CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(dat),
  annotation = annotation,
  min.features=-1
)

peaks_assay<-subset(peaks_assay, cells=colnames(dat))
dat[["peaks"]]<-peaks_assay
DefaultAssay(dat)<-"peaks"

#set up basic filters
dat<-subset(dat, nCount_RNA>=1000)
dat<-subset(dat, nCount_peaks>=1000)

dat <- NucleosomeSignal(dat)
dat <- TSSEnrichment(dat)


plt<-VlnPlot(
  object = dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

ggsave(plt,file=paste0(outdir,"/",outname,".qc.pdf"))



# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #read in data
  outname=strsplit(x,"[.]")[[1]][1]
  dat<-readRDS(x)
  dat$sample<-outname #set up sample metadata
  dat <- SCTransform(dat,vst.flavor = 'v1')
  print(paste("Finished sample:",outname))
  return(dat)}

seurat_obj_list=strsplit(seurat_obj_list," ")[[1]]
out<-lapply(unlist(seurat_obj_list),merge_seurat)
sample_names=unlist(lapply(strsplit(seurat_obj_list,"[.]"),"[",1))
dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = sample_names, project = "all_data")

dat$cellID<-row.names(dat@meta.data)
dat_met<-as.data.frame(dat@meta.data)
met_merged<-merge(dat_met,met,by.x="sample",by.y="Manuscript_Name",all.x=T)
row.names(met_merged)<-met_merged$cellID
dat<-AddMetaData(dat,met_merged[,c("phase","Original_Sample","sample_ID","sample_weight","Diagnosis","Mol_Diagnosis","sampled_site","batch","outcome")])

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
#opt$object_input="merged.geneactivity.SeuratObject.rds"

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

#define scrublet score based on cell counts per sample (doublets less likely with lower cell counts)

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
