
```bash
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell \
--bind /home/groups/CEDAR/scATACcnv/Hisham_data \
--bind /home/groups/CEDAR/mulqueen/ref \
--bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(Seurat)
library(Signac)
library(ggplot2)
library(optparse)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.passqc.SeuratObject.rds"
dat=readRDS(opt$object_input)

#/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/[sample]_clone_val/[sample]_[ATAC/RNA]_clones.csv
#CNV data processed by TM and stored in directory listed above.

merged_cluster_idents_samples<-lapply(
    list.files(path="/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate",recursive=TRUE,pattern="*ATAC_anno.csv",full.names=T), 
    function(merge){
    in_clus=read.csv(merge,header=T)
    in_clus$sample=gsub(basename(merge),pattern="_ATAC_anno.csv",replace="")
    in_clus$cellID=paste(in_clus$sample,in_clus$X,sep="_")
    in_clus$cellID<-gsub('\\.', '-', in_clus$cellID)
    row.names(in_clus)<-in_clus$cellID
    in_clus$atac_cluster<-paste(in_clus$sample,in_clus$atac_cluster,sep=".")
    in_clus$rna_cluster<-paste(in_clus$sample,in_clus$rna_cluster,sep=".")
    in_clus$merge_cluster<-paste(in_clus$sample,in_clus$merge_cluster,sep=".")
    return(in_clus)}
)

cnv_cluster_in<-do.call("rbind",merged_cluster_idents_samples)
cnv_cluster_in<-cnv_cluster_in[c("atac_cluster","rna_cluster","merge_cluster")]

dat<-AddMetaData(dat,cnv_cluster_in)

dat$ploidy<-"diploid"
diploid_clone<-c("DCIS_03.1","IDC_01.4","IDC_02.1","IDC_05.1","IDC_06.3","IDC_08.2","IDC_09.2","IDC_10.3","IDC_11.2","IDC_12.4","ILC_02.1","ILC_04.2")
dat@meta.data[!(dat$merge_cluster %in% diploid_clone),]$ploidy<-"aneuploid"
dat@meta.data[is.na(dat$merge_cluster),]$ploidy<-NA

dat@meta.data[dat$assigned_celltype=="luminal_epithelial",]$assigned_celltype<-ifelse(dat@meta.data[dat$assigned_celltype=="luminal_epithelial",]$ploidy=="aneuploid","cancer_luminal_epithelial","luminal_epithelial")
saveRDS(dat,"merged.clone_annot.passqc.SeuratObject.rds")
dat<-readRDS("merged.clone_annot.passqc.SeuratObject.rds")

p1<-DimPlot(dat, group.by = 'rna_cluster',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F,na.value=NULL) + theme(legend.position="none") +scale_fill_manual(na.value=NULL) 
p2<-DimPlot(dat, group.by = 'atac_cluster',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F,na.value=NA) + theme(legend.position="none")
p3<-DimPlot(dat, group.by = 'merge_cluster',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F,na.value=NA) + theme(legend.position="none")
p4<-DimPlot(dat, group.by = 'sample',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F,na.value=NA) + theme(legend.position="none")

ggsave((p1+p2)/(p3+p4),file="allcells.umap.clones.passqc.pdf",height=20,width=20,limitsize=F)


