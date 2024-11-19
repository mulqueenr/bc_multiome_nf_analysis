
```bash
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell \
--bind /home/groups/CEDAR/scATACcnv/Hisham_data \
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

RNA_cluster_idents_samples<-lapply(
    list.files(path="/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate",recursive=TRUE,pattern="*RNA_clones.csv",full.names=T), 
    function(RNA){
    in_clus=read.csv(RNA,header=F,col.names=c("cell","RNA_cluster"))
    in_clus$sample=gsub(basename(RNA),pattern="_RNA_clones.csv",replace="")
    in_clus$cellID=paste(in_clus$sample,in_clus$cell,sep="_")
    in_clus$cellID<-gsub('\\.', '-', in_clus$cellID)
    row.names(in_clus)<-in_clus$cellID
    in_clus$RNA_cluster<-paste(in_clus$sample,in_clus$RNA_cluster,sep=".")
    return(in_clus)}
)


ATAC_cluster_idents_samples<-lapply(
    list.files(path="/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate",recursive=TRUE,pattern="*ATAC_clones.csv",full.names=T), 
    function(RNA){
    in_clus=read.csv(RNA,header=F,col.names=c("cell","ATAC_cluster"))
    in_clus$sample=gsub(basename(RNA),pattern="_ATAC_clones.csv",replace="")
    in_clus$cellID=paste(in_clus$sample,in_clus$cell,sep="_")
    in_clus$cellID<-gsub('\\.', '-', in_clus$cellID)
    row.names(in_clus)<-in_clus$cellID
    in_clus$ATAC_cluster<-paste(in_clus$sample,in_clus$ATAC_cluster,sep=".")
    return(in_clus)}
)

RNA_cluster_in<-do.call("rbind",RNA_cluster_idents_samples)
ATAC_cluster_in<-do.call("rbind",ATAC_cluster_idents_samples)
RNA_cluster_in<-RNA_cluster_in["RNA_cluster"]
ATAC_cluster_in<-ATAC_cluster_in["ATAC_cluster"]

dat<-AddMetaData(dat,RNA_cluster_in)
dat<-AddMetaData(dat,ATAC_cluster_in)

saveRDS(dat,"merged.clone_annot.passqc.SeuratObject.rds")

dat<-subset(dat,cells=row.names(dat@meta.data[!is.na(dat@meta.data$RNA_cluster) & !is.na(dat@meta.data$ATAC_cluster),]))
p1<-DimPlot(dat, group.by = 'RNA_cluster',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
p2<-DimPlot(dat, group.by = 'ATAC_cluster',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F)
p3<-DimPlot(dat, group.by = 'sample',label = TRUE, repel = TRUE,reduction = "allcells.wnn.umap",raster=F) 
ggsave(p1/p2/p3,file="allcells.umap.clones.passqc.pdf",height=30,width=10,limitsize=F)

