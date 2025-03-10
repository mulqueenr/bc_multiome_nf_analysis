
### Using EMBO paper for transfer of signatures as well. 
Pal et al. https://doi.org/10.15252/embj.2020107333

Full code here: https://www.nature.com/articles/s41597-022-01236-2 
Data available here https://doi.org/10.6084/m9.figshare.17058077

Download data from GEO FTP server


```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref = "${proj_dir}/ref"
mkdir -p $ref
```

```bash
mkdir -p ${ref}/embo
cd ${ref}/embo
wget https://figshare.com/ndownloader/articles/17058077/versions/1
unzip 1
```


Using common SIF for processing
```bash
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"

singularity shell \
--bind /home/groups/CEDAR/mulqueen/bc_multiome \
--bind /home/groups/CEDAR/mulqueen/ref \
--bind /home/groups/CEDAR/mulqueen/src/miniconda3/bin \
$sif
```

Set up cell types by seurat cluster ID based on main figures.

```R
library(Seurat)
library(ggplot2)
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
setwd(paste0(proj_dir,"/ref/embo"))

#match seurat clusters to assigned cell types in Fig 7C
##ER+ nonepi celltypes##
dat<-readRDS("SeuratObject_ERTotalSub.rds") #ER+ tumor non-epithelial cells
er_nonepi<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("T cells","TAMs","CAFs","Pericytes","NA","Endothelial","TAMs_2","B cells","Myeloid","CAFs","Plasma cells","NA","NA"))
er_nonepi_cells<-setNames(names(er_nonepi[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_nonepi_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalSub.umap.pdf")
saveRDS(dat,file="SeuratObject_ERTotalSub.rds") #overwrite with cell types added to metadata

#match seurat clusters to assigned cell types in Fig EV4
dat<-readRDS("SeuratObject_ERTotalTC.rds") #ER+ tumor T-cells
er_nonepi_tcells<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("CD8+ effector","naive/resting","Treg","plasma","NK","NA"))
er_nonepi_tcells_cells<-setNames(names(er_nonepi_tcells[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_nonepi_tcells_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalTC.umap.pdf")
saveRDS(dat,file="SeuratObject_ERTotalTC.rds") #overwrite with cell types added to metadata

#match seurat clusters to assigned cell types in Fig 6E
dat<-readRDS("SeuratObject_ERTotalTum.rds") #ER+ tumor epithelial
er_epi<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("epithelial","cycling epithelial","epithelial"))
er_epi_cells<-setNames(names(er_epi[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_epi_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalTum.umap.pdf")
saveRDS(dat,"SeuratObject_ERTotalTum.rds")

#ER+ All Cells
dat1<-readRDS("SeuratObject_ERTotalSub.rds") #ER+ tumor non-epithelial cells
dat2<-readRDS("SeuratObject_ERTotalTum.rds") #ER+ tumor epithelial
dat_tc<-readRDS("SeuratObject_ERTotalTC.rds") #ER+ tumor T-cells
dat<-merge(dat1,dat2)
dat<-AddMetaData(dat,dat_tc$celltype,col.name="TCell_Subtype")

DefaultAssay(dat)<-"RNA"
dat<-NormalizeData(dat)
dat<-FindVariableFeatures(dat)
dat<-ScaleData(dat)
saveRDS(dat,"SeuratObject_ERProcessed.rds")

```
