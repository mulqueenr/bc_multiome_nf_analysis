
## Public Datasets for Comparison 
<!-- Done -->

### Using Transfer Anchors for Cell identification.
<!-- Done -->

Using Swarbrick paper labels for transfer. https://pubmed.ncbi.nlm.nih.gov/34493872/

Download data
```bash
cd /home/groups/CEDAR/mulqueen/ref/swarbrick
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
tar -xvf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
```
Make Seurat Object with Metadata

```R
library(Seurat)

setwd("/home/groups/CEDAR/mulqueen/ref/swarbrick")
counts<-ReadMtx(mtx="count_matrix_sparse.mtx",cells="count_matrix_barcodes.tsv",features="count_matrix_genes.tsv",feature.column=1) #sparse matrix of counts
metadata<-read.csv("metadata.csv") #metadata
row.names(metadata)<-metadata$X
# create a Seurat object containing the RNA adata
swarbrick <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)
swarbrick<-AddMetaData(swarbrick,metadata=metadata)
saveRDS(swarbrick,"/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
```

### Using EMBO paper for transfer of signatures as well. 
<!-- Done -->

https://doi.org/10.15252/embj.2020107333
Full code here: https://www.nature.com/articles/s41597-022-01236-2 data available here https://doi.org/10.6084/m9.figshare.17058077

Download data from GEO FTP server

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
wget https://figshare.com/ndownloader/articles/17058077/versions/1
unzip 1
```

Set up cell types by seurat cluster ID based on main figures.

```R
library(Seurat)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/ref/embo")
#match suerat clusters to assigned cell types in Fig 7C
##ER+ nonepi celltypes##
dat<-readRDS("SeuratObject_ERTotalSub.rds") #ER+ tumor non-epithelial cells
er_nonepi<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("T cells","TAMs","CAFs","Pericytes","NA","Endothelial","TAMs_2","B cells","Myeloid","CAFs","Plasma cells","NA","NA"))
er_nonepi_cells<-setNames(names(er_nonepi[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_nonepi_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalSub.umap.pdf")
system("slack -F ERTotalSub.umap.pdf ryan_todo")
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
system("slack -F ERTotalTC.umap.pdf ryan_todo")
saveRDS(dat,file="SeuratObject_ERTotalTC.rds") #overwrite with cell types added to metadata


#match suerat clusters to assigned cell types in Fig 6E
dat<-readRDS("SeuratObject_ERTotalTum.rds") #ER+ tumor epithelial
er_epi<-setNames(
  seq(0,max(as.numeric(unique(dat$seurat_clusters))))
  ,nm=c("epithelial","cycling epithelial","epithelial"))
er_epi_cells<-setNames(names(er_epi[dat$seurat_clusters]),nm=names(dat$seurat_clusters))
dat<-AddMetaData(dat,er_epi_cells,col.name="celltype")
plt<-DimPlot(dat,group.by="celltype")
ggsave(plt,file="ERTotalTum.umap.pdf")
system("slack -F ERTotalTum.umap.pdf ryan_todo")
saveRDS(dat,"SeuratObject_ERTotalTum.rds")

#ER+ All Cells
dat1<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalSub.rds") #ER+ tumor non-epithelial cells
dat2<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalTum.rds") #ER+ tumor epithelial
dat_tc<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERTotalTC.rds") #ER+ tumor T-cells
dat<-merge(dat1,dat2)
dat<-AddMetaData(dat,dat_tc$celltype,col.name="TCell_Subtype")
saveRDS(dat,"SeuratObject_ERProcessed.rds")
```

### Using PBMC Data set for Immune Cell Subtyping
Files downloaded from UCSC Cell Browser (https://cells.ucsc.edu/?ds=multimodal-pbmc+sct) and this manuscript https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3
These will be used later for higher resolution immune cell subtyping.

```bash
mkdir /home/groups/CEDAR/mulqueen/ref/hao
```

```R
library(Seurat)
library(ggplot2)
library(data.table)

setwd("/home/groups/CEDAR/mulqueen/ref/hao")

metadata<-fread("https://cells.ucsc.edu/multimodal-pbmc/sct/meta.tsv") #download metadata
metadata<-as.data.frame(metadata)
row.names(metadata)<-metadata$V1

mat<-fread("https://cells.ucsc.edu/multimodal-pbmc/sct/exprMatrix.tsv.gz") #download counts
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
genes_list<-which(!duplicated(genes)) #remove duplicate gene names
mat = data.frame(mat[genes_list,-1], row.names=genes[genes_list])

hao <- CreateSeuratObject(
  counts = mat,
  assay = "RNA"
)
hao<-AddMetaData(hao,metadata=metadata)

saveRDS(hao,"hao.SeuratObject.Rds")

```
## Swarbrick Paper Label Transfer
<!-- Done -->

### Transfer Swarbrick cell types
<!-- Done -->

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Using Label transfer to label cell types by Swarbrick paper
#seurat object made by AD
swarbrick<-readRDS("/home/groups/CEDAR/mulqueen/ref/swarbrick/swarbrick.SeuratObject.Rds")
swarbrick<-NormalizeData(swarbrick)
swarbrick<-FindVariableFeatures(swarbrick)
swarbrick<-ScaleData(swarbrick)

##########Apply to single samples as well##################

single_sample_label_transfer<-function(x){
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat$sample<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat$sample<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
  file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  dat<-readRDS(file_in)
  dat@assays$peaks<-dat@assays$ATAC
  dat$sample<-paste0(x)
  }
  DefaultAssay(dat)<-"SoupXRNA"
  dat<-NormalizeData(dat)
  dat<-FindVariableFeatures(dat)
  dat<-ScaleData(dat)
  saveRDS(dat,file=file_in)

  transfer.anchors <- FindTransferAnchors(
    reference = swarbrick,
    reference.assay="RNA",
    query = dat,
    query.assay="SoupXRNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = swarbrick$celltype_major,
  )

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=file_in)
  plt1<-FeaturePlot(dat,features=c('prediction.score.Endothelial','prediction.score.CAFs','prediction.score.PVL','prediction.score.B.cells','prediction.score.T.cells','prediction.score.Myeloid','prediction.score.Normal.Epithelial','prediction.score.Plasmablasts','prediction.score.Cancer.Epithelial'),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by='predicted.id',pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_label_transfer)
#
```

### Transfer EMBO Cell Types Per Sample
<!-- Done -->


```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Using Label transfer to label cell types by Embo Paper
#seurat object made by AD
embo_er<-readRDS("/home/groups/CEDAR/mulqueen/ref/embo/SeuratObject_ERProcessed.rds")
DefaultAssay(embo_er)<-"RNA"
embo_er<-NormalizeData(embo_er)
embo_er<-FindVariableFeatures(embo_er)
embo_er<-ScaleData(embo_er)

##########Apply to single samples as well##################

single_sample_label_transfer<-function(x){
  if(x %in% 1:12){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
    dat$sample<-paste0("sample_",x)
  }else if(x %in% 13:20){
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
    dat$sample<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    out_plot<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".predictions.umap.pdf")
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
    dat<-readRDS(file_in)
    dat@assays$peaks<-dat@assays$ATAC
    dat$sample<-paste0(x)
  }
  DefaultAssay(dat)<-"SoupXRNA"

  transfer.anchors <- FindTransferAnchors(
    reference = embo_er,
    reference.assay="RNA",
    query = dat,
    query.assay="SoupXRNA",
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = embo_er$celltype,
  )
  colnames(predictions)<-paste0("EMBO_",colnames(predictions))

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=file_in)
  plt1<-FeaturePlot(dat,features=c(                     
  "EMBO_prediction.score.Endothelial",       
  "EMBO_prediction.score.TAMs",              
  "EMBO_prediction.score.Pericytes",         
  "EMBO_prediction.score.CAFs",              
  "EMBO_prediction.score.T.cells",           
  "EMBO_prediction.score.Plasma.cells",      
  "EMBO_prediction.score.TAMs_2",            
  "EMBO_prediction.score.B.cells",           
  "EMBO_prediction.score.Myeloid",           
  "EMBO_prediction.score.epithelial",        
"EMBO_prediction.score.cycling.epithelial"),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by='EMBO_predicted.id',pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  system(paste0("slack -F ",out_plot," ryan_todo"))
  }

lapply(c(1,3,4,5,6,7,8,9,10,11,12,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),single_sample_label_transfer)


```
