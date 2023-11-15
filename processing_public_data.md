# Public Datasets for Comparison 

### Using Transfer Anchors for Cell identification.

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

https://doi.org/10.15252/embj.2020107333

Full code here: https://www.nature.com/articles/s41597-022-01236-2 
Data available here https://doi.org/10.6084/m9.figshare.17058077

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

### Human Breast Atlas Data

https://www.nature.com/articles/s41586-023-06252-9

```bash
cd /home/groups/CEDAR/mulqueen/ref/hbca
```

Note the curl command below expires after a time. So you will need to get a new https link through the website.

```bash
#https://cellxgene.cziscience.com/collections/4195ab4c-20bd-4cd3-8b3d-65601277e731
curl -o hbca.rds "https://corpora-data-prod.s3.amazonaws.com/cd74e95e-6583-4875-a0ba-f2eae5a1e5a6/local.rds?AWSAccessKeyId=ASIATLYQ5N5XRLQ7SW4J&Signature=m1IH5mfc%2BgYSTndRvMLOcSyxWS0%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEGoaCXVzLXdlc3QtMiJIMEYCIQDjh9nXAWfz302Ako9qR7jh8iGjtuPU2k%2FI4Xj%2BqFMu%2BgIhAPIQ5y90eVnRLYewFD1KlKAzUGcTDpF8VU38gJTPd1A6KvQDCIP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igzd0XhQSR2kjPS77SkqyAMJJL0mRsb2YT4pJMMjwxS8G0m9pcr1ry%2FeZdaKWk574tPGxSk9i0MGf2UiYiaADNKk2d1hKJybMcsUJPHq%2FLQUYQ%2FcxzwN6WncxjrgvLz5gB9445oFqjeenKHlrQcooLTnfmsByab034aMJUG0Wdwpldl%2BHZseec37fHhzJSJdMXLYMucK17Qv5A8x5jMIYNuryJ8NuVRXC9vhY65IHHZfCsYTDrhWxl6Y6nFGIgqP%2BZsuThn5xYc4w8IlEkz9j%2FpnNcMJ32uuuLnLr%2FEgg%2BZmoT43av1Rn9wa%2BqUAzCVCk8OOrhfSBjVSaVvK%2BXYdc0%2BE3OnqDyTMUgRTFhb%2ByvWBNHEDNiGYuT810gY%2Fe9uO4GzVVkp9ZdKmUgDOs4wi%2B2ROpY8c4EtyrkOAZ87CAhLMLbSorCRoLGVWhb2z%2F%2FiBYAUyxW7XwXnvMXakwlszBQfIZQpn4EyP4DHjucmATY6j%2FScT2%2BSGkyvBxu2xb3MQzs%2BHMZp1muaOvu9hoGZmQ%2FzSaQxfpDC1QIRN2m%2BABnPT0Gc%2FG9TDyfw9aWkb5wt8JrQ4L3ZBjF5ullR3B8uNRVXuktHxjtFldodVvHbiQrmbLMHNvcTOVZcwp863qQY6pAHAo98DgLJJ3mmOQ2ZQCtpgj4KnRIbpcq4D1S3Vkx8Hwh%2F651dADNY%2Fi1%2FGYR%2BIYI7op%2BTtLjGSrxG2v%2FZVNYHuCLIzjpuMYJRwhr%2B3xJ3oTJwME6tl4I1y4o0djJjuZSxFk1tip8uHDnSur70Ktvu9DvJeqgQzCl3it9THC%2Ftk3o0E7yEzCwB2nkWVULdZV1vTPYB7zkafoS84leGWalCEfZy2DA%3D%3D&Expires=1698118848"
```