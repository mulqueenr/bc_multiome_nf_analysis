

### Human Breast Atlas Data

https://www.nature.com/articles/s41586-023-06252-9


```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref = "${proj_dir}/ref"
mkdir -p $ref
```

```bash
mkdir -p ${ref}/hbca
cd ${ref}/hbca
```

Note the curl command below expires after a time. So you will need to get a new https link through the website.

```bash
#https://cellxgene.cziscience.com/collections/4195ab4c-20bd-4cd3-8b3d-65601277e731
curl -o hbca.rds "https://corpora-data-prod.s3.amazonaws.com/cd74e95e-6583-4875-a0ba-f2eae5a1e5a6/local.rds?AWSAccessKeyId=ASIATLYQ5N5XRLQ7SW4J&Signature=m1IH5mfc%2BgYSTndRvMLOcSyxWS0%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEGoaCXVzLXdlc3QtMiJIMEYCIQDjh9nXAWfz302Ako9qR7jh8iGjtuPU2k%2FI4Xj%2BqFMu%2BgIhAPIQ5y90eVnRLYewFD1KlKAzUGcTDpF8VU38gJTPd1A6KvQDCIP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1Igzd0XhQSR2kjPS77SkqyAMJJL0mRsb2YT4pJMMjwxS8G0m9pcr1ry%2FeZdaKWk574tPGxSk9i0MGf2UiYiaADNKk2d1hKJybMcsUJPHq%2FLQUYQ%2FcxzwN6WncxjrgvLz5gB9445oFqjeenKHlrQcooLTnfmsByab034aMJUG0Wdwpldl%2BHZseec37fHhzJSJdMXLYMucK17Qv5A8x5jMIYNuryJ8NuVRXC9vhY65IHHZfCsYTDrhWxl6Y6nFGIgqP%2BZsuThn5xYc4w8IlEkz9j%2FpnNcMJ32uuuLnLr%2FEgg%2BZmoT43av1Rn9wa%2BqUAzCVCk8OOrhfSBjVSaVvK%2BXYdc0%2BE3OnqDyTMUgRTFhb%2ByvWBNHEDNiGYuT810gY%2Fe9uO4GzVVkp9ZdKmUgDOs4wi%2B2ROpY8c4EtyrkOAZ87CAhLMLbSorCRoLGVWhb2z%2F%2FiBYAUyxW7XwXnvMXakwlszBQfIZQpn4EyP4DHjucmATY6j%2FScT2%2BSGkyvBxu2xb3MQzs%2BHMZp1muaOvu9hoGZmQ%2FzSaQxfpDC1QIRN2m%2BABnPT0Gc%2FG9TDyfw9aWkb5wt8JrQ4L3ZBjF5ullR3B8uNRVXuktHxjtFldodVvHbiQrmbLMHNvcTOVZcwp863qQY6pAHAo98DgLJJ3mmOQ2ZQCtpgj4KnRIbpcq4D1S3Vkx8Hwh%2F651dADNY%2Fi1%2FGYR%2BIYI7op%2BTtLjGSrxG2v%2FZVNYHuCLIzjpuMYJRwhr%2B3xJ3oTJwME6tl4I1y4o0djJjuZSxFk1tip8uHDnSur70Ktvu9DvJeqgQzCl3it9THC%2Ftk3o0E7yEzCwB2nkWVULdZV1vTPYB7zkafoS84leGWalCEfZy2DA%3D%3D&Expires=1698118848"
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

```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(plyr)
library(patchwork)
library(org.Hs.eg.db)
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
setwd(paste0(proj_dir,"/ref/hbca/"))

#CONVERTING ENSEMBL TO GENE SYMBOLS
gene_names<-select(org.Hs.eg.db, keys = row.names(hbca@assays$RNA), keytype = 'ENSEMBL', columns = 'SYMBOL')
genes.filter <- gene_names[!is.na(gene_names$SYMBOL),]$ENSEMBL 
counts <- hbca@assays$RNA@counts[genes.filter,]
row.names(counts)<-gene_names[!is.na(gene_names$SYMBOL),]$SYMBOL
hbca_sub <- CreateSeuratObject(counts=counts)
hbca_sub<-AddMetaData(hbca_sub,metadata=hbca@meta.data)
DefaultAssay(hbca_sub)<-"RNA"
Idents(hbca_sub)<-hbca_sub$cell_type
hbca_sub$celltype<-hbca_sub$cell_type

#PREPROCESSING
hbca_sub<-NormalizeData(hbca_sub)
hbca_sub<-FindVariableFeatures(hbca_sub)
hbca_sub<-ScaleData(hbca_sub)
hbca<-saveRDS(hbca_sub,file="hbca.rds")

```