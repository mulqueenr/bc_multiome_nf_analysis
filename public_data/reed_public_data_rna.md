
### Reed HBCA
https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637


```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref = "${proj_dir}/ref"
mkdir -p $ref
```

```bash
mkdir -p ${ref}/reed
cd ${ref}/reed
wget https://datasets.cellxgene.cziscience.com/0b158a5d-74aa-49d0-99ff-4aa1173dbeea.rds
mv 0b158a5d-74aa-49d0-99ff-4aa1173dbeea.rds stroma.rds

wget https://datasets.cellxgene.cziscience.com/97e96fb1-8caf-4f08-9174-27308eabd4ea.rds
mv 97e96fb1-8caf-4f08-9174-27308eabd4ea.rds immune.rds

wget https://datasets.cellxgene.cziscience.com/4c2f72fb-81f0-4215-8409-143c8932a15d.rds
mv 4c2f72fb-81f0-4215-8409-143c8932a15d.rds epithelial.rds

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
library(stringr)
library(plyr)
library(org.Hs.eg.db)

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
setwd(paste0(proj_dir,"/ref/reed/"))

preprocess_reed<-function(x){
  #CONVERTING ENSEMBL TO GENE SYMBOLS
  met<-x@meta.data
  gene_names<-select(org.Hs.eg.db, keys = row.names(x@assays$RNA), keytype = 'ENSEMBL', columns = 'SYMBOL')
  genes.filter <- gene_names[!is.na(gene_names$SYMBOL),]$ENSEMBL 
  counts <- x@assays$RNA@counts[genes.filter,]
  row.names(counts)<-gene_names[!is.na(gene_names$SYMBOL),]$SYMBOL
  counts<-counts[!duplicated(row.names(counts)),]
  x <- CreateSeuratObject(counts=counts)
  x<-AddMetaData(x,metadata=met)
  DefaultAssay(x)<-"RNA"

  #PREPROCESSING
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x)
  x<-ScaleData(x)
  return(x)
}

immune<-readRDS("immune.rds")
immune<-preprocess_reed(x=immune)
saveRDS(immune,file="immune.rds")

stromal<-readRDS("stroma.rds")
stromal<-preprocess_reed(x=stromal)
saveRDS(stromal,file="stroma.rds")


epithelial<-readRDS("epithelial.rds")
epithelial<-preprocess_reed(x=epithelial)
saveRDS(epithelial,file="epithelial.rds")
```