
Wu et al. https://pubmed.ncbi.nlm.nih.gov/34493872/
Setting up seurat objects here:

```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref = "${proj_dir}/ref"
mkdir -p $ref
```
### Using Transfer Anchors for Cell identification.

Using Swarbrick paper labels for transfer. 

Download data
```bash
mkdir -p ${ref}/swarbrick
cd ${ref}/swarbrick
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
tar -xvf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
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

Make Seurat Object with Metadata
```R
library(Seurat)
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
setwd(paste0(proj_dir,"/ref/swarbrick"))

counts<-ReadMtx(mtx="count_matrix_sparse.mtx",cells="count_matrix_barcodes.tsv",features="count_matrix_genes.tsv",feature.column=1) #sparse matrix of counts
metadata<-read.csv("metadata.csv") #metadata
row.names(metadata)<-metadata$X
# create a Seurat object containing the RNA adata
swarbrick <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)
swarbrick<-AddMetaData(swarbrick,metadata=metadata)
DefaultAssay(swarbrick)<-"RNA"
swarbrick<-NormalizeData(swarbrick)
swarbrick<-FindVariableFeatures(swarbrick)
swarbrick<-ScaleData(swarbrick)
swarbrick$celltype<-swarbrick$celltype_major

saveRDS(swarbrick,"swarbrick.SeuratObject.Rds")
```