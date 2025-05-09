

### Download of Multiome Normal Breast Data from Cellxgene
Bhat-Nakshatri et al.

From https://cellxgene.cziscience.com/e/61af564d-e5ea-4d34-a0f3-2668a00db376.cxg/
https://pubmed.ncbi.nlm.nih.gov/39122969/

```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref="${proj_dir}/ref"
mkdir -p $ref
```

```bash
mkdir -p ${ref}/nakshatri_multiome
cd ${ref}/nakshatri_multiome
wget https://datasets.cellxgene.cziscience.com/63a485bc-cac7-49d2-83ed-8e07ca4efa2a.rds

#curl download of fastq files
curl --location --fail https://service.azul.data.humancellatlas.org/manifest/files/ksQwlKVkY3A0NKRjdXJsxBDVqi-KPW5cZIBunIwUGnbJxBDNWgLvW59Vl54jxXbB21J3xCAWcw-YKqAgo8jjf-6V_OLtqiQl2o01JAye0IwhVslyhg | curl --fail-early --continue-at - --retry 15 --retry-delay 10 --config -

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244585

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244585/suppl/GSE244585_RAW.tar
tar -xvf GSE244585_RAW.tar
gzip -d *tbi.gz
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
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(ggplot2)
set.seed(1234)
library(BiocParallel)
library(universalmotif)
library(GenomicRanges)

peaks<- as.data.frame(read.table("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/peaks/merged.nf.bed",header=F,sep="\t",col.names=c("chr","start","end")))
peaks<-makeGRangesFromDataFrame(peaks)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri_multiome")

ref<-readRDS("63a485bc-cac7-49d2-83ed-8e07ca4efa2a.rds")
#added cell id names are pools 1-16
#assuming cellid names are [STD 10x Barcodes-1]_[Pool1:16]
#levels(ref@meta.data$Pool)
# [1] "Pool1"  "Pool3"  "Pool4"  "Pool5"  "Pool6"  "Pool7"  "Pool8"  "Pool9" 
# [9] "Pool11" "Pool12" "Pool13" "Pool14" "Pool15" "Pool17" "Pool23" "Pool25"
#Plan is split by pool, rename cells to fit existing fragment files, then recombine object
pool_frag<-setNames(nm=levels(ref@meta.data$Pool),list.files(pattern="fragments.tsv.gz$")) #ordered by GSM###

ref<-SplitObject(ref, split.by = "Pool")
build_atac<-function(x)  {
  test<-ref[[x]]
  working_pool<-test@meta.data$Pool[1]
  fragpath<-pool_frag[working_pool]
  test<-RenameCells(test, new.names = unlist(lapply(strsplit(Cells(test),"_"),"[",1)))
  fragments <- CreateFragmentObject(
    path = unname(fragpath),
    cells = colnames(test),
    validate.fragments = FALSE
  )
  #start with our set of peaks, also call peaks on merged set later
  counts<-FeatureMatrix(fragments=fragments,features=peaks)
  # create ATAC assay and add it to the object
  test[["our_peaks"]] <- CreateChromatinAssay(
    sep = c(":", "-"),
    counts=counts,
    fragments = fragpath,
    annotation = annotation
  )
  return(test)
}

frag_fixed_ref<-lapply(1:length(ref),build_atac)

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = frag_fixed_ref[[1]],
  y =unlist(frag_fixed_ref[2:length(frag_fixed_ref)]),
  add.cell.ids = names(pool_frag)
)
combined[["our_peaks"]]
DefaultAssay(combined)<-"our_peaks"

#filter peaks not seen in data
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
main.chroms <- main.chroms[!(main.chroms %in% c("chrY","chrM"))] 
keep.peaks <- which(as.character(seqnames(granges(combined[["our_peaks"]]))) %in% main.chroms)
combined[["our_peaks"]] <- subset(combined[["our_peaks"]], features = rownames(combined[["our_peaks"]])[keep.peaks])

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi',reduction.name="our_peaks_umap")

#RNA preprocessing, set counts to gene SYMBOLS
gene_names<-select(org.Hs.eg.db, keys = row.names(combined@assays$RNA), keytype = 'ENSEMBL', columns = 'SYMBOL')
genes.filter <- gene_names[!is.na(gene_names$SYMBOL),]$ENSEMBL 
counts <- combined@assays$RNA@counts[genes.filter,]
row.names(counts)<-gene_names[!is.na(gene_names$SYMBOL),]$SYMBOL
counts<-counts[!duplicated(row.names(counts)),]
combined[["RNA"]]<-NULL
combined[["RNA"]]<-CreateAssayObject(counts=counts)
DefaultAssay(combined)<-"RNA"

#PREPROCESSING RNA
combined<-SCTransform(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30, reduction = 'pca',reduction.name="RNA_umap")


plt1<-DimPlot(combined,group.by="author_cell_type",reduction="our_peaks_umap")+ggtitle("Our Peaks")
plt2<-DimPlot(combined,group.by="author_cell_type",reduction="RNA_umap")+ggtitle("RNA")
ggsave(plt1/plt2,file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.atac_umap.pdf"),width=12)
saveRDS(combined,file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.rds"))
```

```R
# RUN CHROMVAR

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
library(plyr)
library(org.Hs.eg.db)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
library(BiocParallel)
library(universalmotif)
library(GenomicRanges)
register(SerialParam()) #using single core mode

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
combined<-readRDS(file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.rds"))

DefaultAssay(combined)<-"our_peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
main.chroms <- main.chroms[!(main.chroms %in% c("chrY","chrM"))] 
keep.peaks <- which(as.character(seqnames(granges(combined[["our_peaks"]]))) %in% main.chroms)
combined[["our_peaks"]] <- subset(combined[["our_peaks"]], features = rownames(combined[["our_peaks"]])[keep.peaks])

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
peaks<-granges(combined[["our_peaks"]])

motif.matrix.hg38 <- CreateMotifMatrix(
  features = peaks, 
  pwm = pfm, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  use.counts = FALSE)

motif.hg38 <- CreateMotifObject(
  data = motif.matrix.hg38, 
  pwm = pfm)

combined <- SetAssayData(object = combined, 
  assay = 'our_peaks', 
  layer = 'motifs', 
  new.data = motif.hg38)

combined <- RegionStats(object = combined, 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="our_peaks")

combined <- RunChromVAR( object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="our_peaks")

saveRDS(combined,file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.chromvar.rds"))

# compute gene activities
gene.activities <- GeneActivity(combined)

# add the gene activity matrix to the Seurat object as a new assay
combined[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_GeneActivity)
)


saveRDS(combined,file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.geneactivity.rds"))

run_top_TFs(combined,prefix="ref_cell_type",i="author_cell_type",n_markers=5) #generate top TF markers per cell type

```

Call peaks per cell type
```R
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
library(plyr)
library(org.Hs.eg.db)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
library(BiocParallel)
library(universalmotif)
library(GenomicRanges)
register(SerialParam()) #using single core mode

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
combined<-readRDS(file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.rds"))

peaks <- CallPeaks(
  object = combined,
  group.by = "author_cell_type",
  macs2.path = "/home/users/mulqueen/macs3",
  outdir="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri",
  fragment.tempdir="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri"
)

saveRDS(peaks,file="nakshatri_celltype_peaks.rds")
```