
## Download HTAN multiome data Terekhanova et al.
Selected data using the HTAN portal. 
Downloading level 3 and 4 data (open source) from synapse


```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref = "${proj_dir}/ref"
mkdir -p $ref
```

```bash
pip install --upgrade synapseclient
#login with authentification token

cd $ref
mkdir -p ${ref}/terekhanova
mkdir -p ${ref}/terekhanova/rna_obj
mkdir -p ${ref}/terekhanova/atac_frag

#rna
cd ${ref}/terekhanova/rna_obj
for i in "syn53214887 syn53214703 syn53214932 syn53214720 syn53214931 syn53214683 syn53214973 syn53214802 syn53214920 syn53214695 syn53214950 syn53214812 syn53214757 syn53214656 syn53214605 syn53214771 syn53214677 syn53214594";
do synapse get $i; done

#atac
cd ${ref}/terekhanova/atac_frag
for i in syn52176627 syn52175906 syn52175971 syn52176152 syn52176835 syn52175918 syn52175960 syn52176281 syn52176599 syn52175948 syn52175949 syn52176099 syn52176830 syn52175958 syn52175962 syn52176354 syn52176791 syn52175969 syn52175982 syn52176198 syn52176842 syn52175970 syn52175997 syn52176397 syn52176839 syn52175993 syn52176003 syn52176217 syn52176825 syn52175992 syn52175991 syn52176141 syn52176769 syn52176004 syn52176005 syn52176145 syn52176693 syn52176006 syn52176012 syn52176156 syn52176747 syn52176007 syn52176010 syn52176140 syn52176633 syn52176008 syn52176015 syn52176100 syn53215775 syn53215789 syn53215774 syn53215798 syn53215811 syn53215784 syn53214546 syn53214579 syn53214453 syn53214571 syn53214519 syn53214649 syn53214791 syn53214450 syn53214493 syn53214547 syn53214543 syn53214471; 
do synapse get $i; done

#index fragment files
for i in *atac_fragments.tsv.gz;
do tabix -p bed $i; done

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
library(Seurat)
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
library(EnsDb.Hsapiens.v86)
library(Signac)

#read in RNA object, use atac_fragments.tsv.gz to make a feature matrix using our peaks

# get gene annotations for hg38
annot<- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annot),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annot) <- ucsc.levels #standard seq level change threw error, using a string replace instead

peaks<- as.data.frame(read.table("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/peaks/merged.nf.bed",header=F,sep="\t",col.names=c("chr","start","end")))
peaks<-makeGRangesFromDataFrame(peaks)

setwd("/home/groups/CEDAR/mulqueen/bc_multiome/ref/terekhanova")
rna_objs<-list.files("./rna_obj",pattern=".rds$",full.names=TRUE,include.dirs=TRUE)
atac_fragments<-list.files("./atac_frag",pattern="fragments.tsv.gz$",full.names=TRUE,include.dirs=TRUE)

rna_objs<-setNames(rna_objs,nm=gsub(".rds","",basename(rna_objs)))
atac_fragments<-setNames(atac_fragments,nm=gsub("-atac_fragments.tsv.gz$","",basename(atac_fragments)))

rna_objs<-rna_objs[names(rna_objs) %in% names(atac_fragments)]
atac_fragments<-atac_fragments[names(atac_fragments) %in% names(rna_objs)]


build_seurat_obj<-function(name_in){
  obj<-readRDS(rna_objs[name_in])
  fragments <- CreateFragmentObject(
    path = unname(atac_fragments[name_in]),
    cells=unique(obj$Original_barcode),
    validate.fragments = FALSE)
  obj$cell_id<-Cells(obj)
  obj<-RenameCells(obj,new.names=obj$Original_barcode) #reset to old names for fragment files
  counts<-FeatureMatrix(fragments=fragments,features=peaks)
  obj[["our_peaks"]] <- CreateChromatinAssay(
      sep = c(":", "-"),
      counts=counts,
      fragments = fragments,
      annotation = annot)
  obj<-RenameCells(obj,new.names=obj$cell_id) #reset to new names now that fragments are calculated
  return(obj)
}

obj_list<-lapply(names(rna_objs),build_seurat_obj)

combined <- merge(obj_list[[1]], y = as.list(obj_list[2:length(out)]), project = "terekhanova")

#merge all objects and then preprocess a bit
#cell labels in $cell_type
combined<-NormalizeData(combined)
combined<-FindVariableFeatures(combined)
combined<-ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30, reduction = 'pca',reduction.name="RNA_umap")

plt<-DimPlot(combined,group.by="cell_type",reduction="RNA_umap")+ggtitle("RNA")
ggsave(plt,file="terekhanova_multiome.atac_umap.pdf")
saveRDS(combined,file="terekhanova_multiome.rds")
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
combined<-readRDS(file=paste0(proj_dir,"/ref/nakshatri/","nakshatri_multiome.rds"))

run_top_TFs(combined,prefix="ref_cell_type",i="author_cell_type",n_markers=5) #generate top TF markers per cell type


```
