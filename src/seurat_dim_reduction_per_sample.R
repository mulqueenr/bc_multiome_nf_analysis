library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(GenomicRanges)
args = commandArgs(trailingOnly=TRUE)

peaks=read.csv(file=args[1],sep="\t",col.names=c("chr","start","end"))
peaks<-peaks[peaks$chr %in% c(paste0("chr",1:22),"chrX"),]
peaks<-peaks[peaks$start>0,]
peaks<-makeGRangesFromDataFrame(peaks)

outname<-args[2]
nf_dir=getwd()
wd=paste0(nf_dir,"/",outname,"/","outs")

outdir=args[3]
system(paste0("mkdir -p ",outdir))


# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

counts <- Read10X_h5(paste0(wd,"/filtered_feature_bc_matrix.h5")) #count data
fragpath <- paste0(wd,"/atac_fragments.tsv.gz") #atac fragments
metadata_cellranger<-read.csv(paste0(wd,"/per_barcode_metrics.csv")) #metadata
row.names(metadata_cellranger)<-metadata_cellranger$barcode
soupx_output<-readRDS(paste0(wd,"/soupx_corrected_counts.rds")) #load SoupX contamination corrected output
scrublet_output<-read.table(paste0(wd,"/",outname,".scrublet.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
#clean up scrublet output to add to metadata columns
#just a hold over from a python output that I'm correcting.
if(startsWith(scrublet_output$cellid[1],"b")){ 
scrublet_output$cellID<-unlist(lapply(scrublet_output$cellid, function(x) substr(x,2,nchar(x))))}
row.names(scrublet_output)<-scrublet_output$cellID
scrublet_output<-scrublet_output[,c("doublet_scores","predicted_doublets")]

# create a Seurat object containing the RNA data
dat <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
dat[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
#Create corrected RNA data and add to object
dat[["SoupXRNA"]]<-CreateAssayObject(
  counts=soupx_output)

#QC cells
DefaultAssay(dat) <- "ATAC"
dat <- NucleosomeSignal(dat)
dat <- TSSEnrichment(dat)
dat<-AddMetaData(dat,metadata=metadata_cellranger)
dat<-AddMetaData(dat,metadata=scrublet_output)

plt<-VlnPlot(
  object = dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

ggsave(plt,file=paste0(outdir,"/",outname,".qc.pdf"))

# quantify counts in each peak (using merged peak set)
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = dat@assays$ATAC@fragments,
  annotation = annotation
)

#set up colors for samples
my_cols = brewer.pal(1,"Spectral")
alpha_val=0.33

#RNA Processing
DefaultAssay(dat) <- "SoupXRNA"
dat <- SCTransform(dat)
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)
p1<-DimPlot(dat,reduction="rna_umap")+ggtitle("RNA UMAP")

#DNA Accessibility processing
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)
p2<-DimPlot(dat,reduction="atac_umap")+ggtitle("ATAC UMAP")


# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #I think the ATAC UMAP does a better job integrating samples, maybe skip dim 1 for RNA also?
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat <- RunUMAP(
  object = dat,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE
)
p3<-DimPlot(dat,
  reduction="multimodal_umap",
  group.by="predicted_doublets")+ggtitle("Multimodal UMAP Doublets")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-FeaturePlot(dat,
  reduction="multimodal_umap",
  features="doublet_scores")+ggtitle("Multimodal UMAP Scublet Scores")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file=paste0(outdir,"/",outname,".umap.pdf"))
table(dat$predicted_doublets)


saveRDS(dat,file=paste0(outname,".SeuratObject.rds"))