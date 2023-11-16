library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)

peaks=readRDS(file=arg[1])
dat=readRDS(arg[2])
outname<-strsplit(arg[2],"[.]")[1]
outdir=arg[3]

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead


#perform initial clustering, and remove scrublet detected doublets
# call peaks using MACS2
DefaultAssay(dat) <- "ATAC"

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

saveRDS(dat,file=arg[2])