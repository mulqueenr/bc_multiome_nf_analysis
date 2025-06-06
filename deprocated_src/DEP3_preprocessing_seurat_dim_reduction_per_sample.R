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
library(optparse)

option_list = list(
  make_option(c("-s", "--sample_dir"), type="character", default="NAT_1", 
              help="Sample directory from cellranger output.", metavar="character"),
    make_option(c("-p", "--peaks_bed"), type="character", default="NULL", 
              help="Bed file formated peaks.", metavar="character"),
        make_option(c("-o", "--output_directory"), type="character", default=NULL, 
              help="Output directory, defined in nextflow parameters.", metavar="character")

); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample_in=opt$sample_dir
#outname<-"NAT_11"
outname<-sample_in
nf_dir=getwd()
wd=paste0(nf_dir,"/",outname,"/","outs")

#######testing######
#module load singularity/3.8.0 #load singularity
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif
#proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
#opt$output_directory="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3"
#opt$sample_dir="IDC_7"
#opt$peaks_bed="/home/groups/CEDAR/mulqueen/bc_multiome/merged.nf.bed"
####################

#peaks=read.csv(file="merged.nf.bed",sep="\t",col.names=c("chr","start","end"))
peaks=read.csv(file=opt$peaks_bed,sep="\t",col.names=c("chr","start","end"))
peaks<-peaks[peaks$chr %in% c(paste0("chr",1:22),"chrX"),]
peaks<-peaks[peaks$start>0,]
peaks<-makeGRangesFromDataFrame(peaks)

#outdir="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis"
outdir=paste0(opt$output_directory,"/plots")
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
#scrublet is broken due to an annoy versioning error interaction with docker/singularity. I'm just skipping for now
scrublet_output<-read.table(paste0(wd,"/","scrublet_results.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
scrublet_output$Barcode<-substr(scrublet_output$Barcode,2,nchar(scrublet_output$Barcode))
row.names(scrublet_output)<-scrublet_output$Barcode

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
dat<-AddMetaData(dat,metadata=metadata_cellranger)
dat<-AddMetaData(dat,metadata=scrublet_output)

# quantify counts in each peak (using merged peak set)
#i actually used macs3
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat),
  features = peaks,
  cells = colnames(dat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object

peaks_assay <-CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(dat),
  annotation = annotation,
  min.features=-1
)

peaks_assay<-subset(peaks_assay, cells=colnames(dat))
dat[["peaks"]]<-peaks_assay
DefaultAssay(dat)<-"peaks"

#set up basic filters
dat<-subset(dat, nCount_RNA>=1000)
dat<-subset(dat, nCount_peaks>=1000)

dat <- NucleosomeSignal(dat)
dat <- TSSEnrichment(dat)


plt<-VlnPlot(
  object = dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

ggsave(plt,file=paste0(outdir,"/",outname,".qc.pdf"))


if(ncol(dat)>200){

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
  group.by="scrublet_DropletType")+ggtitle("Multimodal UMAP Doublets")

#Cluster on multimodal graph
dat <- FindClusters(dat, resolution = 0.8, verbose = FALSE,graph="wknn")
p4<-FeaturePlot(dat,
  reduction="multimodal_umap",
  features="scrublet_Scores")+ggtitle("Multimodal UMAP Scublet Scores")

#Finally Plot results
plt<-(p1 | p2)/(p3 | p4)
ggsave(plt,file=paste0(outdir,"/",outname,".umap.pdf"))
table(dat$scrublet_DropletType)
}

saveRDS(dat,file=paste0(outname,".SeuratObject.rds"))