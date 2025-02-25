library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
nf_dir=getwd()

outname=args[1]
wd=paste0("./",args[2],"/","outs")
outdir=args[3]

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

setwd(wd)
counts <- Read10X_h5("filtered_feature_bc_matrix.h5") #count data
fragpath <- "atac_fragments.tsv.gz" #atac fragments
metadata_cellranger<-read.csv("per_barcode_metrics.csv") #metadata
row.names(metadata_cellranger)<-metadata_cellranger$barcode
soupx_output<-readRDS("soupx_corrected_counts.rds") #load SoupX contamination corrected output
scrublet_output<-read.table(paste0(outname,".scrublet.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
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

setwd(nf_dir)#save to working directory for nextflow
ggsave(plt,file=paste0(outdir,"/",outname,".qc.pdf"))
saveRDS(dat,file=paste0(outname,".SeuratObject.rds"))