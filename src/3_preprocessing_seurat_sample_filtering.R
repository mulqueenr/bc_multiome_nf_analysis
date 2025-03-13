#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(stringr)
library(reshape2)
library(optparse)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(EnsDb.Hsapiens.v86)
set.seed(123)

option_list = list(
  make_option(c("-s", "--sample_dir"), type="character", default=".", 
              help="Sample directory from cellranger output.", metavar="character"),
  make_option(c("-p", "--peaks_bed"), type="character", default="merged.nf.bed", 
              help="Bed file formated peaks.", metavar="character"),
  make_option(c("-m","--metadata"), type="character", default="sample_metadata.csv",
              help="Comma separated (CSV) metadata file of cell information to be appended.", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample_in=list.dirs(opt$sample_dir,full.names=FALSE,recursive=FALSE)

#peaks=read.csv(file="merged.nf.bed",sep="\t",col.names=c("chr","start","end"))
peaks=read.csv(file=opt$peaks_bed,sep="\t",col.names=c("chr","start","end"))
peaks<-peaks[peaks$chr %in% c(paste0("chr",1:22),"chrX"),]
peaks<-peaks[peaks$start>0,]
peaks<-makeGRangesFromDataFrame(peaks)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

make_seurat_object<-function(sample_input){
  wd=paste0(sample_input,"/","outs")
  
  #add a 0 in front of single digits for samples, just for sorting sake
  outname<-unlist(lapply(sample_input,function(x){
  if(nchar(strsplit(x,"_")[[1]][2])==1){
    paste(strsplit(x,"_")[[1]][1],paste0("0",strsplit(x,"_")[[1]][2]),sep="_")
  } else{x}}))

  #read in data
  counts <- Read10X_h5(paste0(wd,"/filtered_feature_bc_matrix.h5")) #count data
  fragpath <- paste0(wd,"/atac_fragments.tsv.gz") #atac fragments
  metadata_cellranger<-read.csv(paste0(wd,"/per_barcode_metrics.csv")) #metadata
  row.names(metadata_cellranger)<-metadata_cellranger$barcode
  soupx_output<-readRDS(paste0(wd,"/soupx_corrected_counts.rds")) #load SoupX contamination corrected output

  scrublet_output<-read.table(paste0(wd,"/","scrublet_results.tsv"),sep="\t",header=T) #load scrublet output for doublet detection
  scrublet_output$Barcode<-substr(scrublet_output$Barcode,2,nchar(scrublet_output$Barcode))
  row.names(scrublet_output)<-scrublet_output$Barcode

  # create a Seurat object containing the RNA data
  dat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA")

  # quantify counts in each peak (using merged peak set)
  #i actually used macs3
  fragments <- CreateFragmentObject(
    path = fragpath,
    cells = colnames(dat),
    validate.fragments = FALSE)

  atac_counts <- FeatureMatrix(
    fragments = fragments,
    features = peaks,
    cells = colnames(dat))

  # create ATAC assay and add it to the object
  dat[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation,
    min.features=-1)

  #Create corrected RNA data and add to object
  dat[["SoupXRNA"]]<-CreateAssayObject(counts=soupx_output)

  #QC cells
  DefaultAssay(dat) <- "ATAC"
  dat<-AddMetaData(dat,metadata=metadata_cellranger)
  dat<-AddMetaData(dat,metadata=scrublet_output)

  #Rename cells for sample specificity
  #set up basic filters
    dat <- NucleosomeSignal(dat)
    dat <- TSSEnrichment(dat)

    plt<-VlnPlot(
      object = dat,
      features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal","scrublet_Scores"),
      ncol = 4,
      pt.size = 0)

    ggsave(plt,file=paste0(outname,".qc.pdf"))
    dat$sample<-outname
    dat<-RenameCells(dat,new.names=paste(outname,Cells(dat),sep="_"))
  return(dat)
}

out<-lapply(sample_in,make_seurat_object)
dat <- merge(out[[1]], y = as.list(out[2:length(out)]), project = "all_data")

#metadata 
met<-read.table(opt$metadata,sep=",",header=T)
#add a 0 in front of single digits for samples, just for sorting sake
met$sample<-unlist(lapply(met$Manuscript_Name,function(x){
if(nchar(strsplit(x,"_")[[1]][2])==1){
  paste(strsplit(x,"_")[[1]][1],paste0("0",strsplit(x,"_")[[1]][2]),sep="_")
} else{x}}))

row.names(met)<-met$sample
meta_per_cell<-data.frame(met[dat$sample,],row.names=row.names(dat@meta.data))
dat<-AddMetaData(dat,met[dat$sample,])
saveRDS(dat,file="1_merged.unfiltered.SeuratObject.rds")

#Filter data
#First by read counts
dat<-subset(dat,subset = nCount_ATAC>=1000 & nCount_RNA>=1000)

#Match quantile scores to expected doublet rates by 10x
doublet_rate<-setNames(c(0.4,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8),
              nm=c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))

#Then by expected doublet rates
dat<-SplitObject(dat,split.by = "sample")
doublet_score_thresholding<-function(dat_sample){
  cell_count<-nrow(dat_sample@meta.data)
  expected_doublet_rate<-doublet_rate[which.min(as.numeric(names(doublet_rate))> cell_count)]/100
  doublet_score_threshold<-stats::quantile(dat_sample$scrublet_Scores,probs=1-expected_doublet_rate)
  dat_sample$scrublet_DropletType<-ifelse(dat_sample$scrublet_Scores<doublet_score_threshold,"singleton","doublet")
  return(dat_sample)
}
out<-lapply(dat,doublet_score_thresholding)
dat <- merge(out[[1]], y = as.list(out[2:length(out)]), project = "all_data")
dat<-subset(dat, subset = scrublet_DropletType=="singleton")
saveRDS(dat,file="2_merged.scrublet_count_filtered.SeuratObject.rds")
