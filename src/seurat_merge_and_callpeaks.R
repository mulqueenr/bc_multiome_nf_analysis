library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(RColorBrewer)

#args[1]=list of seurat files
wd=args[1]

setwd(paste0(wd,"/","outs"))

# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #function to handle different sample directories##################
  if(x %in% 1:12){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else if(x %in% 13:20){
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
  outname<-paste0("sample_",x)
  }else{
  wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
  outname<-x
  }
  ####################################################################
  #read in data
  dat<-readRDS(paste0(wd,"/",outname,".SeuratObject.rds"))
  dat$sample<-outname #set up sample metadata
  return(dat)}

out<-lapply(c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20,"RM_1","RM_2","RM_3","RM_4"),merge_seurat)


dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = c(paste0("sample_",c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,19,20)),"RM_1","RM_2","RM_3","RM_4"), project = "all_data")
saveRDS(dat,file="/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/phase2.SeuratObject.rds")

dat<-readRDS("phase2.SeuratObject.rds")
dat
table(dat$sample)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

# call peaks using MACS2
DefaultAssay(dat)<-"ATAC"
peaks <- CallPeaks(dat, macs2.path = "/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2")
#use this set of peaks for all samples

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

DefaultAssay(dat) <- "ATAC"
saveRDS(peaks,file="combined.peakset.rds")

# quantify counts in each peak
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

saveRDS(dat,file="phase2.SeuratObject.rds")