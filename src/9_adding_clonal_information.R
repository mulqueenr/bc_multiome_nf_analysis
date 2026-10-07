
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(dplyr)
library(optparse)
library(parallel)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(grid)
library(ggrepel)
library(seriation)
library(org.Hs.eg.db)
library(dendextend)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GeneNMF)
library(rGREAT)
library(msigdbr,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(fgsea,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(presto,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
setwd("/home/groups/MohammedLab/bc_multiome/seurat_object")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

############# Add clone and filtering info to object
read_clone_annot<-function(i){
clone_annot<-list.files(i,pattern="_ATAC_anno.csv",full.names=T)
  if(length(clone_annot)!= 0 ){
    sample_name<-gsub(basename(clone_annot),pattern="_ATAC_anno.csv",replace="")
    clones<-read.csv(clone_annot,row.names=1)
    row.names(clones)<-gsub(row.names(clones),pattern="[.]",replace="-")
    row.names(clones)<-paste(sample_name,row.names(clones),sep="_")
    clones$merge_cluster<-paste(sample_name,clones$merge_cluster,sep="_")
    clones$atac_cluster<-paste(sample_name,clones$atac_cluster,sep="_")
    clones$rna_cluster<-paste(sample_name,clones$rna_cluster,sep="_")
    return(clones)
  }
}

############# Add clone and filtering info to object
read_cnv_plots<-function(i){
  sample<-gsub(basename(i),pattern="_consensus_cnvs.bed",replacement="")
  cnv<-read.table(i)
  cnv<-cnv[5:ncol(cnv)]
  colnames(cnv)<-gsub(colnames(cnv),pattern="[.]",replacement="-")
  colnames(cnv)<-paste(sample,colnames(cnv),sep="_")
  return(cnv)
}

clone_annot_list<-list.files("/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/",pattern="*clone_val",full.name=T)
clone_annot<-do.call("rbind",lapply(clone_annot_list,read_clone_annot))
colnames(clone_annot)<-c("merged_assay_clones","atac_clones","rna_clones","clones_celltype","clones_log2reads")
dat<-AddMetaData(dat,clone_annot)

cnv_list<-list.files("/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/",recursive=T,pattern="*_consensus_cnvs.bed",full.name=T)
cnv_list<-cnv_list[!grepl(cnv_list,pattern="unfiltered")]
cnv_windows<-read.table(cnv_list[1])
cnv_windows<-cnv_windows[1:4]
colnames(cnv_windows)<-c("chr","start","end","win")

cnv_out<-do.call("cbind",lapply(cnv_list,read_cnv_plots))
row.names(cnv_out)<-paste(cnv_windows$chr,cnv_windows$start,cnv_windows$end,sep="-")
cnv_out<-cnv_out[colnames(cnv_out) %in% row.names(dat@meta.data)]
cnv_assay <- CreateAssayObject(counts = cnv_out)
dat[["cnv"]]<-cnv_assay

#normal clones based on lack of CNV calls and co-clustering
normal_clones<-c(
"DCIS_03_1",
"IDC_01_4",
"IDC_02_1",
"IDC_05_1",
"IDC_06_3",
"IDC_08_2",
"IDC_09_2",
"IDC_10_3",
"IDC_11_2",
"IDC_12_4",
"ILC_02_1",
"ILC_04_2",
"NAT_14_1",
"IDC_15_1")

contamination_clones<-c(
"ILC_04_3",
"ILC_04_4"
)

dat@meta.data[dat@meta.data$merged_assay_clones %in% normal_clones,]$merged_assay_clones<-"normal"

#reassign any epithelial cells not assigned to "normal" group as cancer
dat@meta.data[!is.na(dat@meta.data$merged_assay_clones) & dat@meta.data$merged_assay_clones != "normal",]$assigned_celltype<-"cancer"
dat@meta.data[dat@meta.data$merged_assay_clones %in% contamination_clones,]$merged_assay_clones<-"contamination"
dat@meta.data[dat@meta.data$merged_assay_clones %in% c("normal"),]$merged_assay_clones <-paste(dat@meta.data[dat@meta.data$merged_assay_clones %in% c("normal"),]$sample,"normal",sep="_")

saveRDS(dat,file="8_merged.cnv_clones.SeuratObject.rds")

