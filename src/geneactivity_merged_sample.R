library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(SeuratObjects)
library(EnsDb.Hsapiens.v86)
args = commandArgs(trailingOnly=TRUE)

dat=readRDS(args[1])
outname<-strsplit(args[1],"[.]")[1]

gene_activity<-GeneActivity(dat,process_n=10000)
saveRDS(gene_activity,file=paste0(outname,".GeneActivity.rds"))

dat[["GeneActivity"]]<-CreateAssayObject(counts=gene_activity)
dat<- NormalizeData(
  object = dat,
  assay = "GeneActivity",
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_GeneActivity)
)
saveRDS(dat,file=paste0(outname,".geneactivity.SeuratObject.rds"))


