library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(plyr)
library(patchwork)

args = commandArgs(trailingOnly=TRUE)
obj_in=args[1]
outdir=args[2]
ref_dir=args[3] #/home/groups/CEDAR/mulqueen/bc_multiome/ref


single_sample_label_transfer<-function(x,ref_obj,ref_prefix){
  outname=strsplit(x,"[.]")[[1]][1]
  out_plot<-paste0(outdir,"/",outname,".",ref_prefix,".predictions.umap.pdf")
  dat<-readRDS(x)
  DefaultAssay(dat)<-"RNA" #can use SoupXRNA here also
  dat<-NormalizeData(dat)
  dat<-FindVariableFeatures(dat)
  dat<-ScaleData(dat)

  transfer.anchors <- FindTransferAnchors(
    reference = ref_obj,
    reference.assay="RNA",
    query = dat,
    query.assay="RNA",
    features=VariableFeatures(ref_obj),
    verbose=T
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = ref_obj$celltype,
  )
  colnames(predictions)<-paste0(ref_prefix,"_",colnames(predictions))

  dat<-AddMetaData(dat,metadata=predictions)
  saveRDS(dat,file=x)
  plt1<-FeaturePlot(dat,features=colnames(predictions),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by=paste0(ref_prefix,'_predicted.id'),pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  }

swarbrick<-readRDS(paste0(ref_dir,"/swarbrick/swarbrick.SeuratObject.Rds"))#swarbrick types
single_sample_label_transfer(obj_in,swarbrick,"swarbrick")
rm(swarbrick)

embo_er<-readRDS(paste0(ref_dir,"/embo/SeuratObject_ERProcessed.rds"))#EMBO cell types
single_sample_label_transfer(obj_in,embo_er,"EMBO")
rm(embo_er)

hbca<-readRDS(paste0(ref_dir,"/hbca/hbca.rds")) #HBCA cell types
single_sample_label_transfer(obj_in,hbca,"HBCA")


