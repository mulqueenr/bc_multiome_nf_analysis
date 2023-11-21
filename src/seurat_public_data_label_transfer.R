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
seurat_obj_list=unlist(strsplit(args[1]," "))
outdir=args[2]
ref_dir=args[3] #/home/groups/CEDAR/mulqueen/bc_multiome/ref


# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #read in data
  outname=strsplit(x,"[.]")[[1]][1]
  dat<-readRDS(x)
  dat$sample<-outname #set up sample metadata
  print(paste("Finished sample:",outname))
  return(dat)}

out<-lapply(seurat_obj_list,merge_seurat)
sample_names=unlist(lapply(strsplit(seurat_obj_list,"[.]"),"[",1))
dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = sample_names, project = "all_data")


saveRDS(dat,file="merged.SeuratObject.rds")


single_sample_label_transfer<-function(dat,ref_obj,ref_prefix){
  out_plot<-paste0(outdir,"/",ref_prefix,".predictions.umap.pdf")
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
  plt1<-FeaturePlot(dat,features=colnames(predictions),pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by=paste0(ref_prefix,'_predicted.id'),pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  return(dat)
  }

swarbrick<-readRDS(paste0(ref_dir,"/swarbrick/swarbrick.SeuratObject.Rds"))#swarbrick types
dat<-single_sample_label_transfer(dat,swarbrick,"swarbrick")
rm(swarbrick)

embo_er<-readRDS(paste0(ref_dir,"/embo/SeuratObject_ERProcessed.rds"))#EMBO cell types
dat<-single_sample_label_transfer(dat,embo_er,"EMBO")
rm(embo_er)

hbca<-readRDS(paste0(ref_dir,"/hbca/hbca.rds")) #HBCA cell types
dat<-single_sample_label_transfer(dat,hbca,"HBCA")

saveRDS(dat,file="merged.SeuratObject.rds")

