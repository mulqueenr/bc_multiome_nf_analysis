library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(plyr)

args = commandArgs(trailingOnly=TRUE)
seurat_obj_list=args[1]
ref_dir=args[2] #/home/groups/CEDAR/mulqueen/bc_multiome/ref

#add metadata.tsv per sample if supplied
if(length(args>2)){
  met<-read.csv(args[3],header=T,sep=",") #/home/groups/CEDAR/mulqueen/bc_multiome/sample_metadata.csv

}
sample_metadata=args[3]
# set up sample loop to load the RNA and ATAC data, save to seurat object
merge_seurat<-function(x){
  #read in data
  outname=strsplit(x,"[.]")[[1]][1]
  dat<-readRDS(x)
  dat$sample<-outname #set up sample metadata
  print(paste("Finished sample:",outname))
  return(dat)}

seurat_obj_list=strsplit(seurat_obj_list," ")[[1]]
out<-lapply(unlist(seurat_obj_list),merge_seurat)
sample_names=unlist(lapply(strsplit(seurat_obj_list,"[.]"),"[",1))
dat <- merge(out[[1]], y = as.list(out[2:length(out)]), add.cell.ids = sample_names, project = "all_data")
dat_met<-as.data.frame(dat@meta.data)
met<-merge(dat_met,met,by.x="sample",by.y="Manuscript_Name")
row.names(met)<-row.names(dat_met)
dat<-AddMetaData(dat,met)

saveRDS(dat,file="merged.SeuratObject.rds")
#dat<-readRDS(file="merged.SeuratObject.rds"))

#prepare data
DefaultAssay(dat)<-"RNA" #can use SoupXRNA here also
dat<-NormalizeData(dat)
dat<-FindVariableFeatures(dat)
dat<-ScaleData(dat)
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

single_sample_label_transfer<-function(dat,ref_obj,ref_prefix){

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
  
  return(dat)
  }

swarbrick<-readRDS(paste0(ref_dir,"/swarbrick/swarbrick.SeuratObject.Rds"))#swarbrick types
dat<-single_sample_label_transfer(dat,ref_obj=swarbrick,"swarbrick")
rm(swarbrick)

embo_er<-readRDS(paste0(ref_dir,"/embo/SeuratObject_ERProcessed.rds"))#EMBO cell types
dat<-single_sample_label_transfer(dat,ref_obj=embo_er,"EMBO")
rm(embo_er)

hbca<-readRDS(paste0(ref_dir,"/hbca/hbca.rds")) #HBCA cell types
dat<-single_sample_label_transfer(dat,ref_obj=hbca,"HBCA")

saveRDS(dat,file="merged.public_transfer.SeuratObject.rds")

