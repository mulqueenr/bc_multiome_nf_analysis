library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(plyr)
library(optparse)

option_list = list(
  make_option(c("-s", "--sample_list"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character"),
    make_option(c("-r", "--ref_dir"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref", 
              help="Reference directory containing genome information. default: %default]", metavar="character"),
    make_option(c("-m","--metadata"), type="character", default=NULL,
              help="Comma separated (CSV) metadata file of cell information to be appended.", metavar="character"),
    make_option(c("-o","--plot_output_directory"), type="character", default=NULL,
              help="Directory to publish output plots to.", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
seurat_obj_list=opt$sample_list
ref_dir=opt$ref_dir
met<-read.csv(opt$metadata,header=T,sep=",")
outdir<-opt$plot_output_directory

#######testing#########
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
#seurat_obj_list="IDC_12.SeuratObject.rds IDC_4.SeuratObject.rds DCIS_2.SeuratObject.rds IDC_7.SeuratObject.rds IDC_2.SeuratObject.rds IDC_5.SeuratObject.rds IDC_9.SeuratObject.rds IDC_6.SeuratObject.rds NAT_4.SeuratObject.rds NAT_11.SeuratObject.rds IDC_3.SeuratObject.rds ILC_1.SeuratObject.rds IDC_10.SeuratObject.rds IDC_1.SeuratObject.rds DCIS_1.SeuratObject.rds IDC_11.SeuratObject.rds DCIS_3.SeuratObject.rds NAT_14.SeuratObject.rds IDC_8.SeuratObject.rds"
#ref_dir="/home/groups/CEDAR/mulqueen/bc_multiome/ref" 
#met=read.csv("sample_metadata.csv",header=T,sep=",")
#outdir<-"/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/plots"
###################


#saveRDS(dat,file="merged.SeuratObject.rds")
#dat<-readRDS(file="merged.SeuratObject.rds"))

#prepare data
DefaultAssay(dat)<-"RNA" #can use SoupXRNA here also
#dat<-NormalizeData(dat)
#dat<-FindVariableFeatures(dat)
#dat<-ScaleData(dat)
dat <- SCTransform(dat,vst.flavor = 'v2')
dat <- RunPCA(dat)
dat<- RunUMAP(
  object = dat,
  reduction.name="rna_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

single_sample_label_transfer<-function(dat,ref_obj,ref_prefix,celltype="celltype"){
  ref_obj<-UpdateSeuratObject(ref_obj)
  ref_obj <- SCTransform(ref_obj,vst.flavor = 'v2')

  transfer.anchors <- FindTransferAnchors(
    reference = ref_obj,
    reference.assay="SCT",
    query = dat,
    query.assay="SCT",
    features=VariableFeatures(ref_obj),
    verbose=T,
    dims=1:30
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = ref_obj[,celltype],
  )
  colnames(predictions)<-paste0(ref_prefix,"_",colnames(predictions))

  dat<-AddMetaData(dat,metadata=predictions)
  
  return(dat)
  }

swarbrick<-readRDS(paste0(ref_dir,"/swarbrick/swarbrick.SeuratObject.Rds"))#swarbrick types
dat<-single_sample_label_transfer(dat,ref_obj=swarbrick,ref_prefix="swarbrick")
rm(swarbrick)

embo_er<-readRDS(paste0(ref_dir,"/embo/SeuratObject_ERProcessed.rds")) #EMBO cell types
dat<-single_sample_label_transfer(dat,ref_obj=embo_er,ref_prefix="EMBO")
rm(embo_er)

hbca<-readRDS(paste0(ref_dir,"/hbca/hbca.rds")) #HBCA cell types
dat<-single_sample_label_transfer(dat,ref_obj=hbca,ref_prefix="HBCA")


nakshatri<-readRDS(file=paste0(ref_dir,"/nakshatri/","nakshatri_multiome.rds"))
dat<-single_sample_label_transfer(dat,celltype="author_cell_type",ref_obj=nakshatriref_prefix="nakshatri")

saveRDS(dat,file="merged.public_transfer.SeuratObject.rds")

