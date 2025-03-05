#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

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
  make_option(c("-i", "--object_input"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects/merged.geneactivity.SeuratObject.rds", 
              help="Input seurat object", metavar="character"),
    make_option(c("-r", "--ref_dir"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref", 
              help="Reference directory containing genome information. default: %default]", metavar="character"),
    make_option(c("-m","--metadata"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/sample_metadata.csv",
              help="Comma separated (CSV) metadata file of cell information to be appended.", metavar="character"),
    make_option(c("-o","--plot_output_directory"), type="character", default=NULL,
              help="Directory to publish output plots to.", metavar="character")

); 
 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
ref_dir=opt$ref_dir
met<-read.csv(opt$metadata,header=T,sep=",")
outdir<-opt$plot_output_directory
dat<-readRDS(file=opt$object_input)

#prepare data, use SCT to normalize
DefaultAssay(dat)<-"RNA" 
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

swarbrick_file<-paste0(ref_dir,"/swarbrick/swarbrick.SeuratObject.Rds")
if(file.exists(swarbrick_file)){
swarbrick<-readRDS(swarbrick_file)#swarbrick types
dat<-single_sample_label_transfer(dat,ref_obj=swarbrick,ref_prefix="swarbrick2")
rm(swarbrick)
}

embo_er_file<-paste0(ref_dir,"/embo/SeuratObject_ERProcessed.rds")
if(file.exists(embo_er_file)){
embo_er<-readRDS(embo_er_file) #EMBO cell types
dat<-single_sample_label_transfer(dat,ref_obj=embo_er,ref_prefix="wu2")
rm(embo_er)
}

hbca_file<-paste0(ref_dir,"/hbca/hbca.rds")
if(file.exists(hbca_file)){
hbca<-readRDS(hbca_file) #HBCA cell types
dat<-single_sample_label_transfer(dat,ref_obj=hbca,ref_prefix="HBCA2")
rm(hbca)
}

nakshatri_file<-paste0(paste0(ref_dir,"/nakshatri/","nakshatri_multiome.rds"))
if(file.exists(nakshatri_file)){
nakshatri<-readRDS(nakshatri_file)
dat<-single_sample_label_transfer(dat,ref_obj=nakshatri,ref_prefix="nakshatri")
rm(nakshatri)
}

saveRDS(dat,file="merged.public_transfer.SeuratObject.rds")

