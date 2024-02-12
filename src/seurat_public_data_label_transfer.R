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
dat$cellID<-row.names(dat@meta.data)
dat_met<-as.data.frame(dat@meta.data)
met_merged<-merge(dat_met,met,by.x="sample",by.y="Manuscript_Name",all.x=T)
row.names(met_merged)<-met_merged$cellID
dat<-AddMetaData(dat,met_merged[,c("phase","Original_Sample","sample_ID","sample_weight","Diagnosis","Mol_Diagnosis","sampled_site","batch","outcome")])

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
  assay = "RNA",
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

