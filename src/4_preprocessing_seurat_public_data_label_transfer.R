#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome,/home/groups/CEDAR/mulqueen/bc_multiome/ref:/ref $sif

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(plyr)
library(optparse)
library(ggplot2)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="2_merged.scrublet_count_filtered.SeuratObject.rds", 
              help="Input seurat object", metavar="character"),
    make_option(c("-r", "--ref_dir"), type="character", default="/ref", 
              help="Reference directory containing genome information. default: %default]", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
ref_dir=opt$ref_dir
dat<-readRDS(file=opt$object_input)
dat <- SCTransform(dat,vst.flavor='v2')
DefaultAssay(dat)<-"SCT"

#dat<-JoinLayers(dat)
#using only SCT for all label transfers (not atac for the multimodal ones)

single_sample_label_transfer<-function(dat,ref_obj,ref_prefix,celltype="celltype"){
  ref_obj<-UpdateSeuratObject(ref_obj)
  DefaultAssay(ref_obj)<-"RNA"
  ref_obj <- SCTransform( ref_obj,vst.flavor='v2')
  ref_obj <- RunPCA(ref_obj)

  transfer.anchors <- FindTransferAnchors(
    reference = ref_obj,
    normalization.method = "SCT",
    query = dat,
    features=VariableFeatures(ref_obj),
    verbose=T,
    dims=1:30,
    reference.reduction = "pca"
  )

  predictions <- TransferData(
    anchorset = transfer.anchors,
    refdata = ref_obj@meta.data[,celltype],
  )
  colnames(predictions)<-paste0(ref_prefix,"_",colnames(predictions))

  dat<-AddMetaData(dat,metadata=predictions)

  #plot umap of reference RNA and projection of our data in theirs
  ref_obj <- RunUMAP(ref_obj, dims = 1:30, reduction = "pca", return.model = TRUE)

  dat <- MapQuery(anchorset = transfer.anchors, 
    reference = ref_obj, 
    query = dat,
    refdata = list(celltype_ref = celltype), 
    reference.reduction = "pca", 
    reduction.model = "umap")
    
  p1 <- DimPlot(ref_obj, reduction = "umap", group.by = celltype, 
    label = TRUE, label.size = 3, repel = TRUE,raster=FALSE) + 
    NoLegend() + 
    ggtitle(paste(ref_prefix,"Reference annotations"))

  p2 <- DimPlot(dat, reduction = "ref.umap", group.by = 'predicted.celltype_ref', 
    label = TRUE, label.size = 3, repel = TRUE,raster=FALSE) + 
    NoLegend() + ggtitle(paste(ref_prefix,"Query transferred labels"))
  
  
  ggsave(p1 + p2,file=paste0(ref_prefix,"_label_transfer.rna.umap.pdf"),width=10)
  return(dat)
  }

swarbrick_file<-paste0(ref_dir,"/swarbrick/swarbrick.SeuratObject.Rds")
if(file.exists(swarbrick_file)){
swarbrick<-readRDS(swarbrick_file)#swarbrick types
dat<-single_sample_label_transfer(dat,ref_obj=swarbrick,ref_prefix="wu")
rm(swarbrick)
} else {print("Wu et al. not found.")}

embo_er_file<-paste0(ref_dir,"/embo/SeuratObject_ERProcessed.rds")
if(file.exists(embo_er_file)){
embo_er<-readRDS(embo_er_file) #EMBO cell types
dat<-single_sample_label_transfer(dat,ref_obj=embo_er,ref_prefix="pal")
rm(embo_er)
} else {print("Pal et al. not found.")}

hbca_file<-paste0(ref_dir,"/hbca/hbca.rds")
if(file.exists(hbca_file)){
hbca<-readRDS(hbca_file) #HBCA cell types
dat<-single_sample_label_transfer(dat,ref_obj=hbca,ref_prefix="kumar")
rm(hbca)
} else {print("Kumar et al. not found.")}

nakshatri_file<-paste0(paste0(ref_dir,"/nakshatri/","nakshatri_multiome.geneactivity.rds"))
if(file.exists(nakshatri_file)){
nakshatri<-readRDS(nakshatri_file)
dat<-single_sample_label_transfer(dat,ref_obj=nakshatri,ref_prefix="nakshatri")
rm(nakshatri)
} else {print("Bhat-Nakshatri et al. not found.")}

terekhanova_file<-paste0(paste0(ref_dir,"/terekhanova/","terekhanova_multiome.geneactivity.rds"))
if(file.exists(terekhanova_file)){
terekhanova<-readRDS(terekhanova_file)
dat<-single_sample_label_transfer(dat,ref_obj=terekhanova,ref_prefix="terekhanova",celltype="cell_type")
rm(terekhanova)
} else {print("Terekhanova et al. not found.")}

saveRDS(dat,file="3_merged.public_transfer.SeuratObject.rds")

