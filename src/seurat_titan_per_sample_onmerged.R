library(TITAN)
library(Seurat)
library(tidyverse)
library(Signac)
set.seed(1234)
library(optparse)

option_list = list(
  make_option(c("-s", "--sample"), type="character", default="NAT_1", 
              help="Sample name to subset to", metavar="character"),
    make_option(c("-i", "--input_merged_object"), type="character", default="merged.public_transfer.SeuratObject.rds", 
              help="Reference directory containing genome information. default: %default]", metavar="character"),
        make_option(c("-o", "--output_directory"), type="character", default=NULL, 
              help="Output directory, defined in nextflow parameters.", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample_in=strsplit(opt$sample,"[.]")[[1]][1]
dat<-readRDS(file=opt$input_merged_object)
outdir=paste0(opt$output_directory,"/titan_objects") #/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/plots
system(paste0("mkdir -p ",outdir))

add_topics_to_seurat_obj <- function(model, Object) {
    modelMat <- t(scale(model$document_expects, center = TRUE, scale = TRUE))
    rownames(modelMat) <- paste(1:ncol(Object), colnames(Object), sep = "_")
    colnames(modelMat) <- paste("Topic", 1:ncol(modelMat), sep = "_")
    Object@meta.data <- cbind(Object@meta.data, modelMat)
    rownames(modelMat) = colnames(Object)
    Object[["lda"]] <- CreateDimReducObject(embeddings = modelMat, 
        key = "lda_", assay = "RNA", global = TRUE)
    return(Object)
}

#filter to cells with RNA reads > 500
single_sample_titan_generation<-function(x,outdir,epithelial_only=TRUE,sample_in){
    if(epithelial_only){
      dat<-subset(x,sample==sample_in)
      dat<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
      out_seurat_object<-paste0(sample_in,".titan_epithelial.SeuratObject.rds")
      out_titan_obj<-paste0(sample_in,".titan_epithelial.titanObject.rds")
      umap_out<-paste0(outdir,"/",sample_in,".titan_epithelial.umap.pdf")
    }
    else {
      dat<-subset(x,sample==sample_in)
      out_seurat_object<-paste0(sample_in,".titan.SeuratObject.rds")
      out_titan_obj<-paste0(sample_in,".titan.cistopicObject.rds")
      umap_out<-paste0(outdir,"/",sample_in,".titan_epithelial.umap.pdf")

    }

  #skip titan if cell count too low
  if(sum(dat$nCount_RNA>500)<500){
      saveRDS(dat,file=out_seurat_object)
  }else{
  dat<-subset(dat,nCount_RNA>500)
  DefaultAssay(dat)<-"RNA"
  LDA_model <- runLDA(dat, ntopics = 20, normalizationMethod = "CLR", outDir=paste0(outdir,"/",sample_in,"/"))
  #using 20 topics for each sample
  
  #setting top model to 20 topics for all (this can be changed for a range of topics)
  top_model <- LDA_model
  dat <- addTopicsToSeuratObject(model = LDA_model, Object = dat)
  TopicGenes <- TopTopicGenes(top_model, ngenes = 50)
  LDA_topics <- GetTopics(top_model, dat)

  print("Running UMAP")
  dat<-RunUMAP(dat,reduction="lda",dims=1:20)
  #dat <- FindNeighbors(object = dat, reduction = 'lda', dims = 1:20 ) 
  #dat <- FindClusters(object = dat, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(dat,reduction="umap",group.by=c("HBCA_predicted.id"))
  pdf(umap_out,width=10)
  print(plt1)
  dev.off()

  saveRDS(LDA_model, file=out_titan_obj)
  saveRDS(dat,file=out_seurat_object)
  }
}


single_sample_titan_generation(x=obj_in,outdir=outdir,sample_in=sample_in,epithelial_only=TRUE)
single_sample_titan_generation(x=obj_in,outdir=outdir,sample_in=sample_in,epithelial_only=FALSE) #all cells per sample


