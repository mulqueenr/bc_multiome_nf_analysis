library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
set.seed(1234)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
library(rtracklayer)
library(ggplot2)
library(optparse)
#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
option_list = list(
  make_option(c("-s", "--sample_array_in"), type="character", default="NAT_1", 
              help="Sample array to subset to (number 1 - total sample count", metavar="character"),
    make_option(c("-i", "--object_input"), type="character", default="merged.public_transfer.SeuratObject.rds", 
              help="Reference directory containing genome information. default: %default]", metavar="character"),
        make_option(c("-o", "--output_directory"), type="character", default=NULL, 
              help="Output directory, defined in nextflow parameters.", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.passqc.SeuratObject.rds"
dat<-readRDS(file=opt$object_input)
opt$output_directory="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3"
outdir=paste0(opt$output_directory,"/cistopic_objects") #/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
system(paste0("mkdir -p ",outdir))
opt$sample_array_in<-1
sample_in=unique(dat$sample)[opt$sample_array_in]

#filter to cells with ATAC features > 1000
#for samples with not enough ATAC data, skip
single_sample_cistopic_generation<-function(Object,sample_in,outdir,epithelial_only=TRUE){
  print(paste0("Running cistopic on ",sample_in," ..."))
  if(epithelial_only){
      dat<-subset(Object,sample==sample_in)
      dat<-subset(dat,reclust %in% c("luminal_epithelial","basal_epithelial"))
      out_seurat_object<-paste0(outdir,"/",sample_in,".cistopic_epithelial.SeuratObject.rds")
      model_selection_out<-paste0(outdir,"/",sample_in,".cistopic_epithelial.model_selection.pdf")
      out_cistopic_obj<-paste0(outdir,"/",sample_in,".cistopic_epithelial.cistopicObject.rds")
      umap_out<-paste0(outdir,"/",sample_in,".cistopic_epithelial.umap.pdf")
    }
    else {
      dat<-subset(Object,sample==sample_in)
      out_seurat_object<-paste0(outdir,"/",sample_in,".cistopic.SeuratObject.rds")
      model_selection_out<-paste0(outdir,"/",sample_in,".cistopic.model_selection.pdf")
      out_cistopic_obj<-paste0(outdir,"/",sample_in,".cistopic.cistopicObject.rds")
      umap_out<-paste0(outdir,"/",sample_in,".cistopic_epithelial.umap.pdf")

    }

  #skip cistopic if cell count too low
  if(sum(dat$nCount_ATAC>1000)<500){
      print("Cell count for ATAC seq is too low...")
      saveRDS(dat,file=out_seurat_object)
  }else{
  dat<-subset(dat,nCount_ATAC>1000)

  cistopic_counts_frmt<-dat@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("Made cistopic object")
  sub_cistopic_models<-cisTopic::runModels(sub_cistopic,
    topic=seq(from=10, to=30, by=5),
    nCores=1,
    addModels=FALSE) #using v2 of cistopic (we are only using single core anyway)

  sub_cistopic_models<-cisTopic::addCellMetadata(sub_cistopic_models, cell.data =dat@meta.data)
  sub_cistopic_models<- cisTopic::selectModel(sub_cistopic_models, type='derivative')
  
  print("Finshed running cistopic")

  #Add cell embeddings into seurat
  cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
  colnames(cell_embeddings)<-sub_cistopic_models@cell.names
  n_topics<-nrow(cell_embeddings)
  row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
  cell_embeddings<-as.data.frame(t(cell_embeddings))

  #Add feature loadings into seurat
  feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
  row.names(feature_loadings)<-paste0("topic_",1:n_topics)
  feature_loadings<-as.data.frame(t(feature_loadings))

  #combined cistopic results (cistopic loadings and umap with seurat object)
  cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
  print("Cistopic Loading into Seurat")
  dat@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(dat,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  dat<-RunUMAP(dat,reduction="cistopic",dims=1:n_topics)
  dat <- FindNeighbors(object = dat, reduction = 'cistopic', dims = 1:n_topics ) 
  dat <- FindClusters(object = dat, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(dat,reduction="umap",group.by=c("HBCA_predicted.id"))
  pdf(umap_out,width=10)
  print(plt1)
  dev.off()
  saveRDS(sub_cistopic_models,file=out_cistopic_obj)
  saveRDS(dat,file=out_seurat_object)
  }
}


lapply(1:length(unique(dat$sample)), function(x) {
sample_in=unique(dat$sample)[x]
single_sample_cistopic_generation(Object=dat,outdir=outdir,sample_in=sample_in,epithelial_only=TRUE) #only epithelial cells per sample
single_sample_cistopic_generation(Object=dat,outdir=outdir,sample_in=sample_in,epithelial_only=FALSE) #all cells per sample
})

#IDC_09

