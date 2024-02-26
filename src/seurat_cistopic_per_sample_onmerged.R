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
outdir=paste0(opt$output_directory,"/cistopic_objects") #/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/plots
system(paste0("mkdir -p ",outdir))

#filter to cells with ATAC features > 1000
#for samples with not enough ATAC data, skip
single_sample_cistopic_generation<-function(x,sample_in,outdir,epithelial_only=TRUE){
  if(epithelial_only){
      dat<-subset(dat,sample==sample_in)
      dat<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
      out_seurat_object<-paste0(sample_in,".cistopic_epithelial.SeuratObject.rds")
      model_selection_out<-paste0(outdir,"/",sample_in,".cistopic_epithelial.model_selection.pdf")
      out_cistopic_obj<-paste0(sample_in,".cistopic_epithelial.cistopicObject.rds")
      umap_out<-paste0(outdir,"/",sample_in,".cistopic_epithelial.umap.pdf")
    }
    else {
      dat<-subset(dat,sample==sample_in)
      out_seurat_object<-paste0(sample_in,".cistopic.SeuratObject.rds")
      model_selection_out<-paste0(outdir,"/",sample_in,".cistopic.model_selection.pdf")
      out_cistopic_obj<-paste0(sample_in,".cistopic.cistopicObject.rds")
      umap_out<-paste0(outdir,"/",sample_in,".cistopic_epithelial.umap.pdf")

    }

  #skip cistopic if cell count too low
  if(sum(dat$nCount_ATAC>1000)<500){
      saveRDS(dat,file=out_seurat_object)
  }else{
  dat<-subset(dat,nCount_ATAC>1000)

  cistopic_counts_frmt<-dat@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runModels(sub_cistopic,
    topic=seq(from=10, to=30, by=10),
    nCores=1,
    addModels=FALSE) #using v2 of cistopic (we are only using single core anyway)

  sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =dat@meta.data)
  sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
  
  print("finshed running cistopic")

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


single_sample_cistopic_generation(x=dat,outdir=outdir,sample_in=sample_in,epithelial_only=TRUE) #only epithelial cells per sample
single_sample_cistopic_generation(x=dat,outdir=outdir,sample_in=sample_in,epithelial_only=FALSE) #all cells per sample


