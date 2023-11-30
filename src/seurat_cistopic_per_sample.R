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

args = commandArgs(trailingOnly=TRUE)
obj_in=args[1] #NAT_11.SeuratObject.rds
outdir=args[2] #/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/plots
obj_out=args[3] #/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects


#filter to cells with ATAC features > 1000
#for samples with not enough ATAC data, skip
single_sample_cistopic_generation<-function(x,outdir,obj_out){
  outname<-strsplit(x,"[.]")[[1]][1]
  atac_sub<-readRDS(x)
  #skip cistopic if cell count too low
  if(sum(atac_sub$nCount_ATAC>1000)<500){
      saveRDS(atac_sub,paste0(outname,".cistopic.SeuratObject.rds"))
  }else{
  atac_sub<-subset(atac_sub,nCount_ATAC>1000)

  cistopic_counts_frmt<-atac_sub@assays$peaks@counts
  row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
  sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
  print("made cistopic object")
  sub_cistopic_models<-cisTopic::runWarpLDAModels(sub_cistopic,
    topic=seq(from=10, to=30, by=5),
    nCores=1,
    addModels=FALSE)

  sub_cistopic_models<-addCellMetadata(sub_cistopic_models, cell.data =atac_sub@meta.data)
  pdf(paste0(outdir,"/",outname,"_model_selection.pdf"))
  par(mfrow=c(3,3))
  sub_cistopic_models<- selectModel(sub_cistopic_models, type='derivative')
  dev.off()
  
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
  atac_sub@reductions$cistopic<-cistopic_obj
  n_topics<-ncol(Embeddings(atac_sub,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
  print("Running UMAP")
  atac_sub<-RunUMAP(atac_sub,reduction="cistopic",dims=1:n_topics)
  atac_sub <- FindNeighbors(object = atac_sub, reduction = 'cistopic', dims = 1:n_topics ) 
  atac_sub <- FindClusters(object = atac_sub, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
  print("Plotting UMAPs")
  plt1<-DimPlot(atac_sub,reduction="umap",group.by=c("seurat_clusters"))
  #plt2<-FeaturePlot(atac_sub,reduction="umap",features=c("nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks"))
  pdf(paste0(outdir,"/",outname,".cistopic.umap.pdf"),width=10)
  print(plt1)
  #print(plt2)
  dev.off()
  saveRDS(sub_cistopic_models,file=paste0(outname,".CisTopicObject.rds"))
  saveRDS(atac_sub,paste0(outname,".cistopic.SeuratObject.rds"))
  }
  }


single_sample_cistopic_generation(x=obj_in,outdir=outdir,obj_out=obj_out)


