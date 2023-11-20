library(TITAN)
library(Seurat)
library(tidyverse)
library(Signac)
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)
obj_in=args[1] #IDC_7.SeuratObject.rds
outdir=args[2]
obj_out=args[3]

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
single_sample_titan_generation<-function(x,outdir,obj_out=obj_out){
  outname<-strsplit(x,"[.]")[[1]][1]
  dat<-readRDS(x)
  dat<-subset(dat,nCount_RNA>500)
  DefaultAssay(dat)<-"RNA"
  LDA_model <- runLDA(dat, ntopics = 20, normalizationMethod = "CLR", outDir=paste0(obj_out,"/",outname,"/"))
  #using 20 topics for each sample
  
  #setting top model to 20 topics for all (this can be changed for a range of topics)
  top_model <- LDA_model
  dat <- addTopicsToSeuratObject(model = LDA_model, Object = dat)
  TopicGenes <- TopTopicGenes(top_model, ngenes = 50)
  LDA_topics <- GetTopics(top_model, dat)

  #pdf(paste0(outdir,"/",outname,".titan.elbowplot.pdf"),width=10)
  #HeatmapTopic(Object = dat,
  #        topics =  LDA_topics,
  #        AnnoVector = dat@meta.data$HBCA_predicted.id,
  #        AnnoName = "HBCA_predicted.id")
  #dev.off()
  saveRDS(LDA_model, paste0(obj_out,"/",outname,"TITANObject.rds"))
  saveRDS(dat,paste0(outname,".SeuratObject.rds"))
  }


single_sample_titan_generation(x=obj_in,outdir=outdir,obj_out=obj_out)


