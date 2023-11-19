library(TITAN)
library(Seurat)
library(tidyverse)
library(Signac)

args = commandArgs(trailingOnly=TRUE)
obj_in=args[1]
outdir=args[2]

single_sample_cistopic_generation<-function(x,outdir){
  outname<-strsplit(x,"[.]")[1]
  dat<-readRDS(x)
  DefaultAssay(dat)<-"RNA"
  LDA_model <- runLDA(dat, ntopics = 20, normalizationMethod = "CLR", seed.number = 8)
  #using 20 topics for each sample
  
  #setting top model to 20 topics for all (this can be changed for a range of topics)
  top_model <- LDA_model
  dat <- addTopicsToSeuratObject(model = top_model, Object = dat)
  TopicGenes <- TopTopicGenes(top_model, ngenes = 50)
  LDA_topics <- GetTopics(top_model, dat)

  pdf(paste0(outdir,"/",outname,".titan.elbowplot.pdf"),width=10)
  HeatmapTopic(Object = dat,
          topics =  LDA_topics,
          AnnoVector = dat@meta.data$HBCA_predicted.id,
          AnnoName = "HBCA_predicted.id")

  saveRDS(LDA_model, paste0(outname,"TITANObject.rds"))
  saveRDS(dat,paste0(outname,".SeuratObject.rds"))
  }


single_sample_titan_generation(x=obj_in,outdir=outdir)


