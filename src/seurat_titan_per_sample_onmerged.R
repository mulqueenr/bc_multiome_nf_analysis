library(TITAN)
library(Seurat)
library(parallel)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(Signac)
set.seed(1234)
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
#setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
#opt$object_input="merged.geneactivity.SeuratObject.rds"
dat<-readRDS(file=opt$object_input)
#opt$output_directory="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3"
outdir=paste0(opt$output_directory,"/titan_objects") #/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
system(paste0("mkdir -p ",outdir))
#opt$sample_array_in<-1
sample_in=unique(dat$sample)[opt$sample_array_in]


model_maker <- function(topics,outDir,cellList,Genes,iterations,burnin,alpha,beta) {
    selected.Model <- lda::lda.collapsed.gibbs.sampler(
      cellList, topics, Genes, 
      num.iterations = iterations, 
      alpha = alpha, eta = beta, 
      compute.log.likelihood = TRUE, 
      burnin = burnin)[-1]
      saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
    }

RPC_calculation<-function(model_file,outDir,data.use) {
    topic_num <- as.numeric(gsub("[^0-9]+([0-9]+).*", "\\1", model_file))
    topic_numbers <- c(topic_numbers, topic_num)
    model <- readRDS(paste0(outDir, "/", model_file))
    docterMat <- t(as.matrix(data.use))
    docterMat <- as(docterMat, "sparseMatrix")
    topworddist <- normalize(model$topics, byrow = T)
    doctopdist <- normalize(t(model$document_sums), byrow = T)
    perp <- text2vec::perplexity(docterMat, topworddist, doctopdist)
    return(c(topic_num,perp))
}
   
single_sample_titan_generation<- function(Object,
  assay="SCT",
  seed.number=123,
  iterations=500,
  burnin=250,
  alpha_val=50,
  beta_val=0.1,
  nfeat=10000,
  outDir="./TITAN_LDA",
  topic_counts=seq(from=10, to=30, by=5),
  epithelial_only=FALSE,
  sample_in){
  
  print(paste0("Running TITAN on ",sample_in," ..."))
  set.seed(seed.number)
  print(paste0("Setting ",assay," as assay..."))
  if(epithelial_only){
    Object<-subset(Object,sample==sample_in)
    Object<-subset(Object,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
    out_seurat_object<-paste0(outDir,"/",sample_in,".titan_epithelial.SeuratObject.rds")
    out_titan_obj<-paste0(outDir,"/",sample_in,".titan_epithelial.titanObject.rds")
    elbow_out<-paste0(outDir,"/",sample_in,".titan_epithelial.elbow.pdf")
    umap_out<-paste0(outDir,"/",sample_in,".titan_epithelial.umap.pdf")
    model_outdir<-paste0(outDir,"/",sample_in,"_epithelial")
  }
  else {
    Object<-subset(Object,sample==sample_in)
    out_seurat_object<-paste0(outDir,"/",sample_in,".titan.SeuratObject.rds")
    out_titan_obj<-paste0(outDir,"/",sample_in,".titan.titanObject.rds")
    umap_out<-paste0(outDir,"/",sample_in,".titan.umap.pdf")
    elbow_out<-paste0(outDir,"/",sample_in,".titan.elbow.pdf")
    model_outdir<-paste0(outDir,"/",sample_in)
  }

  #skip titan if cell count too low
  if(sum(Object$nCount_RNA>500)<500){
    print("Cell count for RNA seq is too low...")
    saveRDS(Object,file=out_seurat_object)


  }else{
    system(paste0("mkdir -p ",model_outdir))
    Object<-subset(Object,nCount_RNA>500)
    Object[[assay]]<- as(object = Object[[assay]], Class = "Assay") #enforcing v3 assay style

    print(paste0("Finding ",as.character(nfeat)," variable features..."))
    Object <- FindVariableFeatures(
      Object, 
      selection.method = "vst",
      nfeatures = nfeat, 
      assay=assay)

    print(paste0("Setting up data for TITAN LDA..."))
    Object.sparse <- GetAssayData(Object, slot = "data", assay = assay)
    Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assay), ]
    data.use <- Matrix::Matrix(Object.sparse, sparse = T)
    data.use <- data.use * 10 #not sure if needed
    data.use <- round(data.use) #not sure if needed
    data.use <- Matrix::Matrix(data.use, sparse = T)
    sumMat <- Matrix::summary(data.use)
    cellList <- split(as.integer(data.use@i), sumMat$j)
    ValueList <- split(as.integer(sumMat$x), sumMat$j)
    cellList <- mapply(rbind, cellList, ValueList, SIMPLIFY = F)
    Genes <- rownames(data.use)
    cellList <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1, ] + 1]; x})

    print(paste0("Running ",as.character(length(topic_counts)), " topic models..."))
    lda_out<-mclapply(topic_counts, 
      function(x) 
      model_maker(topics=x,outDir=model_outdir,
        cellList=cellList,Genes=Genes,
        iterations=iterations,burnin=burnin,alpha=alpha_val,beta=beta_val), 
      mc.cores = length(topic_counts))

    files <- list.files(path = model_outdir, pattern = "Model_")
    perp_list <- NULL
    topic_numbers <- NULL
    RPC <- NULL
    files <- files[order(nchar(files), files)]

    print(paste0("Generating perplexity estimate for model selection..."))
    perp_out<-as.data.frame(do.call("rbind",
      lapply(files,function(x) RPC_calculation(model_file=x,outDir=model_outdir,data.use=data.use))))
    colnames(perp_out) <- c("Topics", "RPC")
    perp_out$Topics<-as.numeric(as.character(perp_out$Topics))
    perp_out$RPC<-as.numeric(as.character(perp_out$RPC))

    rpc_dif<-diff(perp_out$RPC)
    topic_dif<-diff(perp_out$Topics)
    perp_out$perp<-c(NA,abs(rpc_dif)/topic_dif)

    #select topics from model based on elbow
    elbow_topic<-perp_out$Topics[which(min(diff(perp_out$perp),na.rm=T)==diff(perp_out$perp))+1]
    print(paste0("Found ",as.character(elbow_topic), " topics as best topic model based on elbow plot..."))
    plt1 <- ggplot(data = perp_out, aes(x = Topics, y = RPC, group = 1)) + geom_line() + geom_point() +geom_vline(xintercept=elbow_topic,color="red")
    plt2 <- ggplot(data = perp_out, aes(x = Topics, y = perp, group = 1)) + geom_line() + geom_point()+geom_vline(xintercept=elbow_topic,color="red")
    ggsave(plt1/plt2,file=elbow_out)

    top_topics<-elbow_topic #set this up as autoselect based on elbow of plt2 in future
    top_model<-readRDS(paste0(model_outdir, "/", "Model_",as.character(top_topics),"topics.rds"))

    print(paste0("Adding topic model to Seurat Object as lda reduction..."))
    Object <- addTopicsToSeuratObject(model = top_model, Object = Object)
    #GeneDistrubition <- GeneScores(top_model)

    print(paste0("Running UMAP and clustering..."))
    Object<-RunUMAP(Object,reduction="lda",dims=1:top_topics,reduction.name="lda_umap")
    Object <- FindNeighbors(object = Object, reduction = 'lda', dims = 1:top_topics ,graph.name="lda_snn")
    Object <- FindClusters(object = Object, verbose = TRUE, graph.name="lda_snn", resolution=0.2 ) 
    print("Plotting UMAPs...")
    plt1<-DimPlot(Object,reduction="lda_umap",group.by=c("HBCA_predicted.id"))
    ggsave(plt1,file=umap_out,width=10)
    print("Done!")
    saveRDS(Object,out_seurat_object)
    saveRDS(top_model,out_titan_obj)
  }
}

lapply(25:length(unique(dat$sample)), function(x) {
sample_in=unique(dat$sample)[x]
single_sample_titan_generation(Object=dat,outDir=outdir,sample_in=sample_in,epithelial_only=TRUE) #only epithelial cells per sample
single_sample_titan_generation(Object=dat,outDir=outdir,sample_in=sample_in,epithelial_only=FALSE) #all cells per sample
})

IDC_10
