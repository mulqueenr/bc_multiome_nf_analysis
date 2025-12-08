Read in signac object and output files necessary for cistopic/scenic
```bash
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell \
--bind /home/groups/CEDAR/mulqueen/bc_multiome \
$sif

cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

```


```R

library(Seurat)
library(Signac)
library(optparse)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
dat<-subset(dat,cells=row.names(dat@meta.data)[isNA(dat@meta.data$merged_assay_clones) | dat@meta.data$merged_assay_clones != "contamination"])
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])

#run cistopic on cells with cnv clone assignment
cistopic_outdir=paste0(getwd(),"/cistopic_clones")
system(paste0("mkdir -p ", cistopic_outdir))

### RUN FOR ALL CELLS ASSIGNED CLONES###
#write out fragment paths per cell
fragment_paths<-lapply(1:length(dat@assays$ATAC@fragments),function(x) cbind(dat@assays$ATAC@fragments[[x]]@path,unlist(names(dat@assays$ATAC@fragments[[x]]@cells)),unlist(unname(dat@assays$ATAC@fragments[[x]]@cells))))
fragment_paths<-as.data.frame(do.call("rbind",fragment_paths))
row.names(fragment_paths)<-fragment_paths$V2
write.table(fragment_paths,file=paste0(cistopic_outdir,"/","frag_paths.csv"),sep=",",col.names=F,row.names=T)

#write out metadata
write.table(dat@meta.data,file=paste0(cistopic_outdir,"/","metadata.atac.csv"),col.names=T,row.names=T,sep=",")
write.table(dat@meta.data,file=paste0(cistopic_outdir,"/","metadata.rna.csv"),col.names=T,row.names=T,sep=",")

#write out ATAC counts matrix
write.table(colnames(dat@assays$ATAC@counts),file=paste0(cistopic_outdir,"/","atac_counts.cells.csv"),sep=",")
write.table(row.names(dat@assays$ATAC@counts),file=paste0(cistopic_outdir,"/","atac_counts.peaks.csv"),sep=",")
Matrix::writeMM(dat@assays$ATAC@counts,file=paste0(cistopic_outdir,"/","atac_counts.mtx"))

#write out RNA counts
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])
dat[["RNA"]]<-as(object = dat[["RNA"]], Class = "Assay")
write.table(colnames(dat@assays$RNA@counts),file=paste0(cistopic_outdir,"/","rna_counts.cells.csv"),sep=",")
write.table(row.names(dat@assays$RNA@counts),file=paste0(cistopic_outdir,"/","rna_counts.genes.csv"),sep=",")
Matrix::writeMM(dat@assays$RNA@counts,file=paste0(cistopic_outdir,"/","rna_counts.mtx"))
```

Load scenicplus singularity file for processing
```bash
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/scenicplus.sif
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/cistopic_clones
```


```python
import os
from pycisTopic.cistopic_class import *
from scipy.io import mmread
import pickle
import scanpy as sc
from pycisTopic.lda_models import run_cgs_models_mallet
import argparse


parser = argparse.ArgumentParser(
    description="Function to run single cistopic model.")

parser.add_argument("-f", "--frag_path", default = 'frag_paths.csv', help = "List of fragment locations, csv format")
parser.add_argument("-n", "--atac_counts", default = 'atac_counts.mtx', help = "Raw counts, mtx")
parser.add_argument("-c", "--atac_cells", default = 'atac_counts.cells.csv', help = "List of cells, csv")
parser.add_argument("-p", "--atac_peaks", default = 'atac_counts.peaks.csv', help = "List of peaks, csv")
parser.add_argument("-m", "--meta", default = 'metadata.atac.csv', help = "Metadata table, csv")
parser.add_argument("-o", "--outDir", default ="./", help = "Output Directory")

args = parser.parse_args()

# Project directory and files
frag_path=args.frag_path
meta=args.meta
atac_counts=args.atac_counts
atac_cells=args.atac_cells
atac_peaks=args.atac_peaks
outDir=args.outDir

# Create cisTopic object
atac_counts = mmread(atac_counts)
atac_peaks =  pd.read_csv(atac_peaks)
atac_peaks = [peak.replace("-",":",1) for peak in atac_peaks['x']]  #reformat name
atac_cells =  pd.read_csv(atac_cells)
cell_data =  pd.read_csv(meta,dtype="string")

cistopic_obj = create_cistopic_object(fragment_matrix=atac_counts.tocsr(),
    cell_names=atac_cells['x'].tolist(),
    region_names=atac_peaks,
    tag_cells=False)

# Adding cell information
cistopic_obj.add_cell_data(cell_data)
pickle.dump(
    cistopic_obj,
    open(os.path.join(outDir, "scenicplus_"+args.atac_counts.split('_')[0]+"_cistopic_obj.pkl"), "wb")
)

parser.add_argument("-m", "--memory", required = False, default = '400G', help = "Memory for MALLET_MEMORY")
parser.add_argument("-t", "--taskCpus", required = False, default = 10, help = "CPUS to run")
parser.add_argument("-M", "--mallet", required = False, default ="/container_mallet/bin/mallet", help = "Mallet path, built into container")

args = parser.parse_args()

#tmp dir
tmpDir = './'

# Run models with mallet
os.environ['MALLET_MEMORY'] = '300G'

# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=list(range(5, 40, 5)), 
    n_cpu=int(10),
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=tmpDir,
    save_path=tmpDir,
    mallet_path='/container_mallet/bin/mallet',
)


```

```R
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(dplyr)
library(optparse)
library(org.Hs.eg.db)
library(dendextend)
library(msigdbr,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(fgsea,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(grImport,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(cisTopic)
library(SeuratWrappers)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
set.seed(1234)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

dat<-subset(dat,merged_assay_clones != "contamination")
dat<-subset(dat,cells=names(which(!is.na(dat$merged_assay_clones))))
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])

#run cistopic on cells with cnv clone assignment
cistopic_outdir=paste0(getwd(),"/cistopic_clones")
system(paste0("mkdir -p ", cistopic_outdir))
out_seurat_object<-paste0(cistopic_outdir,"/","9_merged.cnv_clones_cistopic.SeuratObject.rds")
model_selection_out<-paste0(cistopic_outdir,"/","cistopic.model_selection.pdf")
out_cistopic_obj<-paste0(cistopic_outdir,"/","cistopic.cistopicObject.rds")
umap_out<-paste0(cistopic_outdir,"/","cistopic.umap.pdf")
cistopic_counts_frmt<-dat@assays$ATAC@counts
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


#run titan on cells with cnv clone assignment
model_maker <- function(topics,outDir,cellList,Genes,iterations,burnin,alpha,beta) {
    selected.Model <- lda::lda.collapsed.gibbs.sampler(
      cellList, topics, Genes, 
      num.iterations = iterations, 
      alpha = alpha, eta = beta, 
      compute.log.likelihood = TRUE, 
      burnin = burnin)[-1]
      saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
    }

RPC_calculation<-function(model_file,outDir,data.use,topic_numbers) {
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
  assay="RNA",
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
    Object<-subset(Object,assigned_celltype %in% c("cancer_luminal_epithelial","luminal_epithelial","basal_epithelial"))
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
  if(sum(Object$nCount_RNA>500)<200){
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
      lapply(files,function(x) RPC_calculation(model_file=x,outDir=model_outdir,data.use=data.use,topic_numbers=topic_numbers))))
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

#run pca on cnvs

#correlate
