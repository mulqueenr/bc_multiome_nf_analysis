library(Signac)
library(Seurat)
set.seed(1234)
library(ggplot2)
library(Matrix)
library(rliger)
library(SeuratWrappers)
library(parallel)
library(GenomicRanges)
library(patchwork)
library(optparse)
#Following http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html 
#Using a parallelized Signac GeneActivity function for the scATAC.
#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character"),
    make_option(c("-o","--plot_output_directory"), type="character", default=NULL,
              help="Directory to publish output plots to.", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt$object_input<-"merged.public_transfer.SeuratObject.rds"
#opt$plot_output_directory<-"."
dat=readRDS(opt$object_input)
outname<-strsplit(opt$object_input,"[.]")[1]
outdir<-opt$plot_output_directory


#RNA liger
rna_liger<-function(nfeat=1000,dims=10,k_in=10){
  DefaultAssay(dat)<-"RNA"
  dat <- NormalizeData(dat)
  dat <- FindVariableFeatures(dat,nfeatures=nfeat)
  dat <- ScaleData(dat, split.by = "sample", do.center = FALSE)
  dat <- RunOptimizeALS(dat, k = k_in, lambda = 5, split.by = "sample")
  dat <- RunQuantileNorm(dat, split.by = "sample")
  dat <- FindNeighbors(dat, reduction = "iNMF", dims = seq(1,dims,1))
  dat <- FindClusters(dat, resolution = 0.3)
  dat <- RunUMAP(dat, dims = 1:ncol(dat[["iNMF"]]), reduction = "iNMF")
  plt<-DimPlot(dat, group.by = c("sample", "HBCA_predicted.id","diagnosis","EMBO_predicted.id"), ncol = 2)+ggtitle(paste("nfeat:",as.character(nfeat),"dim:",as.character(dims),"k:",as.character(k_in)))
  ggsave(plt,
    file=paste0(outdir,"/","merged.liger.RNA.",as.character(nfeat),".",as.character(dims),".",as.character(k_in),".pdf"),
    width=20,height=20)
}

#from https://github.com/stuart-lab/signac/blob/HEAD/R/utilities.R
CollapseToLongestTranscript <- function(ranges) {
  range.df <- data.table::as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c("gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name")
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(df = collapsed, keep.extra.columns = TRUE)
  return(gene.ranges)
}

#from https://github.com/stuart-lab/signac/blob/HEAD/R/utilities.R
Extend <- function(x,upstream = 0,downstream = 0,from.midpoint = FALSE) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(test = on_plus, yes = upstream, no = downstream)
    new_end <- midpoints + ifelse(test = on_plus, yes = downstream, no = upstream)
  } else {
    new_start <- start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
    new_end <- end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}

split_gene_count<-function(x,obj=dat){
    DefaultAssay(obj)<-"peaks"
    FeatureMatrix(fragments = Fragments(obj),cells=Cells(obj),
                              features= feat_split[[x]],
                              verbose = TRUE,
                              process_n=20000)
}

#GA liger
GA_liger<-function(nfeat=1000,dims=10,k_in=10){
  DefaultAssay(dat)<-"GeneCount"
  dat <- NormalizeData(dat)
  dat <- FindVariableFeatures(dat,nfeatures=nfeat)
  dat <- ScaleData(dat, split.by = "sample", do.center = FALSE)
  dat <- RunOptimizeALS(dat, k = k_in, lambda = 5, split.by = "sample")
  dat <- RunQuantileNorm(dat, split.by = "sample")
  dat <- FindNeighbors(dat, reduction = "iNMF", dims = seq(1,dims,1))
  dat <- FindClusters(dat, resolution = 0.3)
  dat <- RunUMAP(dat, dims = 1:ncol(dat[["iNMF"]]), reduction = "iNMF")
  plt<-DimPlot(dat, group.by = c("sample", "HBCA_predicted.id","diagnosis","EMBO_predicted.id"), ncol = 2)+ggtitle(paste("nfeat:",as.character(nfeat),"dim:",as.character(dims),"k:",as.character(k_in)))
  ggsave(plt,file=paste0(outdir,"/","merged.liger.GA.",as.character(nfeat),".",as.character(dims),".",as.character(k_in),".pdf"),width=30)
}

peak_liger<-function(nfeat=1000,dims=10,k_in=10){
  #following scIB filtering and no scaling of peaks
  DefaultAssay(dat)<-"peaks"
  dat<-BinarizeCounts(dat) # binarize peak accessibility
  dat<-FindTopFeatures(dat,assay="peaks",min.cutoff=200)
  dat <- FindVariableFeatures(dat,nfeatures=nfeat)
  dat<-SetAssayData(dat,assay="peaks",slot="scale.data",new.data=as.matrix(dat@assays$peaks@data[dat@assays$peaks@var.features,]))
  dat <- RunOptimizeALS(dat, k = k_in, lambda = 5, split.by = "sample")
  dat <- RunQuantileNorm(dat, split.by = "sample")
  dat <- FindNeighbors(dat, reduction = "iNMF", dims = seq(1,dims,1))
  dat <- FindClusters(dat, resolution = 0.3)
  dat <- RunUMAP(dat, dims = 1:ncol(dat[["iNMF"]]), reduction = "iNMF")
  plt<-DimPlot(dat, group.by = c("sample", "HBCA_predicted.id","diagnosis"), ncol = 2)+ggtitle(paste("nfeat:",as.character(nfeat),"dim:",as.character(dims),"k:",as.character(k_in)))
  ggsave(plt,file=paste0(outdir,"/","merged.liger.peaks.",as.character(nfeat),".",as.character(dims),".",as.character(k_in),".pdf"),width=20,height=20)
}


#http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html

RNA_and_GA_liger<-function(dat,nfeat_rna=1000,nfeat_peaks=1000,dim_in=10,k_in=10,epithelial_only=FALSE){
  dat_in<-dat
  if(epithelial_only){
    dat_in<-subset(dat_in,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
    out_liger_umap<-paste0(outdir,"/","merged.liger.epithelial_only",,as.character(dim_in),".",as.character(nfeat_peaks),".",as.character(nfeat_rna),".RNA_and_GA.pdf")
  } else {
    out_liger_umap<-paste0(outdir,"/","merged.liger.",as.character(dim_in),".",as.character(nfeat_peaks),".",as.character(nfeat_rna),".RNA_and_GA.pdf")
  }
  dat_in<-SetAssayData(dat_in,assay="liger_in",slot="scale.data",new.data=as.matrix(dat_in@assays$liger_in@counts))
  DefaultAssay(dat_in)<-"liger_in"
  #filter out batches with too few of cells for correction (require at least 100)
  dat_in <- subset(dat_in, cells=colnames(dat_in)[which(dat_in$sample %in% names(which(table(dat_in$sample)>100)))])
  dat_in <- RunOptimizeALS(dat_in, k = k_in, lambda = 5, split.by = "sample")
  dat_in <- RunQuantileNorm(dat_in, split.by = "sample")
  dat_in <- FindNeighbors(dat_in, reduction = "iNMF", dims = seq(1,dim_in,1))
  dat_in <- FindClusters(dat_in, resolution = 0.3)
  dat_in <- RunUMAP(dat_in, dims = 1:ncol(dat_in[["iNMF"]]), reduction = "iNMF")
  plt<-DimPlot(dat_in, group.by = c("sample", "EMBO_predicted.id","Diagnosis","HBCA_predicted.id"), ncol = 2)+
  ggtitle(paste("nfeat_rna:",as.character(nfeat_rna),
      "nfeat_ga:",as.character(nfeat_peaks),
      "dim:",as.character(dim_in),
      "k:",as.character(k_in)))
  ggsave(plt,file=out_liger_umap,width=20,height=20)
  return(dat_in)
}



#plot label transfer predictions over integration
plot_predictions<-function(dat=dat_in,ref_prefix){
  out_plot<-paste0(ref_prefix,".predictions.umap.pdf")
  feat_predictions=colnames(dat_in@meta.data)
  feat_predictions=feat_predictions[startsWith(feat_predictions,prefix=paste0(ref_prefix,"_prediction."))]
  feat_predictions=feat_predictions[!endsWith(feat_predictions,c("max"))]

  plt1<-FeaturePlot(dat,features=feat_predictions,pt.size=0.1,order=T,col=c("white","red"))
  plt2<-DimPlot(dat,group.by=paste0(ref_prefix,'_predicted.id'),pt.size=0.5)
  plt3<-DimPlot(dat,group.by='sample',pt.size=0.5)

  plt<-(plt2|plt3)/plt1
  ggsave(plt,file=out_plot,width=20,height=30,limitsize=F)
  print(paste("Plotted predictions for:",ref_prefix))
}

#doing feature setting following signac gene activity calculation
#filter to protein coding
#subset to which protein coding gene is longest that has the same name
#extend 2kb upstream for promoter inclusion
#feat=dat@assays$peaks@annotation[dat@assays$peaks@annotation$gene_biotype=="protein_coding",]
#feat<-mclapply(unique(feat$gene_name),function(x) CollapseToLongestTranscript(feat[feat$gene_name==x,]),mc.cores=5) #collapse to longest transcripts
#feat<-unlist(as(feat, "GRangesList"))
#feat<-setNames(feat,feat$gene_name)#set row names as gene names
#feat<-feat[feat@ranges@width<500000,]#filter extra long transcripts
#transcripts <- Extend(x = feat,upstream = 2000,downstream = 0)# extend to include promoters
#feat_split<-split(transcripts, rep_len(1:300, length(transcripts)))#parallelize gene count to speed up feature matrix generation

#dat_atac_counts<-mclapply(1:length(feat_split),split_gene_count,mc.cores=5)
#x<-do.call("rbind",dat_atac_counts)
#dat_atac_counts<-x
#saveRDS(dat_atac_counts,file=paste0(outdir,"/","merged.genecounts.rds"))
#dat[['GeneCount']] <- CreateAssayObject(counts = dat_atac_counts)

DefaultAssay(dat)<-"RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat,nfeatures=nfeat_rna)
dat <- ScaleData(dat, split.by = "sample", do.center = FALSE)

DefaultAssay(dat)<-"GeneActivity"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat,nfeatures=nfeat_peaks)
dat <- ScaleData(dat, split.by = "sample", do.center = FALSE)

ga<-dat[["GeneActivity"]]@scale.data
dat[["RNA"]] <- as(object = dat[["RNA"]], Class = "Assay")
rna<-dat@assays$RNA$scale.data
row.names(ga)<-paste0("GA_",row.names(ga))
row.names(rna)<-paste0("RNA",row.names(rna))
merged_dat<-as.matrix(rbind(ga,rna))

dat[["liger_in"]]<-CreateAssayObject(counts = merged_dat)

for(i in c(5000,10000)){
  for(j in c(5000,10000)){
  for(k in c(20,30,50)){
    RNA_and_GA_liger(dat=dat,nfeat_rna=i,nfeat_peaks=j,dim_in=k,k_in=k,epithelial_only=FALSE)
    }
  }
}
#loops for testing
#rna liger: nfeat 5000, dim 30 and k 30 seems to have the best cell type separation
#for(i in c(1000,2000,5000,10000)){for(j in c(10,20,30)){for(k in c(10,20,30,50)){if(k>=j){rna_liger(nfeat=i,dims=j,k_in=k)}}}}
#ga liger: #nfeat 10000, dim 20 and k 20 seems to have the best cell type separation (not as good as RNA)
#for(i in c(1000,2000,5000,10000)){for(j in c(10,20,30)){for(k in c(10,20,30,50)){if(k>=j){GA_liger(nfeat=i,dims=j,k_in=k)}}}}
#for(i in c(1000,2000,5000,10000)){for(k in c(10,20,30,50)){peak_liger(nfeat=i,dims=k,k_in=k)}}


k=50
dat_in<-RNA_and_GA_liger(nfeat_rna=10000,nfeat_peaks=10000,dim_in=k,k_in=k,epithelial_only=FALSE)
saveRDS(dat_in,file="merged.liger.SeuratObject.rds")
lapply(c("swarbrick","EMBO","HBCA"), function(x) plot_predictions(dat_in,ref_prefix=x))

dat_in<-RNA_and_GA_liger(nfeat_rna=10000,nfeat_peaks=10000,dim_in=k,k_in=k,epithelial_only=TRUE)
saveRDS(dat_in,file="merged.liger.epithelial.SeuratObject.rds")
lapply(c("swarbrick","EMBO","HBCA"), function(x) plot_predictions(dat_in,ref_prefix=x))
