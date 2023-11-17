####Run InferCNV
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(parallel)
args = commandArgs(trailingOnly=TRUE)

dat=readRDS(arg[1])
outdir=args[2]

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels #standard seq level change threw error, using a string replace instead

infercnv_per_sample<-function(dat,outname){
  outdir_sample=paste0(outdir,"/",outname,"/infercnv/")
  system(paste0("mkdir -p ",outdir_sample))
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this

  dat$cnv_ref<-"FALSE"
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")),]$cnv_ref<-"TRUE" #set cnv ref by cell type

  #write out gene order list
  gene_order<-annotation[!duplicated(annotation$gene_name),]
  gene_order<-as.data.frame(gene_order[gene_order$gene_name %in% row.names(dat),])
  gene_order<-gene_order[c("gene_name","seqnames","start","end")]
  chrorder<-paste0("chr",c(1:22,"X","Y","M"))
  gene_order$seqnames<-factor(gene_order$seqnames,levels=chrorder) # set chr order
  gene_order<-with(gene_order, gene_order[order(seqnames, start),]) #order by chr and start position

  counts=as.matrix(dat@assays$RNA@counts[,colnames(dat)])
  cell_annotation=as.data.frame(cbind(row.names(dat@meta.data),dat@meta.data["cnv_ref"]))

  write.table(gene_order,file="inferCNV.gene_order.txt",sep="\t",col.names=F,row.names=F,quote=F)
  write.table(counts,file="inferCNV.counts.txt",sep="\t",col.names=T,row.names=T,quote=F)
  write.table(cell_annotation,file="inferCNV.annotation.txt",sep="\t",col.names=F,row.names=F,quote=F)

  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix="inferCNV.counts.txt",
    annotations_file="inferCNV.annotation.txt",
    delim="\t",
    gene_order_file="inferCNV.gene_order.txt",
    ref_group_names="TRUE")

  infercnv_obj = infercnv::run(infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=outdir_sample, 
    cluster_by_groups=TRUE, 
    denoise=TRUE,
    HMM=TRUE,
    HMM_report_by="cell",
    resume_mode=F,
    HMM_type='i3',
    num_threads=1)

  write.table(gene_order,file=paste0(outdir_sample,"/inferCNV.gene_order.txt"),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(counts,file=paste0(outdir_sample,"/inferCNV.counts.txt"),sep="\t",col.names=T,row.names=T,quote=F)
  write.table(cell_annotation,file=paste0(outdir_sample,"/inferCNV.annotation.txt"),sep="\t",col.names=F,row.names=F,quote=F)
  saveRDS(infercnv_obj,paste0(outdir_sample,"/",outname,".inferCNV.Rds"))
}

mclapply(unique(dat$sample),function(x) infercnv_per_sample(dat=dat,outname=x),mc.cores=10)
