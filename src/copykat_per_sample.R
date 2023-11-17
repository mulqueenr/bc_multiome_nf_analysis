library(Signac)
library(Seurat)
library(copykat)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
dat=readRDS(args[1])
outdir=args[2]

copykat_per_sample<-function(dat=dat,outname=x){
  outdir_sample=paste0(outdir,"/",outname,"/copykat")
  system(paste0("mkdir -p ",outdir_sample))

  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
  exp.rawdata <- as.matrix(dat@assays$RNA@counts)

  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")),]$cnv_ref<-"TRUE" #set cnv ref by cell type

  cnv_ref<-row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",])

  copykat_out <- copykat(rawmat=exp.rawdata, 
    KS.cut=0.15,LOW.DR=0.05,
    UP.DR=0.2,id.type="S", 
    ngene.chr=0, win.size=25, 
    sam.name=sample_name, distance="euclidean", 
    norm.cell.names=cnv_ref,output.seg="FALSE", 
    plot.genes="FALSE", genome="hg20",n.cores=1)

  saveRDS(copykat_out,paste0(outdir_sample,"/",outname,".copykat.RDS"))
}

mclapply(unique(dat$sample),function(x) copykat_per_sample(dat=dat,outname=x),mc.cores=10)

#to set CNV discrete changes, as per correspondence suggetions with Ruli Gao, 1.5x SD threshold, 1.5 absolute distance, or use +/-0.25 as cutoff