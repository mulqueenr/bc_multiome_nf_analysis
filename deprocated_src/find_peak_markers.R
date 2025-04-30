#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/harmony.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(reshape2)
library(parallel)
library(optparse)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.SeuratObject.rds"

dat=readRDS(opt$object_input)
dat<-subset(dat,Diagnosis %in% c("NAT","DCIS"))
#remove unnecessary assays, so I can better parallelize
dat[["RNA"]]<-NULL;dat[["SoupXRNA"]]<-NULL;dat[["chromvar"]]<-NULL;dat[["GeneActivity"]]<-NULL;dat[["SCT"]]<-NULL;
DefaultAssay(dat)<-"peaks"
Idents(dat)<-dat$HBCA_predicted.id

parallelized_markers<-function(obj,x){
  markers<-FindMarkers(
    obj, assay="peaks",
    slot = "data",
    ident.1 = x,
    test.use = "LR",
    only.pos=T
  )
  markers$ident<-x
return(markers)
}

out<-mclapply(unique(as.character(Idents(dat))),function(x) parallelized_markers(dat,x),mc.cores=5)


saveRDS(out,"multiome_peak_markers.rds")



