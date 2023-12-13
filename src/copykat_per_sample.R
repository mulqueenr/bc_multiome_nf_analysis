library(Signac)
library(Seurat)
library(copykat)

args = commandArgs(trailingOnly=TRUE)

dat=readRDS(args[1]) #dat=readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")
sample_arr=as.numeric(as.character(args[2])) #sample_arr=as.numeric(as.character(8)) 
copykat_per_sample<-function(dat=dat,outname=x){
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
  exp.rawdata <- as.matrix(dat@assays$RNA@counts)

  DefaultAssay(dat)<-"RNA"
  dat$cnv_ref<-"FALSE"
  dat@meta.data[!(dat$HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")),]$cnv_ref<-"TRUE" #set cnv ref by cell type

  cnv_ref<-row.names(dat@meta.data[dat@meta.data$cnv_ref=="TRUE",])

  copykat_out <- copykat(rawmat=exp.rawdata, 
    KS.cut=0.15,LOW.DR=0.05,
    UP.DR=0.2,id.type="S", 
    ngene.chr=0, win.size=25, 
    sam.name=outname, distance="euclidean", 
    norm.cell.names=cnv_ref,output.seg="FALSE", 
    plot.genes="FALSE", genome="hg20",n.cores=5)

  saveRDS(copykat_out,paste0(outname,".copykat.rds"))
}

copykat_per_sample(dat=dat,outname=unique(dat$sample)[sample_arr])
#to set CNV discrete changes, as per correspondence suggetions with Ruli Gao, 1.5x SD threshold, 1.5 absolute distance, or use +/-0.25 as cutoff

"""
#alternative on single node (three jobs at once):
arr_in=$(seq 1 19)
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
cd ${proj_dir}/nf_analysis/cnv_analysis/copykat
src_dir=${proj_dir}"/src"
obj="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"
parallel -j 1 Rscript ${src_dir}/copykat_per_sample.R $obj {} ::: $arr_in
"""