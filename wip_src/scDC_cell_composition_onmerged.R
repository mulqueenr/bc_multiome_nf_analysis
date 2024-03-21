#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
library(scDC)
library(Signac)
library(Seurat)
library(optparse)

setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
dat<-readRDS("merged.public_transfer.SeuratObject.rds")

subject<-dat$sample
cellTypes<-dat$HBCA_predicted.id

res_scDC_noClust <- scDC_noClustering(cellTypes, subject, calCI = TRUE, calCI_method = c("percentile", "BCa", "multinom"),nboot = 1000,ncores=20)

cond_dat<-unique(data.frame(subject=dat$sample,diagnosis=dat$Diagnosis,mol_diagnosis=dat$Mol_Diagnosis))
cond_dat<-cond_dat[match(unlist(unique(res_scDC_noClust$info[2])),cond_dat$subject),]
cond<-paste(cond_dat$diagnosis,cond_dat$mol_diagnosis)
out<-fitGLM(res_scDC_noClust,cond,subject_effect = TRUE, pairwise = TRUE,fixed_only = FALSE, verbose = TRUE)


#function 'chm_factor_ldetL2' not provided by package 'Matrix'