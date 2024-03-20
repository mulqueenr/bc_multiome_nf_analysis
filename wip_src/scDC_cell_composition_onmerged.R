#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
#cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects

library(scDC)
library(Signac)
library(Seurat)
library(optparse)
library(dplyr) 
library(reshape2)
library(ggplot2)

dat<-readRDS("merged.public_transfer.SeuratObject.rds")

subject<-dat$sample
cellTypes<-dat$HBCA_predicted.id

res_scDC_noClust <- scDC_noClustering(cellTypes, subject, calCI = TRUE, calCI_method = c("percentile", "BCa", "multinom"),nboot = 1000,ncores=20)

cond_dat<-unique(data.frame(subject=dat$sample,diagnosis=dat$Diagnosis,mol_diagnosis=dat$Mol_Diagnosis))
cond_dat<-cond_dat[match(unlist(unique(res_scDC_noClust$info[2])),cond_dat$subject),]
cond<-paste(cond_dat$diagnosis,cond_dat$mol_diagnosis)
out<-fitGLM(res_scDC_noClust,cond,subject_effect = TRUE, pairwise = TRUE,fixed_only = FALSE, verbose = TRUE)


pdf("scDC_test.pdf",width=10,height=10)
barplotCI(res_scDC_noClust, cond)
dev.off()



#Set up metadata and set up facet labels as factors for ordering
metadat<-as.data.frame(dat@meta.data)
metadat$diagnosis = factor(metadat$Diagnosis, levels=c("NAT","DCIS","IDC","ILC"), labels=c("NAT","DCIS","IDC","ILC")) 
metadat$molecular_type = factor(metadat$Mol_Diagnosis, levels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2+","ER+/PR-/HER2-"), labels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2+","ER+/PR-/HER2-")) 

#Cells PF
metadat$epi<-"Nonepi"
metadat[metadat$HBCA_predicted.id %in% c("basal cell" ,"luminal epithelial cell of mammary gland"),]$epi<-"Epi"
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,epi) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=epi,y=n))+geom_bar(stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free") #+ scale_y_continuous(trans='log10')
ggsave(plt1,file="barplot_qc_cellcount.pdf")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,HBCA_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=HBCA_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="HBCA_barplot_qc_celltype.pdf")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,swarbrick_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=swarbrick_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="swarbrick_barplot_qc_celltype.pdf")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(diagnosis, molecular_type,sample,EMBO_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=EMBO_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="Embo_barplot_qc_celltype.pdf")

#Cell types (Epi excluded) (stacked bar)
DF<-as.data.frame(metadat %>% filter(!(HBCA_predicted.id %in% c("basal cell" ,"luminal epithelial cell of mammary gland"))) %>% group_by(diagnosis, molecular_type,sample,HBCA_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=HBCA_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="HBCA_barplot_qc_celltype.nonepi.pdf")

#Cell types (Epi excluded) (stacked bar)
DF<-as.data.frame(metadat %>% filter(!(EMBO_predicted.id %in% c("epithelial", "cycling.epithelial"))) %>% group_by(diagnosis, molecular_type,sample,EMBO_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=EMBO_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~diagnosis+molecular_type,scales="free_x",space="free")
ggsave(plt1,file="Embo_barplot_qc_celltype.nonepi.pdf")

