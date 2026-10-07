#CONVERT THIS TO CLONE LEVEL ANALYSIS INSTEAD
#mulqueen@arc-infra-3
#srun --partition=interactive --cpus-per-task=30 --time=12:00:00 --mem=100G --nodes=1 --pty /bin/bash
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
#singularity shell \
#--bind /home/groups/CEDAR/mulqueen/bc_multiome \
#$sif
#cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

library(Seurat)
library(Signac)
library(ggplot2)
library(optparse)
library(dplyr)
library(ComplexHeatmap)
library(dendextend)
library(ggdendro)
library(circlize)
library(ggtern)
set.seed(1234)
setwd("/home/groups/MohammedLab/bc_multiome/seurat_object")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)


dat<-JoinLayers(dat,assay="RNA")
dat_sub<-subset(dat,cells=names(which(!is.na(dat$merged_assay_clones))))
dat_sub<-subset(dat_sub,cells=names(which(dat_sub$merged_assay_clones != "contamination")))
dat_sub<-subset(dat_sub,cells=names(which(dat_sub$assigned_celltype %in% c("cancer"))))
dat_sub<-subset(dat_sub,cells=row.names(dat_sub@meta.data)[which(!endsWith(dat_sub$merged_assay_clones,suffix="_normal"))])

clone_filter<-names(which(table(dat_sub$merged_assay_clones)>=30))
dat_sub<-subset(dat_sub,cells=row.names(dat_sub@meta.data)[dat_sub$merged_assay_clones %in% clone_filter]) #limit to cells passing cnv
table(dat_sub$merged_assay_clones,dat_sub$sample)

outdir="/home/groups/MohammedLab/bc_multiome/suppfig2"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}


scsubtype_col=c(
  "SC_Subtype_Basal_SC"="#da3932",
  "SC_Subtype_Her2E_SC"="#f0c2cb",
  "SC_Subtype_LumA_SC"="#2b2c76",
  "SC_Subtype_LumB_SC"="#86cada"
)

Phase_col=c(
  "G1"="#e5f5e0",
  "S"="#a1d99b",
  "G2M"="#00441b"
)

#barplot of scsubtype,  
scsubtype_freq<-dat_sub@meta.data %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,assigned_celltype,merged_assay_clones) %>%
  dplyr::count(scsubtype,.drop=FALSE)

plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ theme(axis.text.x = element_text(angle = 90))+
  facet_grid(assigned_celltype~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")
  ggsave(plt,file=paste0(outdir,"/","scsubtype_percbarplot.pdf"),width=10,height=10)

#barplot of S/G2M cells
cellcycle_freq<-dat_sub@meta.data %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,assigned_celltype,Phase,merged_assay_clones) %>% 
  dplyr::count(Phase,.drop=FALSE)

plt<-ggplot(cellcycle_freq, aes(fill=factor(Phase,levels=names(Phase_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=Phase_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file=paste0(outdir,"/","phase_percbarplot.pdf"),width=10,height=10)


#and stem cell and/or TICs features (CD44, CD24, ALDH1A1, EPCAM) across the cell line database.

#box plot of wu scores
emt_and_d_scores<-dat_sub@meta.data 

plt<-ggplot(emt_and_d_scores, aes(y=Wu_EMTScores, x=merged_assay_clones, fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file=paste0(outdir,"/","emtscore_perclone_boxplot.pdf"),width=10,height=10)


plt<-ggplot(emt_and_d_scores, aes(y=Wu_DScores, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file=paste0(outdir,"/","dscore_perclone_boxplot.pdf"),width=10,height=10)


#check esr1 expression on clones
dat_sub<-SCTransform(dat_sub)

motif<-data.frame(gene_name=ConvertMotifID(object = dat_sub, assay="ATAC",id = row.names(dat_sub[["chromvar"]]@data)), motif_name=row.names(dat_sub[["chromvar"]]@data))
gene_expression<-FetchData(object = dat_sub, vars = c("ESR1","KRT5","PGR"), assay="SCT")
chromatin_expression_esr1<-FetchData(object = dat_sub, vars = c("MA0112.3"), assay="chromvar")
chromatin_expression_pgr<-FetchData(object = dat_sub, vars = c("MA0113.3"), assay="chromvar") #NR3C1 is a AR motif, standard proxy for PGR
dat_tmp<-AddMetaData(dat_sub,gene_expression)
dat_tmp<-AddMetaData(dat_tmp,chromatin_expression_esr1,col.name="ESR1_TF")
dat_tmp<-AddMetaData(dat_tmp,chromatin_expression_pgr,col.name="PGR_TF")

esr1_clones<-dat_tmp@meta.data %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) 


plt1<-ggplot(esr1_clones, aes(y=ESR1, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

plt2<-ggplot(esr1_clones, aes(y=ESR1_TF, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")


plt3<-ggplot(esr1_clones, aes(y=PGR, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

plt4<-ggplot(esr1_clones, aes(y=PGR_TF, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt1/plt2/plt3/plt4,file=paste0(outdir,"/","ESR1_perclone_boxplot.pdf"),width=10,height=10)

