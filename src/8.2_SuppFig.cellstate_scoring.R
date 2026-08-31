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

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

outdir="~/supp_fig"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

#plot stacked barplot of assigned celltype per sample
scsubtype_sample <- dat@meta.data %>% group_by(sample,scsubtype,Diag_MolDiag) %>% summarize(value=n()) %>% as.data.frame()
plt<-ggplot(scsubtype_sample,aes(x=sample,y=value,fill=scsubtype))+geom_bar(position="fill", stat="identity")+facet_wrap(~Diag_MolDiag,scale="free_x")
ggsave(plt,file=paste0(outdir,"/","scsubtype_assignment.pdf"))

#i think this assignment is pretty in line with the regner paper
cc.genes.updated.2019$s.genes %in% unlist(module_feats)
cc.genes.updated.2019$g2m.genes %in% unlist(module_feats)
#neither cell cycle scoring gene sets are listed in module genes

#Add cell cycle scoring
DefaultAssay(dat)<-"SCT"
dat<-CellCycleScoring(dat,
  s.features=cc.genes.updated.2019$s.genes,
  g2m.features=cc.genes.updated.2019$g2m.genes)

#wu et al.
#Add "D score" for differentiation, expecting basal-like to be less differentiated
#Expression of selected genes associated with luminal differentiation (KRT8, KRT5, KRT14, KRT19, ESR1, ERBB2)
d_score_genelist<-c("KRT8", "KRT5", "KRT14", "KRT19", "ESR1", "ERBB2") 
d_score_genelist %in% Features(dat_epi@assays$SCT)
d_scores<-base::colMeans(as.data.frame(rna_dat[row.names(rna_dat) %in% d_score_genelist,]),na.rm=TRUE)

#Add EMT score for EMT
#EMT (CDH1, CLDN3, CLDN4, CLDN7, VIM, TWIST1, SNAI1, SNAI2, ZEB1, ZEB2) 
emt_score_genelist<-c("CDH1", "CLDN3", "CLDN4", "CLDN7", "VIM", "TWIST1", "SNAI1", "SNAI2", "ZEB1", "ZEB2") 
emt_score_genelist %in% Features(dat_epi@assays$RNA)
emt_scores<-base::colMeans(as.data.frame(rna_dat[row.names(rna_dat) %in% emt_score_genelist,]),na.rm=TRUE)

dat<-AddMetaData(dat,d_scores,col.name="Wu_DScores")
dat<-AddMetaData(dat,emt_scores,col.name="Wu_EMTScores")

saveRDS(dat,file="7_merged.scsubtype.SeuratObject.rds")
dat<-readRDS(file="7_merged.scsubtype.SeuratObject.rds")

epicell_filter<-names(which(table(subset(dat,assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs"))$sample)>=30))
cancercell_filter<-names(which(table(subset(dat,assigned_celltype=="cancer")$sample)>=30))

celltype_col=c("cancer"="#c93c96",
"luminal_hs"="#603a96",
"luminal_asp"="#402b7e",
"basal_myoepithelial"="#c483b8",
"adipocyte"="#781118",
"endothelial_vascular"="#d8572a",
"endothelial_lymphatic"="#da7c27",
"pericyte"="#f7b535",
"fibroblast"="#c32f27",
"myeloid"="#6073b7",
"bcell"="#6acad5",
"plasma"="#8fd1bf",
"tcell"="#0e5169")

Idents(dat)<-factor(dat$assigned_celltype,levels=c("cancer","luminal_hs","luminal_asp","basal_myoepithelial",
"adipocyte","endothelial_vascular","endothelial_lymphatic","pericyte","fibroblast",
"myeloid","bcell","plasma","tcell"))

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

celltype_freq<-dat@meta.data %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,assigned_celltype) %>% 
  count(assigned_celltype,.drop=FALSE)

plt<-ggplot(celltype_freq, aes(fill=factor(assigned_celltype,levels=names(celltype_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=celltype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")+ theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file=paste0(outdir,"/","celltype_percbarplot.pdf"),width=10,height=10)

#barplot of scsubtype,  
scsubtype_freq<-dat@meta.data %>% 
  filter(assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs")) %>% 
  filter(sample %in% epicell_filter) %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,assigned_celltype) %>% 
  count(scsubtype,.drop=FALSE)

plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ theme(axis.text.x = element_text(angle = 90))+
  facet_grid(assigned_celltype~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")
  ggsave(plt,file=paste0(outdir,"/","scsubtype_percbarplot.pdf"),width=10,height=10)

#barplot of S/G2M cells
cellcycle_freq<-dat@meta.data %>% 
  filter(assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs")) %>% 
  filter(sample %in% epicell_filter) %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,assigned_celltype,Phase) %>% 
  count(Phase,.drop=FALSE)

plt<-ggplot(cellcycle_freq, aes(fill=factor(Phase,levels=names(Phase_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=Phase_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file=paste0(outdir,"/","phase_percbarplot.pdf"),width=10,height=10)


#and stem cell and/or TICs features (CD44, CD24, ALDH1A1, EPCAM) across the cell line database.
