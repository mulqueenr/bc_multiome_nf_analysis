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
setseed(123)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="6_merged.celltyping.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character"),
  make_option(c("-r", "--ref_object"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri/nakshatri_multiome.geneactivity.rds", 
              help="Nakshatri reference object for epithelial comparisons", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

dat_epi<-subset(dat,assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs"))
dat_epi[["RNA"]]<-JoinLayers(dat_epi[["RNA"]])

#downloading gene list from SCSubtype git repo
system("wget https://raw.githubusercontent.com/Swarbricklab-code/BrCa_cell_atlas/main/scSubtype/NatGen_Supplementary_table_S4.csv")
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv",col.names=c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))

module_feats<-list()
module_feats[["Basal_SC"]]=as.vector(sigdat[,"Basal_SC"])
module_feats[["Her2E_SC"]]=as.vector(sigdat[,"Her2E_SC"])
module_feats[["LumA_SC"]]=as.vector(sigdat[,"LumA_SC"])
module_feats[["LumB_SC"]]=as.vector(sigdat[,"LumB_SC"])
module_feats<-lapply(module_feats,function(x) {x[x!=""]}) #remove empty
module_feats<-lapply(module_feats,function(x) {unlist(lapply(x, function(gene) gsub(gene,pattern=".",replace="-",fixed=TRUE)))}) #correct syntax

#find genes with name changes
module_feats[["Basal_SC"]][!module_feats[["Basal_SC"]] %in% Features(dat_epi@assays$RNA)]
module_feats[["Her2E_SC"]][!module_feats[["Her2E_SC"]] %in% Features(dat_epi@assays$RNA)]
module_feats[["LumA_SC"]][!module_feats[["LumA_SC"]] %in% Features(dat_epi@assays$RNA)]
module_feats[["LumB_SC"]][!module_feats[["LumB_SC"]] %in% Features(dat_epi@assays$RNA)]
#limit to protein coding genes

rename_genes<-c(
  "STRA13"="CENPX",
  "AIM1"="CRYBG1",
  "C17orf89"="NDUFAF8",
  "ATP5I"="ATP5ME",
  "ENTHD2"="TEPSIN",
  "ATP5C1"="ATP5F1C",
  "C6orf203"="MTRES1",
  "C6orf48"="SNHG32",
  "MYEOV2"="COPS9",
  "MLLT4"="AFDN")
#"RP1-60O19-1",
#"RP11-206M11-7",  
#"GPX1"
#"AP000769-1"


#rename some gene symbols
module_feats<-lapply(module_feats,function(i){
    unlist(lapply(i,function(gene){ifelse(gene %in% names(rename_genes),yes=rename_genes[gene],no=gene)}
))}) #rename select genes



#From Regner et al.
#"To assign a subtype call to a cell, 
#we calculated the average (that is, the mean) 
#read counts for each of the four signatures for each cell. 
#The SC subtype with the highest signature score was then assigned to each cell."

rna_dat<-GetAssayData(dat_epi,assay="RNA",layer="scale.data")

scsubtype_scores<-lapply(module_feats,function(scsubtype){
  base::colMeans(as.data.frame(rna_dat[row.names(rna_dat) %in% scsubtype,]),na.rm=TRUE)
})
names(scsubtype_scores)<-paste0("SC_Subtype_",c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))
scsubtype_scores<-as.data.frame(scsubtype_scores)

scsubtype_scores <- scsubtype_scores %>%
  mutate(scsubtype = case_when(
    SC_Subtype_Basal_SC == pmax(SC_Subtype_Basal_SC, SC_Subtype_Her2E_SC, SC_Subtype_LumA_SC, SC_Subtype_LumB_SC) ~ "SC_Subtype_Basal_SC",
    SC_Subtype_Her2E_SC == pmax(SC_Subtype_Basal_SC, SC_Subtype_Her2E_SC, SC_Subtype_LumA_SC, SC_Subtype_LumB_SC) ~ "SC_Subtype_Her2E_SC",
    SC_Subtype_LumA_SC ==  pmax(SC_Subtype_Basal_SC, SC_Subtype_Her2E_SC, SC_Subtype_LumA_SC, SC_Subtype_LumB_SC) ~ "SC_Subtype_LumA_SC",
     SC_Subtype_LumB_SC ==  pmax(SC_Subtype_Basal_SC, SC_Subtype_Her2E_SC, SC_Subtype_LumA_SC, SC_Subtype_LumB_SC) ~ "SC_Subtype_LumB_SC"
  ))
dat<-AddMetaData(dat,scsubtype_scores)

#dat_epi<-AddModuleScore(dat_epi,
#  module_feats,
#  assay = "RNA",
#  name = paste0("SC_Subtype_",c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")),
#  search = TRUE)


#plot stacked barplot of assigned celltype per sample
scsubtype_sample <- dat@meta.data %>% group_by(sample,scsubtype,Diag_MolDiag) %>% summarize(value=n()) %>% as.data.frame()
plt<-ggplot(scsubtype_sample,aes(x=sample,y=value,fill=scsubtype))+geom_bar(position="fill", stat="identity")+facet_wrap(~Diag_MolDiag,scale="free_x")
ggsave(plt,file="scsubtype_assignment.pdf")

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

#use only clones with at least 50 cell 
clone_filter<-names(which(table(dat$merged_assay_clones)>=50))
clone_filter<-clone_filter[clone_filter!="normal"]

cancercell_filter<-names(which(table(subset(dat,assigned_celltype=="cancer")$sample)>=50))

epicell_filter<-names(which(table(subset(dat,assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs"))$sample)>=50))

#barplot of cell types across samples
celltype_col=c("cancer"="#9e889e",
"luminal_hs"="#4c3c97",
"luminal_asp"="#7161ab",
"basal_myoepithelial"="#ee6fa0",
"adipocyte"="#af736d",
"endothelial_vascular"="#72c8f1",
"endothelial_lymphatic"="#b8dca5",
"pericyte"="#edb379",
"fibroblast"="#e12228",
"myeloid"="#239ba8",
"bcell"="#243d97",
"plasma"="#742b8c",
"tcell"="#003147")

celltype_freq<-dat@meta.data %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,assigned_celltype) %>% 
  count(assigned_celltype,.drop=FALSE)

plt<-ggplot(celltype_freq, aes(fill=factor(assigned_celltype,levels=names(celltype_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=celltype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="celltype_percbarplot.pdf",width=10,height=10)

#barplot of scsubtype,  
scsubtype_col=c(
  "SC_Subtype_Basal_SC"="#da3932",
  "SC_Subtype_Her2E_SC"="#f0c2cb",
  "SC_Subtype_LumA_SC"="#2b2c76",
  "SC_Subtype_LumB_SC"="#86cada"
)

#per sample 50 epi cells minimum
scsubtype_freq<-dat@meta.data %>% 
  filter(assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs")) %>% 
  filter(sample %in% epicell_filter) %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,assigned_celltype) %>% 
  count(scsubtype,.drop=FALSE)

plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ 
  facet_grid(assigned_celltype~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")
  ggsave(plt,file="scsubtype_percbarplot.pdf",width=10,height=10)

#per clone, 50 cells per clone minimum
scsubtype_freq<-dat@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) %>% 
  count(scsubtype,.drop=TRUE)

plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="scsubtype_perclone_percbarplot.pdf",width=10,height=10)

#barplot of S/G2M cells
#per clone, 50 cancer cells minimum
Phase_col=c(
  "G1"="#e5f5e0",
  "S"="#a1d99b",
  "G2M"="#00441b"
)

cellcycle_freq<-dat@meta.data %>% 
  filter(assigned_celltype=="cancer") %>%
  filter(sample %in% cancercell_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,Phase,merged_assay_clones) %>% 
  count(Phase,.drop=TRUE)

plt<-ggplot(cellcycle_freq, aes(fill=factor(Phase,levels=names(Phase_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=Phase_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="phase_percbarplot.pdf",width=10,height=10)


cellcycle_freq<-dat@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) %>% 
  count(Phase,.drop=TRUE)

plt<-ggplot(cellcycle_freq, aes(fill=factor(Phase,levels=names(Phase_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=Phase_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="phase_perclone_percbarplot.pdf",width=10,height=10)


#box plot of wu scores
emt_and_d_scores<-dat@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>%
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,merged_assay_clones) 

plt<-ggplot(emt_and_d_scores, aes(y=Wu_EMTScores, x=merged_assay_clones)) + 
  geom_boxplot(width=1) +
  geom_jitter(width = 1) +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="emtscore_perclone_boxplot.pdf",width=10,height=10)


plt<-ggplot(emt_and_d_scores, aes(y=Wu_DScores, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="dscore_perclone_boxplot.pdf",width=10,height=10)

#and stem cell and/or TICs features (CD44, CD24, ALDH1A1, EPCAM) across the cell line database.

