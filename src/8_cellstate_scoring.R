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
  make_option(c("-i", "--object_input"), type="character", default="6_merged.celltyping.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

outdir="/home/groups/MohammedLab/bc_multiome/suppfig2"

if (!dir.exists(outdir)) {
  dir.create(outdir)
}

dat_epi<-subset(dat,assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs"))
dat_epi[["RNA"]]<-JoinLayers(dat_epi[["RNA"]])
dat_epi<-SCTransform(dat_epi)

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
module_feats[["Basal_SC"]][!module_feats[["Basal_SC"]] %in% Features(dat_epi@assays$SCT)]
module_feats[["Her2E_SC"]][!module_feats[["Her2E_SC"]] %in% Features(dat_epi@assays$SCT)]
module_feats[["LumA_SC"]][!module_feats[["LumA_SC"]] %in% Features(dat_epi@assays$SCT)]
module_feats[["LumB_SC"]][!module_feats[["LumB_SC"]] %in% Features(dat_epi@assays$SCT)]
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

rna_dat<-GetAssayData(dat_epi,assay="RNA",layer="data")

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

#i think this assignment is pretty in line with the regner paper
cc.genes.updated.2019$s.genes %in% unlist(module_feats)
cc.genes.updated.2019$g2m.genes %in% unlist(module_feats)
#neither cell cycle scoring gene sets are listed in module genes
#which is what we want so the metrics arent conflated


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

