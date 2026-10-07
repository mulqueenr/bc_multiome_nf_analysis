

library(Signac)
library(Seurat)
library(circlize)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
library(plyr)
library(optparse)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(harmony) #local
library(dendextend)
set.seed(1234)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="6_merged.celltyping.SeuratObject.rds", 
              help="Input seurat object", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat <- readRDS(file=opt$object_input)
dir.create("/home/groups/MohammedLab/bc_multiome/fig1")
dir.create("/home/groups/MohammedLab/bc_multiome/seurat_object")


####################################################
#           Fig 1 Colors                           #
###################################################

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

hist_col=c("NAT"="#c2d9ea",
"DCIS"="#cccccb",
"ILC"="#f99c1c",
"IDC"="#c79ec9")

clin_col=c("IDC ER+/PR-/HER2+"="#f37872", 
"DCIS DCIS"="#cccccb", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"NAT NA"="#c2d9ea")

#both dcis and bloom-richardson grading
grade_col=c("3"="#540b0e","2"="#9e2a2b","1"="#e09f3e","Intermediate"="#89a59e","High"="#52796f","NA"="#cccccb")

ethnicity_col=c('NOT HISPANIC OR LATINO'="blue", "HISPANIC OR LATINO"="green","UNKNOWN"="#cccccb")
race_col=c('ASIAN'="#c98ea6", "BLACK"="#7a3394","WHITE"="#c3a5cf")

assay_col=c("0"="white","1"="black","Negative"="white","Positive"="black","NA"="#cccccb")


Idents(dat)<-factor(dat$assigned_celltype,levels=c("cancer","luminal_hs","luminal_asp","basal_myoepithelial",
"adipocyte","endothelial_vascular","endothelial_lymphatic","pericyte","fibroblast",
"myeloid","bcell","plasma","tcell"))


####################################################
#           Fig 1 Sample Heatmap                  #
###################################################
#updating metadata with finalized clinical notes 260924

age=c('DCIS_01'=31, 'DCIS_02'=49, 'DCIS_03'=61, 'IDC_01'=75, 'IDC_02'=51, 'IDC_03'=74, 'IDC_04'=67, 'IDC_05'=34, 'IDC_06'=76, 'IDC_07'=44, 'IDC_08'=63, 'IDC_09'=63, 'IDC_10'=68, 'IDC_11'=37, 'IDC_12'=67, 'IDC_13'=68, 'IDC_14'=40, 'IDC_15'=43, 'IDC_16'=75, 'ILC_01'=57, 'ILC_02'=64, 'ILC_03'=71, 'ILC_04'=65, 'ILC_05'=83, 'NAT_04'=67, 'NAT_11'=37, 'NAT_14'=50)
ethnicity=c('DCIS_01'='NOT HISPANIC OR LATINO', 'DCIS_02'='NOT HISPANIC OR LATINO', 'DCIS_03'='NOT HISPANIC OR LATINO', 'IDC_01'='NOT HISPANIC OR LATINO', 'IDC_02'='NOT HISPANIC OR LATINO', 'IDC_03'='UNKNOWN', 'IDC_04'='NOT HISPANIC OR LATINO', 'IDC_05'='NOT HISPANIC OR LATINO', 'IDC_06'='NOT HISPANIC OR LATINO', 'IDC_07'='NOT HISPANIC OR LATINO', 'IDC_08'='NOT HISPANIC OR LATINO', 'IDC_09'='NOT HISPANIC OR LATINO', 'IDC_10'='NOT HISPANIC OR LATINO', 'IDC_11'='NOT HISPANIC OR LATINO', 'IDC_12'='NOT HISPANIC OR LATINO', 'IDC_13'='NOT HISPANIC OR LATINO', 'IDC_14'='NOT HISPANIC OR LATINO', 'IDC_15'='NOT HISPANIC OR LATINO', 'IDC_16'='NOT HISPANIC OR LATINO', 'ILC_01'='NOT HISPANIC OR LATINO', 'ILC_02'='NOT HISPANIC OR LATINO', 'ILC_03'='HISPANIC OR LATINO', 'ILC_04'='NOT HISPANIC OR LATINO', 'ILC_05'='NOT HISPANIC OR LATINO', 'NAT_04'='NOT HISPANIC OR LATINO', 'NAT_11'='NOT HISPANIC OR LATINO', 'NAT_14'='NOT HISPANIC OR LATINO')
race=c('DCIS_01'='ASIAN', 'DCIS_02'='WHITE', 'DCIS_03'='WHITE', 'IDC_01'='WHITE', 'IDC_02'='WHITE', 'IDC_03'='WHITE', 'IDC_04'='WHITE', 'IDC_05'='WHITE', 'IDC_06'='WHITE', 'IDC_07'='WHITE', 'IDC_08'='WHITE', 'IDC_09'='WHITE', 'IDC_10'='WHITE', 'IDC_11'='WHITE', 'IDC_12'='WHITE', 'IDC_13'='WHITE', 'IDC_14'='BLACK', 'IDC_15'='WHITE', 'IDC_16'='WHITE', 'ILC_01'='WHITE', 'ILC_02'='WHITE', 'ILC_03'='WHITE', 'ILC_04'='WHITE', 'ILC_05'='WHITE', 'NAT_04'='WHITE', 'NAT_11'='WHITE', 'NAT_14'='WHITE')
clinical_diagnosis=c( 'DCIS_01'='DCIS', 'DCIS_02'='DCIS', 'DCIS_03'='DCIS', 'IDC_01'='IDC', 'IDC_02'='IDC', 'IDC_03'='IDC', 'IDC_04'='IDC', 'IDC_05'='IDC', 'IDC_06'='IDC', 'IDC_07'='IDC', 'IDC_08'='IDC', 'IDC_09'='IDC', 'IDC_10'='IDC', 'IDC_11'='IDC', 'IDC_12'='IDC', 'IDC_13'='IDC', 'IDC_14'='IDC', 'IDC_15'='IDC', 'IDC_16'='IDC', 'ILC_01'='ILC', 'ILC_02'='ILC', 'ILC_03'='ILC', 'ILC_04'='ILC', 'ILC_05'='ILC', 'NAT_04'='NAT', 'NAT_11'='NAT', 'NAT_14'='NAT' )
er_status=c('DCIS_01'='NA', 'DCIS_02'='NA', 'DCIS_03'='NA', 'IDC_01'='Positive', 'IDC_02'='Positive', 'IDC_03'='Positive', 'IDC_04'='Positive', 'IDC_05'='Positive', 'IDC_06'='Positive', 'IDC_07'='Positive', 'IDC_08'='Positive', 'IDC_09'='Positive', 'IDC_10'='Positive', 'IDC_11'='Positive', 'IDC_12'='Positive', 'IDC_13'='Positive', 'IDC_14'='Positive', 'IDC_15'='Positive', 'IDC_16'='Positive', 'ILC_01'='Positive', 'ILC_02'='Positive', 'ILC_03'='Positive', 'ILC_04'='Positive', 'ILC_05'='Positive', 'NAT_04'='NA', 'NAT_11'='NA', 'NAT_14'='NA')
pr_status=c('DCIS_01'='NA', 'DCIS_02'='NA', 'DCIS_03'='NA', 'IDC_01'='Negative', 'IDC_02'='Negative', 'IDC_03'='Negative', 'IDC_04'='NA', 'IDC_05'='Negative', 'IDC_06'='Positive', 'IDC_07'='Positive', 'IDC_08'='Positive', 'IDC_09'='Positive', 'IDC_10'='Positive', 'IDC_11'='Positive', 'IDC_12'='NA', 'IDC_13'='Positive', 'IDC_14'='Negative', 'IDC_15'='Positive', 'IDC_16'='Negative', 'ILC_01'='Positive', 'ILC_02'='Positive', 'ILC_03'='Positive', 'ILC_04'='Negative', 'ILC_05'='Positive', 'NAT_04'='NA', 'NAT_11'='NA', 'NAT_14'='NA')
her2_status=c('DCIS_01'='NA', 'DCIS_02'='NA', 'DCIS_03'='NA', 'IDC_01'='Negative', 'IDC_02'='Negative', 'IDC_03'='Negative', 'IDC_04'='Negative', 'IDC_05'='Positive', 'IDC_06'='Negative', 'IDC_07'='Negative', 'IDC_08'='Negative', 'IDC_09'='Negative', 'IDC_10'='Negative', 'IDC_11'='Negative', 'IDC_12'='Negative', 'IDC_13'='Negative', 'IDC_14'='Negative', 'IDC_15'='Negative', 'IDC_16'='Negative', 'ILC_01'='Negative', 'ILC_02'='Negative', 'ILC_03'='Negative', 'ILC_04'='Negative', 'ILC_05'='Negative', 'NAT_04'='NA', 'NAT_11'='NA', 'NAT_14'='NA')
grade=c('DCIS_01'='High', 'DCIS_02'='Intermediate', 'DCIS_03'='NA', 'IDC_01'='3', 'IDC_02'='3', 'IDC_03'='2', 'IDC_04'='2', 'IDC_05'='3', 'IDC_06'='3', 'IDC_07'='1', 'IDC_08'='2', 'IDC_09'='2', 'IDC_10'='2', 'IDC_11'='2', 'IDC_12'='2', 'IDC_13'='2', 'IDC_14'='2', 'IDC_15'='1', 'IDC_16'='3', 'ILC_01'='2', 'ILC_02'='2', 'ILC_03'='2', 'ILC_04'='2', 'ILC_05'='2', 'NAT_04'='NA', 'NAT_11'='NA', 'NAT_14'='NA')
multiome=c('DCIS_01'='1', 'DCIS_02'='1', 'DCIS_03'='1', 'IDC_01'='1', 'IDC_02'='1', 'IDC_03'='1', 'IDC_04'='1', 'IDC_05'='1', 'IDC_06'='1', 'IDC_07'='1', 'IDC_08'='1', 'IDC_09'='1', 'IDC_10'='1', 'IDC_11'='1', 'IDC_12'='1', 'IDC_13'='1', 'IDC_14'='1', 'IDC_15'='1', 'IDC_16'='1', 'ILC_01'='1', 'ILC_02'='1', 'ILC_03'='1', 'ILC_04'='1', 'ILC_05'='1', 'NAT_04'='1', 'NAT_11'='1', 'NAT_14'='1')
plot_order=c('DCIS_01'=1, 'DCIS_02'=2, 'DCIS_03'=3, 'IDC_01'=4, 'IDC_02'=5, 'IDC_03'=6, 'IDC_04'=7, 'IDC_05'=8, 'IDC_06'=9, 'IDC_07'=10, 'IDC_08'=11, 'IDC_09'=12, 'IDC_10'=13, 'IDC_11'=14, 'IDC_12'=15, 'IDC_13'=16, 'IDC_14'=17, 'IDC_15'=18, 'IDC_16'=19, 'ILC_01'=20, 'ILC_02'=21, 'ILC_03'=22, 'ILC_04'=23, 'ILC_05'=24, 'NAT_04'=25, 'NAT_11'=26, 'NAT_14'=27)
paired_bulk_wgs=c('DCIS_01','DCIS_02', 'DCIS_03', 'IDC_01', 'IDC_02', 'IDC_03', 'IDC_04', 'IDC_06', 'IDC_07', 'IDC_08', 'IDC_09', 'IDC_10', 'IDC_11',  'IDC_13', 'IDC_14', 'IDC_15', 'IDC_16', 'ILC_02', 'ILC_03', 'ILC_04', 'ILC_05', 'NAT_11', 'NAT_14')
bulk_wgs<-setNames(nm=names(plot_order),rep("0",length(names(plot_order))))
bulk_wgs[paired_bulk_wgs]<-"1"


dat@meta.data$age<-age[dat@meta.data$sample]
dat@meta.data$ethnicity<-ethnicity[dat@meta.data$sample]
dat@meta.data$race<-race[dat@meta.data$sample]
dat@meta.data$clinical_diagnosis<-clinical_diagnosis[dat@meta.data$sample]
dat@meta.data$er_status<-er_status[dat@meta.data$sample]
dat@meta.data$pr_status<-pr_status[dat@meta.data$sample]
dat@meta.data$her2_status<-her2_status[dat@meta.data$sample]
dat@meta.data$grade<-grade[dat@meta.data$sample]
dat@meta.data$plot_order<-plot_order[dat@meta.data$sample]


met<-dat@meta.data
met<-met[!duplicated(met$sample),]
row.names(met)<-met$sample
met$age<-as.numeric(age[met$sample])
met$multiome<-as.numeric(multiome[met$sample])
met$bulk_wgs<-as.numeric(bulk_wgs[met$sample])
met$plot_order<-as.numeric(plot_order[met$sample])

sample_heatmap<-met[c("plot_order","Manuscript_Name",
                    "age","Diagnosis",
                    "Mol_Diagnosis",
                    "ethnicity",
                    "race",
                    "clinical_diagnosis",
                    "er_status",
                    "pr_status",
                    "her2_status",
                    "grade",
                    "multiome",
                    "bulk_wgs")]
row.names(sample_heatmap)<-sample_heatmap$Manuscript_Name
sample_heatmap$Diag_MolDiag<-paste(sample_heatmap$Diagnosis,sample_heatmap$Mol_Diagnosis)
sample_heatmap<-sample_heatmap[order(sample_heatmap$plot_order),]
age_col=colorRamp2(breaks=c(min(sample_heatmap$age,na.rm=T),max(sample_heatmap$age,na.rm=T)),c("#f0f0f0","#252525"))

#plot metadata
sample_heatmap<-met[c("Manuscript_Name","Diagnosis","Mol_Diagnosis",
                      "age","ethnicity","race","clinical_diagnosis",
                      "er_status","pr_status","her2_status",
                      "grade","multiome","bulk_wgs","plot_order")]

sample_heatmap<-sample_heatmap[order(sample_heatmap$plot_order),]
ha = rowAnnotation(age=sample_heatmap$age,
                  ethnicity=sample_heatmap$ethnicity,
                  race=sample_heatmap$race,
                  clinical_diagnosis_updated=sample_heatmap$clinical_diagnosis,
                  er_status=sample_heatmap$er_status,
                  pr_status=sample_heatmap$pr_status,
                  her2_status=sample_heatmap$her2_status,
                  grade_status=sample_heatmap$grade,
                  col = list(age=age_col,
                                  ethnicity=ethnicity_col,
                                  race=race_col,
                                  clinical_diagnosis_updated=hist_col,
                                  grade_status=grade_col,
                                  er_status=assay_col,
                                  pr_status=assay_col,
                                  her2_status=assay_col,
                                  histological_type =hist_col,
                                  molecular_type=clin_col))

plt<-Heatmap(sample_heatmap[c("multiome","bulk_wgs")],
 cluster_columns=F,cluster_rows=F,
 left_annotation=ha,
 col=assay_col)

pdf(file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_sample_metadata.heatmap.pdf")
print(plt)
dev.off()


####################################################
#           Fig 1 All Cell UMAPS                   #
###################################################

p1<-DimPlot(dat,group.by="seurat_clusters",reduction = "allcells.wnn.umap")
ggsave(p1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_umap_assigned_celltype.seurat_clusters.pdf",width=10,height=10,limitsize=F)
p1<-DimPlot(dat,cols=celltype_col,group.by="assigned_celltype",reduction = "allcells.wnn.umap",col=celltype_col)
ggsave(p1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_umap_assigned_celltype.celltype.pdf",width=10,height=10,limitsize=F)
p1<-DimPlot(dat,cols=hist_col,group.by="Diagnosis",reduction = "allcells.wnn.umap")
ggsave(p1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_umap_assigned_celltype.diagnosis.pdf",width=10,height=10,limitsize=F)
p1<-DimPlot(dat,cols=clin_col,group.by="Diag_MolDiag",reduction = "allcells.wnn.umap")
ggsave(p1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_umap_assigned_celltype.diag_moldiag.pdf",width=10,height=10,limitsize=F)
p1<-DimPlot(dat,group.by="sample",reduction = "allcells.wnn.umap")
ggsave(p1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_umap_assigned_celltype.sample.pdf",width=10,height=10,limitsize=F)


# #~~~~~~~rerun umap by cancer and noncancer split~~~~~~~ 260120 #
# #~~~~~~~added harmony integration for noncancer, decided against it in main figure
# dat_noncancer<-subset(dat,assigned_celltype!="cancer")
# dat_noncancer<-multimodal_cluster(dat_noncancer,harmony_integrate=TRUE)

# p1<-DimPlot(dat_noncancer,group.by="seurat_clusters",reduction = "allcells.wnn.umap")
# p2<-DimPlot(dat_noncancer,cols=celltype_col,group.by="assigned_celltype",reduction = "allcells.wnn.umap",col=celltype_col)
# p3<-DimPlot(dat_noncancer,cols=hist_col,group.by="Diagnosis",reduction = "allcells.wnn.umap")
# p4<-DimPlot(dat_noncancer,cols=clin_col,group.by="Diag_MolDiag",reduction = "allcells.wnn.umap")
# p5<-DimPlot(dat_noncancer,group.by="sample",reduction = "allcells.wnn.umap")
# ggsave(p1/p2/p3/p4/p5,file="FIG1_umap_assigned_celltype.noncancer.pdf",width=10,height=50,limitsize=F)
# dat_cancer<-subset(dat,assigned_celltype=="cancer")
# dat_cancer[["RNA"]] <- JoinLayers(dat_cancer[["RNA"]]) #rejoining layers removes any empty layers (samples without cancer cells)
# dat_cancer<-multimodal_cluster(dat_cancer)
# p1<-DimPlot(dat_cancer,group.by="seurat_clusters",reduction = "allcells.wnn.umap")
# p2<-DimPlot(dat_cancer,cols=celltype_col,group.by="assigned_celltype",reduction = "allcells.wnn.umap",col=celltype_col)
# p3<-DimPlot(dat_cancer,cols=hist_col,group.by="Diagnosis",reduction = "allcells.wnn.umap")
# p4<-DimPlot(dat_cancer,cols=clin_col,group.by="Diag_MolDiag",reduction = "allcells.wnn.umap")
# p5<-DimPlot(dat_cancer,group.by="sample",reduction = "allcells.wnn.umap")
# ggsave(p1/p2/p3/p4/p5,file="FIG1_umap_assigned_celltype.cancer.pdf",width=10,height=50,limitsize=F)

####################################################
#           Fig 1 Stacked Celltype ID             #
###################################################
#Make stacked barplot on identities per cluster
DF<-as.data.frame(dat@meta.data %>% group_by(assigned_celltype,sample) %>% tally())
DF$assigned_celltype<-factor(DF$assigned_celltype,levels=names(celltype_col))
DF$log_count<-log10(DF$n)
plt1<-ggplot(DF,aes(x=sample,fill=assigned_celltype,y=log_count))+geom_bar(position="stack",stat="identity")+theme_minimal()+scale_fill_manual(values=celltype_col)
ggsave(plt1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_allcells.assigned_celltype_barplots.pdf",width=50,limitsize=F)

#plot of cancer cells over total count (to make proportion)
DF<-as.data.frame(dat@meta.data %>% group_by(sample) %>% tally())
DF2<-as.data.frame(dat@meta.data %>% filter(assigned_celltype=="cancer") %>% group_by(sample) %>% tally())
DF$log_count_all<-log10(DF$n)
DF2$log_count_cancer<-log10(DF2$n)
row.names(DF)<-DF$sample
row.names(DF2)<-DF2$sample

DF$log_count_cancer<-0
DF[row.names(DF2),]$log_count_cancer<-DF2$log_count_cancer

plt1<-ggplot(DF,aes(x=sample))+
geom_col(aes(y=log_count_all),fill="#666666",stat="identity")+
geom_col(aes(y=log_count_cancer),fill="#c93c96",stat="identity")+
theme_minimal()

ggsave(plt1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_allcells.assigned_cancer_barplots.pdf",width=50,limitsize=F)


#plot cellcount (no epi distinction)
DF<-as.data.frame(dat@meta.data %>% group_by(sample) %>% tally())
DF$sample<-factor(DF$sample,levels=names(plot_order))
DF$log_count<-log10(DF$n)
plt1<-ggplot(DF)+geom_bar(aes(x=sample,y=log_count),stat="identity",position="dodge")+theme_minimal()
ggsave(plt1,file="/home/groups/MohammedLab/bc_multiome/fig1/FIG1_allcells.cellcount_barplots.pdf",width=50,limitsize=F)

#figure 1 also includes cnv profiles (given by TM)
saveRDS(dat,file="/home/groups/MohammedLab/bc_multiome/seurat_object/6_merged.celltyping.SeuratObject.rds")
