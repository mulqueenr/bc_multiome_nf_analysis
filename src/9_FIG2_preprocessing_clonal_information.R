sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell \
--bind /home/groups/CEDAR/mulqueen/bc_multiome \
--bind /home/groups/CEDAR/scATACcnv/Hisham_data \
$sif
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(dplyr)
library(optparse)
library(parallel)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(grid)
library(ggrepel)
library(seriation)
library(org.Hs.eg.db)
library(dendextend)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GeneNMF)
library(rGREAT)
library(msigdbr,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(fgsea,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(presto,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character"),
  make_option(c("-r", "--ref_object"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri/nakshatri_multiome.geneactivity.rds", 
              help="Nakshatri reference object for epithelial comparisons", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)


clin_col=c(
"DCIS DCIS"="#cccccb", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"IDC ER+/PR-/HER2+"="#f37872", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"NAT NA"="#c2d9ea")

############# Add clone and filtering info to object
read_clone_annot<-function(i){
clone_annot<-list.files(i,pattern="_ATAC_anno.csv",full.names=T)
  if(length(clone_annot)!= 0 ){
    sample_name<-gsub(basename(clone_annot),pattern="_ATAC_anno.csv",replace="")
    clones<-read.csv(clone_annot,row.names=1)
    row.names(clones)<-gsub(row.names(clones),pattern="[.]",replace="-")
    row.names(clones)<-paste(sample_name,row.names(clones),sep="_")
    clones$merge_cluster<-paste(sample_name,clones$merge_cluster,sep="_")
    clones$atac_cluster<-paste(sample_name,clones$atac_cluster,sep="_")
    clones$rna_cluster<-paste(sample_name,clones$rna_cluster,sep="_")
    return(clones)
  }
}

############# Add clone and filtering info to object
read_cnv_plots<-function(i){
  sample<-gsub(basename(i),pattern="_consensus_cnvs.bed",replacement="")
  cnv<-read.table(i)
  cnv<-cnv[5:ncol(cnv)]
  colnames(cnv)<-gsub(colnames(cnv),pattern="[.]",replacement="-")
  colnames(cnv)<-paste(sample,colnames(cnv),sep="_")
  return(cnv)
}

clone_annot_list<-list.files("/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/",pattern="*clone_val",full.name=T)
clone_annot<-do.call("rbind",lapply(clone_annot_list,read_clone_annot))
colnames(clone_annot)<-c("merged_assay_clones","atac_clones","rna_clones","clones_celltype","clones_log2reads")
dat<-AddMetaData(dat,clone_annot)

cnv_list<-list.files("/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/",recursive=T,pattern="*_consensus_cnvs.bed",full.name=T)
cnv_list<-cnv_list[!grepl(cnv_list,pattern="unfiltered")]
cnv_windows<-read.table(cnv_list[1])
cnv_windows<-cnv_windows[1:4]
colnames(cnv_windows)<-c("chr","start","end","win")

cnv_out<-do.call("cbind",lapply(cnv_list,read_cnv_plots))
row.names(cnv_out)<-paste(cnv_windows$chr,cnv_windows$start,cnv_windows$end,sep="-")
cnv_out<-cnv_out[colnames(cnv_out) %in% row.names(dat@meta.data)]
cnv_assay <- CreateAssayObject(counts = cnv_out)
dat[["cnv"]]<-cnv_assay

#normal clones based on lack of CNV calls
normal_clones<-c(
"DCIS_03_1",
"IDC_01_4",
"IDC_02_1",
"IDC_05_1",
"IDC_06_3",
"IDC_08_2",
"IDC_09_2",
"IDC_10_3",
"IDC_11_2",
"IDC_12_4",
"ILC_02_1",
"ILC_04_2")

contamination_clones<-c(
"ILC_04_3",
"ILC_04_4"
)

dat@meta.data[dat@meta.data$merged_assay_clones %in% normal_clones,]$merged_assay_clones<-"normal"
dat@meta.data[dat@meta.data$merged_assay_clones %in% contamination_clones,]$merged_assay_clones<-"contamination"
dat@meta.data[dat@meta.data$merged_assay_clones %in% c("normal"),]$merged_assay_clones <-paste(dat@meta.data[dat@meta.data$merged_assay_clones %in% c("normal"),]$sample,"normal",sep="_")

saveRDS(dat,file="8_merged.cnv_clones.SeuratObject.rds")

dat<-readRDS(file="8_merged.cnv_clones.SeuratObject.rds")
#####Plot of heatmap of all clones ####
dat<-subset(dat,merged_assay_clones != "contamination")

cnv_col<-c("0"="#002C3E", "0.5"="#78BCC4", "1"="#F7F8F3", "1.5"="#F7444E", "2"="#aa1407", "3"="#440803")
#from Curtis et al.

dat_sub<-subset(dat,assigned_celltype=="cancer")
table(dat_sub$Diag_MolDiag,dat_sub$sample)

windows<-data.frame(chr=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",1)),
                    start=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",2)),
                    end=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",3)))
          
windows<-makeGRangesFromDataFrame(windows)

#relevant CNV genes from curtis work
#from https://www.nature.com/articles/s41416-024-02804-6#Sec20
#change RAB7L1 to RAB29
#lost RAB7L1

cnv_genes<-c('ESR1','PGR','DLEU2L', 'TRIM46', 'FASLG', 'KDM5B', 'RAB7L1', 'PFN2', 'PIK3CA', 'EREG', 'AIM1', 'EGFR', 'ZNF703', 'MYC', 'SEPHS1', 'ZMIZ1', 'EHF', 'POLD4', 'CCND1', 'P2RY2', 'NDUFC2-KCTD14', 'FOXM1', 'MDM2', 'STOML3', 'NEMF', 'IGF1R', 'TP53I13', 'ERBB2', 'SGCA', 'RPS6KB1', 'BIRC5', 'NOTCH3', 'CCNE1', 'RCN3', 'SEMG1', 'ZNF217', 'TPD52L2', 'PCNT', 'CDKN2AIP', 'LZTS1', 'PPP2R2A', 'CDKN2A', 'PTEN', 'RB1', 'CAPN3', 'CDH1', 'MAP2K4', 'GJC2', 'TERT', 'RAD21', 'ST3GAL1', 'SOCS1')
cnv_genes_class<-c('amp','amp','amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'amp', 'amp', 'amp', 'amp', 'amp')
cnv_genes<-setNames(cnv_genes_class,cnv_genes)
cnv_genes<-cnv_genes[names(cnv_genes) %in% dat@assays$ATAC@annotation$gene_name]
cnv_genes_windows<-dat@assays$ATAC@annotation[dat@assays$ATAC@annotation$gene_name %in% names(cnv_genes),] #filter annotation to genes we want
cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
windows<-findOverlaps(windows,cnv_genes_windows)

#filter to genes that actually show changes

annot<-data.frame(
  window_loc=queryHits(windows),
  gene=cnv_genes_windows$gene_name,
  cnv_class=unname(cnv_genes[cnv_genes_windows$gene_name]))

annot$col<-ifelse(annot$cnv_class=="amp","red","blue")

cell_cnv<-t(as.data.frame(dat@assays$cnv@counts))

#filter out expected amp and del that has less than 50 cells
annot_filt<-unlist(lapply(1:nrow(annot),function(i){
  annot_expected<-annot$cnv_class[i]
  if(annot_expected=="amp"){
    if(sum(cell_cnv[,annot$window_loc[i]]>1)<30){
      i<-NA
    }
  }
  if(annot_expected=="del"){
    if(sum(cell_cnv[,annot$window_loc[i]]<1)<30){
    i<-NA
    }
  }
  return(i)
}))

annot_filt<-annot_filt[!is.na(annot_filt)]
annot<-annot[annot_filt,]

hc = columnAnnotation(common_cnv = anno_mark(at = annot$window_loc, 
                        labels = annot$gene,
                        which="column",side="bottom",
                        labels_gp=gpar(col=annot$col)))


#### CNV call output folder
system(paste0("mkdir -p ",paste0(dirname(getwd()),"/cnv_calls"))) #paste cnv_call folder one dir up from wd
output_directory=paste0(dirname(getwd()),"/cnv_calls")

row_title_color<-unique(dat@meta.data[row.names(dat@meta.data) %in% colnames(dat@assays$cnv@counts),]$merged_assay_clones)
title_color<-data.frame(Diag_MolDiag=dat@meta.data[!duplicated(dat@meta.data$merged_assay_clones),]$Diag_MolDiag,merged_assay_clones=dat@meta.data[!duplicated(dat@meta.data$merged_assay_clones),]$merged_assay_clones)

title_color$col<-clin_col[title_color$Diag_MolDiag]

pdf(paste0(output_directory,"/","all_samples.cnv.heatmap.pdf"),height=90,width=40)
Heatmap(cell_cnv,
  col=cnv_col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = FALSE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  row_title_gp = gpar(col = title_color$col),
  bottom_annotation=hc,
  row_split=dat@meta.data[row.names(dat@meta.data) %in% colnames(dat@assays$cnv@counts),]$merged_assay_clones,
  column_split=factor(unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",1)),levels=paste0("chr",1:22)),
  border = TRUE)
dev.off()
print(paste0(output_directory,"/","all_samples.cnv.heatmap.pdf"))


#ADD CNV TO HEATMAP
plot_top_tf_markers(x=dat,
                    group_by="merged_assay_clones",
                    plot_by="sample_diag",
                    prefix="pairwise_by_diagnosis",
                    n_markers=20,
                    order_by_idents=TRUE,
                    outdir=paste0(output_directory,"/pairwise_by_diagnosis"))




#clonal level analyses
epicell_filter<-names(which(table(subset(dat,assigned_celltype %in% c("cancer","basal_myoepithelial","luminal_asp","luminal_hs"))$sample)>=30))
cancercell_filter<-names(which(table(subset(dat,assigned_celltype=="cancer")$sample)>=30))



#barplot of S/G2M cells
#per clone, 50 cancer cells minimum
Phase_col=c(
  "G1"="#e5f5e0",
  "S"="#a1d99b",
  "G2M"="#00441b"
)


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

#barplot of scsubtype,  
scsubtype_col=c(
  "SC_Subtype_Basal_SC"="#da3932",
  "SC_Subtype_Her2E_SC"="#f0c2cb",
  "SC_Subtype_LumA_SC"="#2b2c76",
  "SC_Subtype_LumB_SC"="#86cada"
)
#use only clones with at least 30 cell 
clone_filter<-names(which(table(dat$merged_assay_clones)>=30))

cellcycle_freq<-dat@meta.data %>% 
  filter(assigned_celltype=="cancer") %>%
  filter(sample %in% cancercell_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,Phase,merged_assay_clones) %>% 
  count(Phase,.drop=TRUE)


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
  filter(assigned_celltype=="cancer") %>%
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

#check esr1 expression on clones
dat<-SCTransform(dat)

gene_expression<-FetchData(object = dat, vars = c("ESR1","KRT5"), assay="SCT")
chromatin_expression<-FetchData(object = dat, vars = c("MA0112.3"), assay="chromvar")
dat_tmp<-AddMetaData(dat,gene_expression)
dat_tmp<-AddMetaData(dat_tmp,chromatin_expression,col.name="ESR1_TF")

esr1_clones<-dat_tmp@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) 


plt1<-ggplot(esr1_clones, aes(y=ESR1, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")


plt2<-ggplot(esr1_clones, aes(y=KRT5, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

plt3<-ggplot(esr1_clones, aes(y=ESR1_TF, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +theme(axis.text.x = element_text(angle = 90))+

  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")
ggsave(plt1/plt2/plt3,file="ESR1_perclone_boxplot.pdf",width=10,height=10)


#per clone, 30 cells per clone minimum
scsubtype_freq<-dat@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) %>% 
  count(scsubtype,.drop=TRUE)


plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")+ theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="scsubtype_perclone_percbarplot.pdf",width=10,height=10)
