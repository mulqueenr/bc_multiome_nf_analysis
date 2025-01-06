```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(GeneNMF,lib="/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3")
library(Seurat)
library(Signac)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
library(optparse)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ComplexHeatmap)
library(dendextend)
library(ggdendro)
library(circlize)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="merged.clone_annot.passqc.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#cistopic=readRDS(opt$cistopic)
#titan=readRDS(opt$titan)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
dat=readRDS(opt$object_input)

hist_col=c("NAT"="#99CCFF","DCIS"="#CCCCCC","IDC"="#FF9966","ILC"="#006633")
clin_col=c("IDC ER+/PR-/HER2-"="#f9bdbd",
                      "IDC ER+/PR+/HER2-"="#fe549d",
                      "NAT NA"="#c5eae7",          
                      "DCIS DCIS"="#707b90",
                      "ILC ER+/PR+/HER2-"="#fdd503",
                      "IDC ER+/PR-/HER2+"="#328983",
                      "ILC ER+/PR-/HER2-"="#123524")
sampled_col=c("Primary"="#8A4C80","Metastasis"="#4c9173","NAT"="#99CCFF")

inmf<-function(dat,assay="RNA",ndim=20,prefix="allcells",col_in=c("black","pink","white")){
  DefaultAssay(dat)<-assay
  seu.list <- SplitObject(dat, split.by = "sample")
  geneNMF.programs <- multiNMF(seu.list, assay=assay,  k=10:25, slot="data", min.cells.per.sample = 20, min.exp = 0.01)
  geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                          max.genes=200,
                                          nMP=ndim,
                                          metric="cosine",
                                          hclust.method="ward.D2",
                                          weight.explained=0.5)
  ph <- plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0.1,1))
  png(paste0("geneNMF_",prefix,"_heatmap.png"))
  print(ph)
  dev.off()
  sim_mat<-geneNMF.metaprograms$programs.similarity

  #cluster by all marker genes
  sum_da_dend <- as.matrix(1-sim_mat) %>% as.dist %>% hclust(method="ward.D2") %>% as.dendrogram 
  k_in<-find_k(sum_da_dend,krange=5:15)
  sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
  saveRDS(sum_da_dend,file=paste0(prefix,"_metaprograms",assay,".dend.rds"))
  member_split<-cutree(sum_da_dend,k_in$k)
  saveRDS(member_split,file=paste0(prefix,"_metaprograms",assay,".members.rds"))

  #set up metadata
  program_samples<-unlist(lapply(strsplit(row.names(sim_mat),"[.]"),"[",1))
  sample_metadata<-dat@meta.data[!duplicated(dat@meta.data$sample),]
  sample_metadata$Mol_Diagnosis<-paste(sample_metadata$Diagnosis,sample_metadata$Mol_Diagnosis)
  sample_metadata$sampled_site<-unlist(lapply(strsplit(sample_metadata$sampled_site," "),"[",1))
  row.names(sample_metadata)=sample_metadata$sample
  ha = rowAnnotation(
                        histological_type=sample_metadata[program_samples,]$Diagnosis,
                        molecular_type=sample_metadata[program_samples,]$Mol_Diagnosis,
                        sampled_site=sample_metadata[program_samples,]$sampled_site,
                        col = list(
                                      histological_type =hist_col,
                                      molecular_type=clin_col,
                                      sampled_site=sampled_col))
  pdf(paste0(prefix,"_metaprograms_heatmap","_",assay,".pdf"))
  ph<-Heatmap(as.matrix(sim_mat),
  cluster_columns=sum_da_dend,
  cluster_rows=sum_da_dend,
  row_split=k_in$k,column_split=k_in$k,
  left_annotation=ha,
  show_row_names=FALSE,
  show_column_names=FALSE,
  col=colorRamp2(c(0,0.5,1),col_in))
  print(ph)
  dev.off()
  #find enrichment in c8 signatures
  top_p_C8 <- do.call("rbind",
    lapply(1:length(geneNMF.metaprograms$metaprograms.genes), 
    function(i) {
    program_name=names(geneNMF.metaprograms$metaprograms.genes)[i]
    program_genes=unlist(geneNMF.metaprograms$metaprograms.genes[program_name])
    out<-runGSEA(program_genes, universe=rownames(dat), category = "C8")
    out$program<-paste0(program_name,"_",assay)
    return(out)
    }
    ))
  pltdat<-top_p_C8 %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)
  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_c8","_",assay,".pdf"),width=20)
  #find cancer hallmark signatures
  top_p_H <- do.call("rbind",
    lapply(1:length(geneNMF.metaprograms$metaprograms.genes), 
    function(i) {
    program_name=names(geneNMF.metaprograms$metaprograms.genes)[i]
    program_genes=unlist(geneNMF.metaprograms$metaprograms.genes[program_name])
    out<-runGSEA(program_genes, universe=rownames(dat), category = "H")
    out$program<-paste0(program_name,"_",assay)
    return(out)
    }
    ))
  pltdat<-top_p_H %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)

  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_hallmark","_",assay,".pdf"),width=20)

  #find enrichment in c1 signatures (positional)
  top_p_C1 <- do.call("rbind",
    lapply(1:length(geneNMF.metaprograms$metaprograms.genes), 
    function(i) {
    program_name=names(geneNMF.metaprograms$metaprograms.genes)[i]
    program_genes=unlist(geneNMF.metaprograms$metaprograms.genes[program_name])
    out<-runGSEA(program_genes, universe=rownames(dat), category = "C1")
    out$program<-paste0(program_name,"_",assay)
    return(out)
    }
    ))
  pltdat<-top_p_C1 %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)

  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_c1","_",assay,".pdf"),width=20)
  #find enrichment in tft (transcription factor targets)
  top_p_tft <- do.call("rbind",
    lapply(1:length(geneNMF.metaprograms$metaprograms.genes), 
    function(i) {
    program_name=names(geneNMF.metaprograms$metaprograms.genes)[i]
    program_genes=unlist(geneNMF.metaprograms$metaprograms.genes[program_name])
    out<-runGSEA(program_genes, universe=rownames(dat), category = "C3",subcategory="TFT:GTRD")
    out$program<-paste0(program_name,"_",assay)
    return(out)
    }
    ))
  pltdat<-top_p_tft %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)

  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_TFT","_",assay,".pdf"),width=20)
  #Add module scores
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  dat <- AddModuleScore_UCell(dat, features = mp.genes, assay=assay, ncores=1, name = paste0("_",assay))
  #plot per cell type
  plt<-VlnPlot(dat, features=paste0(names(mp.genes),"_",assay), group.by = "assigned_celltype",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bycelltype","_",assay,".pdf"),width=20)
  #plot per diagnosis
  dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
  plt<-VlnPlot(dat, features=paste0(names(mp.genes),"_",assay), group.by = "Diag_MolDiag",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bydiagnosis","_",assay,".pdf"),width=20)

  #plot per diagnosis
  #set merge_cluster annots to just those over 50 cells, and those not diploid
  dat$merge_cluster_50min<-NA
  for(i in names(which(table(dat$merge_cluster)>50))){
    dat@meta.data[which(dat@meta.data$merge_cluster==i),]$merge_cluster_50min<-i
  }
  dat@meta.data[which(dat@meta.data$ploidy=="diploid"),]$merge_cluster_50min<-"diploid"

  plt<-VlnPlot(dat, features=paste0(names(mp.genes),"_",assay), group.by = "merge_cluster_50min",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_byclone","_",assay,".pdf"),width=20)
  return(dat)
}

inmf(dat,assay="RNA")


dat_epi<-subset(dat,cell=row.names(dat@meta.data)[dat$reclust %in% c("cancer_luminal_epithelial","luminal_epithelial","basal_epithelial")])
dat_epi<-inmf(dat_epi,assay="RNA",prefix="epi")
saveRDS(dat_epi,file="merged.nmfmodules.SeuratObject.rds")




```
