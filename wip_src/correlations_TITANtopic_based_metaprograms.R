```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```

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
library(TITAN)
library(ComplexHeatmap)
library(reshape2)
library(dendextend)
library(purrr)
library(circlize)


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character"),
  make_option(c("-c", "--cistopic"), type="character", default=NULL, 
              help="List of sample cisTopic RDS files", metavar="character"),
  make_option(c("-t", "--titan"), type="character", default=NULL, 
              help="List of sample TITAN RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#cistopic=readRDS(opt$cistopic)
#titan=readRDS(opt$titan)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3")
opt$object_input="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects/merged.geneactivity.SeuratObject.rds"
dat=readRDS(opt$object_input)
met<-dat@meta.data[!duplicated(dat@meta.data$sample),]

hist_col=c("NAT"="#99CCFF","DCIS"="#CCCCCC","IDC"="#FF9966","ILC"="#006633")
clin_col=c("IDC ER+/PR-/HER2-"="#f9bdbd",
                      "IDC ER+/PR+/HER2-"="#fe549d",
                      "NAT NA"="#c5eae7",          
                      "DCIS DCIS"="#707b90",
                      "ILC ER+/PR+/HER2-"="#fdd503",
                      "IDC ER+/PR-/HER2+"="#328983",
                      "ILC ER+/PR-/HER2-"="#123524")
sampled_col=c("Primary"="#8A4C80","Metastasis"="#4c9173","NAT"="#99CCFF")

titan_genescores<-function(x){
  titan_in<-readRDS(x)
  sample_name<-strsplit(basename(x),"[.]")[[1]][1]
  if(nchar(strsplit(sample_name,"_")[[1]][2])!=2){
        sample_name<-gsub("_","_0",sample_name) #uncorrect sample names for reading in file
  }
  titan_in<-GeneScores(titan_in)
  colnames(titan_in)<-paste(sample_name,colnames(titan_in),sep="_")
  return(titan_in)
}

prerun_topic_correlations<-function(titan_objs,dat,col_in=c("#CAC1B6","#FDF8C1","#A91E2C"),prefix="all_cells_titan_cor"){
  #run correlations on pairwise gene lists
  titan_gs<-lapply(titan_objs,titan_genescores)
  titan_gs<-purrr::reduce(map(titan_gs, ~ as_tibble(.x, rownames = "rn")), full_join, by = "rn")
  titan_gs<-titan_gs[,colnames(titan_gs)!="rn"]
  topic_gene_overlap<-lapply(colnames(titan_gs),function(x) {
    lapply(colnames(titan_gs), function(y) {
        sum(complete.cases(titan_gs[c(x,y)]))
     })})
  summary(colSums(!is.na(titan_gs)))
  cor_topics<-cor(titan_gs,use="pairwise.complete.obs",method="spearman")
  sum_da_dend <- cor_topics %>% dist() %>% hclust() %>% as.dendrogram 
  saveRDS(sum_da_dend,file=paste0(prefix,"_metaprograms.COR_TITAN.dend.rds"))
  k_in<-find_k(sum_da_dend,krange=5:15)
  sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
  member_split<-cutree(sum_da_dend,k_in$k)
  saveRDS(member_split,file=paste0(prefix,"_metaprograms.COR_TITAN.members.rds"))

  #set up metadata
  program_samples<-unlist(lapply(strsplit(row.names(cor_topics),"_Topic"),"[",1))
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

  #plot
  pdf(paste0(prefix,"_metaprograms_heatmap_TITAN.pdf"))
  ph<-Heatmap(as.matrix(cor_topics),
    cluster_columns=sum_da_dend,
    cluster_rows=sum_da_dend,
    row_split=k_in$k,column_split=k_in$k,
    left_annotation=ha,
    show_row_names=FALSE,
    show_column_names=FALSE,
    name = "TITAN Topics",
    col=colorRamp2(quantile(unlist(cor_topics),c(0.01,0.5,0.99)),col_in))
  print(ph)
  dev.off()

  #then make gene list per topic cluster and output (with weights (i.e. how often they occur))
  member_split_df<-data.frame(
    sample=unlist(lapply(strsplit(as.character(names(member_split)),"_Topic_"),"[",1)),
    topic=paste0("Topic_",unlist(lapply(strsplit(as.character(names(member_split)),"_Topic_"),"[",2))),
    cluster=member_split)

  genes_out<-lapply(
    unique(member_split_df$cluster),
    function(i) report_genes(i=i,member_split_df=member_split_df,titan_objs=titan_objs))

  genes_out<-do.call("rbind",genes_out)
  write.csv(genes_out,file=paste0(prefix,"metaprogram.topgenes.csv"))

  #find enrichment in c8 signatures
  top_p_C8 <- do.call("rbind",
    lapply(unique(genes_out$metaprogram_cluster), 
    function(i) {
    program_name=i
    program_genes=unlist(genes_out[genes_out$metaprogram_cluster==program_name,]$genes)
    out<-runGSEA(program_genes, universe=row.names(dat@assays$RNA), category = "C8")
    out$program<-paste0(program_name)
    return(out)
    }
    ))
  pltdat<-top_p_C8 %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)

  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_c8","_","TITAN",".pdf"),width=20)

  #find cancer hallmark signatures
    top_p_H <- do.call("rbind",
    lapply(unique(genes_out$metaprogram_cluster), 
    function(i) {
    program_name=i
    program_genes=unlist(genes_out[genes_out$metaprogram_cluster==program_name,]$genes)
    out<-runGSEA(program_genes, universe=row.names(dat@assays$RNA), category = "H")
    out$program<-paste0(program_name)
    return(out)
    }
    ))
  pltdat<-top_p_H %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)

  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_H","_","TITAN",".pdf"),width=20)

  #Add module scores
  metaprogram_genes<-list()
  for(i in unique(genes_out$metaprogram_cluster)){
    metaprogram_genes[[paste0("MetaTopic_",i)]]<-unlist(genes_out[genes_out$metaprogram_cluster==i,]$genes)
  }
  dat <- AddModuleScore_UCell(dat, features = metaprogram_genes, assay="SoupXRNA", ncores=1, name = paste0("_","TITAN_metaprogram"))

  #plot per cell type
  assay="TITAN_metaprogram"
  plt<-VlnPlot(dat, features=paste0(names(metaprogram_genes),"_",assay), group.by = "HBCA_predicted.id",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bycelltype","_",assay,".pdf"),width=20)

  #plot per diagnosis
  dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
  plt<-VlnPlot(dat, features=paste0(names(metaprogram_genes),"_",assay), group.by = "Diag_MolDiag",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bydiagnosis","_",assay,".pdf"),width=20)
  return(dat)
}

#change this to read in gene scores for each topic
report_genes<-function(i,titan_objs,member_split_df){
  print("Calculating top gene scores for",i)
  tmp=member_split_df[member_split_df$cluster==i,]
  genes<-lapply(1:nrow(tmp), function(x){
      sample=tmp[x,"sample"]
      sample<-gsub("_0","_",sample) #uncorrect sample names for reading in file
      topic=tmp[x,"topic"]
      print(paste("Reading in ",sample,topic))
      titan<-readRDS(titan_objs[grep(pattern=paste0(sample,"[.]"),x=titan_objs)])
      topic_genes_gene_score<-as.data.frame(GeneScores(titan))
      colnames(topic_genes_gene_score)<-paste(sample,"Topic",1:ncol(topic_genes_gene_score),sep="_")
      topic_genes_gene_score<-as.data.frame(topic_genes_gene_score[paste(sample,topic,sep="_")])
      colnames(topic_genes_gene_score)<-c("gene_score")
      topic_genes_gene_score$genes<-row.names(topic_genes_gene_score)
      topic_genes_gene_score$topic<-paste(sample,topic,sep="_")
      return(topic_genes_gene_score)
      })
  genes<-do.call("rbind",genes)
    genes$gene_score<-as.numeric(genes$gene_score)
  genes_out<- genes %>% group_by(genes) %>% dplyr::summarize(count=n(),mean_genescore=mean(gene_score),median_genescore=median(gene_score))
  genes_out$metaprogram_cluster<-i
  genes_out<-genes_out %>% arrange(desc(mean_genescore)) %>% head(n=100)
  return(genes_out)
}
 

titan_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/titan_objects"


#titan all cells
titan_objs<-list.files(path=titan_path, pattern="*titan.titanObject.rds$",full.names=TRUE)
titan_cluster_genes<-prerun_topic_correlations(dat=dat,
                                                titan_objs=titan_objs,
                                                prefix="allcells_cor")

#titan epithelial
dat<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
met<-dat@meta.data[!duplicated(dat@meta.data$sample),]

titan_objs<-list.files(path=titan_path, pattern="*titan_epithelial.titanObject.rds$",full.names=TRUE)
titan_cluster_genes<-prerun_topic_correlations(dat=dat,
                                                titan_objs=titan_objs,
                                                prefix="epithelial_cor")

