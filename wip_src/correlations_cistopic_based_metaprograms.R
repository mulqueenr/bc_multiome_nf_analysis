```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```

library(GeneNMF,lib="/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3") #local
library(Seurat)
library(Signac)
library(ggplot2)
library(UCell) #local
library(patchwork)
library(Matrix)
library(RcppML) #local
library(viridis)
library(optparse)
library(msigdbr) #local
library(fgsea) #local
library(dplyr)
library(TITAN)
library(cisTopic)
library(ComplexHeatmap)
library(reshape2)
library(dendextend)
library(purrr)
library(circlize)

cistopic_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/cistopic_objects"
cistopic_objs<-list.files(path=cistopic_path, pattern="*.cistopic.cistopicObject.rds$",full.names=TRUE)

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

cistopic_peakscores<-function(x){
  cistopic_in<-readRDS(x)
  sample_name<-strsplit(basename(x),"[.]")[[1]][1]
  if(nchar(strsplit(sample_name,"_")[[1]][2])!=2){
        sample_name<-gsub("_","_0",sample_name) #uncorrect sample names for reading in file
  }
  cistopic_in<-getRegionsScores(cistopic_in, method='NormTop', scale=TRUE)
  cistopic_in<-cistopic_in@region.data
  cistopic_in<-cistopic_in[startsWith(colnames(cistopic_in),prefix="Scores")]
  colnames(cistopic_in)<-paste(sample_name,colnames(cistopic_in),sep="_")
  return(cistopic_in)
}


#change this to read in gene scores for each topic
report_genes<-function(i,cistopic_objs,member_split_df){
  print("Calculating top gene scores for",i)
  tmp=member_split_df[member_split_df$cluster==i,]
  genes<-lapply(1:nrow(tmp), function(x){
      sample=tmp[x,"sample"]
      sample<-gsub("_0","_",sample) #uncorrect sample names for reading in file
      topic=tmp[x,"topic"]
      print(paste("Reading in ",sample,topic))
      cistopic<-readRDS(cistopic_objs[grep(pattern=paste0(sample,"[.]"),x=cistopic_objs)])
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
 
prerun_cistopic_correlations<-function(cistopic_objs,dat,col_in=c("#CAC1B6","#FDF8C1","#A91E2C"),prefix="all_cells_cistopic_cor"){
  #run correlations on pairwise gene lists
  #add variable feature filter (run on dat and then filter cistopic_ps)
  cistopic_ps<-lapply(cistopic_objs,cistopic_peakscores)
  cistopic_ps<-purrr::reduce(map(cistopic_ps, ~ as_tibble(.x, rownames = "rn")), full_join, by = "rn")
  #filter to variable peaks
  DefaultAssay(dat)<-"peaks"
  dat <- RunTFIDF(dat)
  dat <- FindVariableFeatures(dat,nfeatures=20000)
  var_feat<-dat@assays$peaks@var.features
  cistopic_ps$rn<-gsub(cistopic_ps$rn,pattern=":",replacement="-")
  cistopic_ps<-cistopic_ps[cistopic_ps$rn %in% var_feat,]
  cistopic_ps<-cistopic_ps[,colnames(cistopic_ps)!="rn"]
  #calculate topic overlap across topics
  topic_peak_overlap<-lapply(colnames(cistopic_ps),function(x) {
    lapply(colnames(cistopic_ps), function(y) {
        sum(complete.cases(cistopic_ps[c(x,y)]))
     })})
  summary(colSums(!is.na(cistopic_ps)))
  #correlate topics
  cor_topics<-cor(cistopic_ps,use="pairwise.complete.obs",method="spearman")
  sum_da_dend <- cor_topics %>% dist() %>% hclust() %>% as.dendrogram 
  saveRDS(sum_da_dend,file=paste0(prefix,"_metaprograms.COR_cistopic.dend.rds"))
  k_in<-find_k(sum_da_dend,krange=10:20)
  sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
  member_split<-cutree(sum_da_dend,k_in$k)
  saveRDS(member_split,file=paste0(prefix,"_metaprograms.COR_cistopic.members.rds"))


  #set up metadata
  program_samples<-unlist(lapply(strsplit(row.names(cor_topics),"_Score"),"[",1))
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
  pdf(paste0(prefix,"_metaprograms_heatmap_cistopic.pdf"))
  ph<-Heatmap(as.matrix(cor_topics),
    cluster_columns=sum_da_dend,
    cluster_rows=sum_da_dend,
    row_split=k_in$k,column_split=k_in$k,
    left_annotation=ha,
    show_row_names=FALSE,
    show_column_names=FALSE,
    name = "cisTopic Topics",
    col=colorRamp2(quantile(unlist(cor_topics),c(0.01,0.5,0.99)),col_in))
  print(ph)
  dev.off()

  #then make gene list per topic cluster and output (with weights (i.e. how often they occur))
  member_split_df<-data.frame(
    sample=unlist(lapply(strsplit(as.character(names(member_split)),"_Scores_"),"[",1)),
    topic=paste0("Topic_",unlist(lapply(strsplit(as.character(names(member_split)),"_Scores_Topic"),"[",2))),
    cluster=member_split)

  genes_out<-lapply(
    unique(member_split_df$cluster),
    function(i) report_genes(i=i,member_split_df=member_split_df,cistopic_objs=cistopic_objs))

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

#add clustering and k selection for output
cistopic_top_sites<-function(x){
  cistopic<-readRDS(x)
  sample_name=strsplit(basename(x),"[.]")[[1]][1]
  cistopic <- getRegionsScores(cistopic, method='NormTop', scale=TRUE)
  cistopic <- binarizecisTopics(cistopic, thrP=0.975, plot=FALSE)
  cistopic_peaks<-lapply(names(cistopic@binarized.cisTopics), function(i) {
    tmp<-as.data.frame(head(cistopic@binarized.cisTopics[[i]],n=500))
    tmp$sample<-sample_name
    tmp$topic<-i 
    return(tmp)
    })
  return(cistopic_peaks)
}


report_peaks<-function(i,member_split_df){
  tmp=member_split_df[member_split_df$cluster==i,]
  peaks<-unlist(lapply(1:nrow(tmp), function(x){
      sample=tmp[x,"sample"]
      topic=tmp[x,"topic"]
      sample_idx=which(unlist(lapply(cistopic_topic_sites,function(k) k[[1]]$sample[1])) %in% sample)
      topic_idx=which(unlist(lapply(cistopic_topic_sites[[sample_idx]],function(k) k$topic[1]))  %in% topic)
      return(row.names(cistopic_topic_sites[[sample_idx]][[topic_idx]]))
      }))
  peaks_out<-as.data.frame(table(peaks))
  peaks_out$cluster<-i
  return(peaks_out)
}

process_cistopic_topics<-function(cistopic_objs,met,out_prefix,cistopic_path){
  cistopic_topic_sites<-lapply(cistopic_objs,function(x) cistopic_top_sites(x)) 
  cistopic_overlaps<-lapply(cistopic_topic_sites,function(i){
    lapply(cistopic_topic_sites, function(j){
      lapply(i,function(a){
        lapply(j,function(b){
          c(a$sample[1],b$sample[1],a$topic[1],b$topic[1],sum(row.names(a) %in% row.names(b)))
          })
        })
      })
    })

  df<-as.data.frame(matrix(unlist(cistopic_overlaps),ncol=5,byrow=T))
  colnames(df)<-c("dat1","dat2","topic_dat1","topic_dat2","overlap")
  df$overlap<-as.numeric(df$overlap)
  out_mat<-as.data.frame(dcast(dat1+topic_dat1~dat2+topic_dat2,data=df,value.var="overlap"))
  row.names(out_mat)<-unique(paste(out_mat$dat1,out_mat$topic_dat1))#check that row is correct for this
  out_mat<-out_mat[,3:ncol(out_mat)]
  ha = rowAnnotation(mol_diag = met[unlist(lapply(strsplit(row.names(out_mat),split=" "),"[",1)),]$full_diag,
                      diag = unlist(lapply(strsplit(row.names(out_mat),split="_"),"[",1)),
                    col = list(mol_diag = c(
                      "IDC ER+/PR-/HER2-"="#f9bdbd",
                      "IDC ER+/PR+/HER2-"="#fe549d",
                      "NAT NA"="#c5eae7",          
                      "DCIS DCIS"="#707b90",
                      "ILC ER+/PR+/HER2-"="#fdd503",
                      "IDC ER+/PR-/HER2+"="#328983",
                      "ILC ER+/PR-/HER2-"="#123524"),
                    diag = c(
                      "DCIS" = "#CCCCCC", 
                      "IDC" = "#FF6699", 
                      "ILC" = "#006633", 
                      "NAT" = "#99CCFF"))
                    )

  #cluster by all marker genes
  sum_da_dend <- out_mat %>% dist %>% hclust %>% as.dendrogram 
  k_in<-find_k(sum_da_dend,krange=5:20)
  sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
  saveRDS(sum_da_dend,file=paste0(cistopic_path,"/",out_prefix,".dend.rds"))
  member_split<-cutree(sum_da_dend,k_in$k)
  saveRDS(member_split,file=paste0(cistopic_path,"/",out_prefix,".members.rds"))
  pdf(paste0(cistopic_path,"/",out_prefix,".heatmap.pdf"))

  top_ha = columnAnnotation(topic_cluster = as.character(member_split),
   col = list(topic_cluster=setNames(
      nm=unique(as.character(member_split)), colorRampPalette(brewer.pal(8, "Accent"))(length(unique(as.character(member_split))))
    )))

  plt<-Heatmap(out_mat, name = "CisTopic Topics", 
      col=colorRamp2(c(0,250,500),c("black","green","white")),
      row_names_gp = gpar(fontsize = 0.5),
      column_names_gp = gpar(fontsize =0.5),
      width=20,height=20, right_annotation = ha,
      cluster_columns=sum_da_dend,cluster_rows=sum_da_dend)
  print(plt)
  dev.off()

  #then make gene list per topic cluster and output (with weights (i.e. how often they occur))
  member_split_df<-data.frame(
    sample=unlist(lapply(strsplit(as.character(names(member_split))," "),"[",1)),
    topic=unlist(lapply(strsplit(as.character(names(member_split))," "),"[",2)),
    cluster=member_split)

  sites_out<-lapply(
    unique(member_split_df$cluster),
    function(i) report_peaks(i=i,member_split_df=member_split_df))

  sites_out<-do.call("rbind",sites_out)
  saveRDS(sites_out,file=paste0(cistopic_path,"/",out_prefix,".peaks_per_cluster.rds"))
  sites_out_filt<-sites_out[sites_out$Freq>1,]
  write.table(sites_out_filt,sep="\t",col.names=T,row.names=F,quote=F,file=paste0(cistopic_path,"/",out_prefix,".peaks_per_cluster.filt.tsv"))
  return(sites_out)

}


cistopic_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/cistopic_objects"
cistopic_objs<-list.files(path=cistopic_path,pattern="*cistopic.cistopicObject.rds$",full.names=TRUE)
process_cistopic_topics(cistopic_objs=cistopic_objs,met=met,out_prefix="cistopic_allcells",cistopic_path=cistopic_path)


cistopic_epi_objs<-list.files(path=cistopic_path,pattern="*cistopic_epithelial.cistopicObject.rds$",full.names=TRUE)
process_cistopic_topics(cistopic_objs=cistopic_epi_objs,met=met,out_prefix="cistopic_epithelial",cistopic_path=cistopic_path)


```

```R




