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
library(TxDb.Hsapiens.UCSC.hg38.knownGene)#local

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

cistopic_peakscores<-function(cistopic_in_file,topic_cistopic=NULL){
  cistopic_in<-readRDS(cistopic_in_file)
  sample_name<-strsplit(basename(cistopic_in_file),"[.]")[[1]][1]
  if(nchar(strsplit(sample_name,"_")[[1]][2])!=2){
        sample_name<-gsub("_","_0",sample_name)  } #uncorrect sample names for reading in file
  cistopic_in<-getRegionsScores(cistopic_in, method='NormTop', scale=TRUE)
  if(is.null(topic_cistopic)){
    cistopic_out<-cistopic_in@region.data[which(startsWith(colnames(cistopic_in@region.data),"Scores"))]
    colnames(cistopic_out)<-paste(sample_name,colnames(cistopic_out),sep="_")
  }else{
    cistopic_in<-annotateRegions(cistopic_in, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')
    cistopic_in<-cistopic_in@region.data
    cistopic_out<-cistopic_in[colnames(cistopic_in) %in% c("seqnames","start","end","width","nCounts","nCells","annotation","geneId","distanceToTSS","ENSEMBL","SYMBOL")]
    cistopic_out$topic_peak_scores<-unlist(cistopic_in[which(endsWith(colnames(cistopic_in),topic_cistopic))])
    }
    return(cistopic_out)

}


#change this to read in gene scores for each topic
report_peaks<-function(i,cistopic_objs,member_split_df){
  print("Calculating top gene scores for metaprogram",i)
  tmp=member_split_df[member_split_df$cluster==i,]
  peaks<-lapply(1:nrow(tmp), function(x){
      sample=tmp[x,"sample"]
      sample_cistopic<-gsub("_0","_",sample) #uncorrect sample names for reading in file
      topic=tmp[x,"topic"]
      topic_cistopic<-gsub("_","",topic)
      print(paste("Reading in ",sample,topic))
      cistopic_file<-cistopic_objs[grep(cistopic_objs,pattern=paste0(sample_cistopic,".cistopic"))]
      topic_peaks_score<-as.data.frame(cistopic_peakscores(cistopic_file,topic_cistopic))
      topic_peaks_score$peaks<-row.names(topic_peaks_score)
      topic_peaks_score$topic<-paste(sample,topic,sep="_")
      return(topic_peaks_score)
      })
  peaks<-do.call("rbind",peaks)
  peaks$topic_peaks_score<-as.numeric(peaks$topic_peak_score)
  peaks_out<- peaks %>% group_by(peaks) %>% dplyr::summarize(
    count=n(),nCounts=sum(nCounts),nCells=sum(nCells),
    annotation=dplyr::first(annotation),geneId=dplyr::first(geneId),
    distanceToTSS=dplyr::first(distanceToTSS),ENSEMBL=dplyr::first(ENSEMBL),SYMBOL=dplyr::first(SYMBOL),
    mean_peakscore=mean(topic_peaks_score),
    median_peakscore=median(topic_peaks_score))
  peaks_out$metaprogram_cluster<-i
  peaks_out<-peaks_out %>% arrange(desc(median_peakscore)) %>% head(n=500)
  return(peaks_out)
}
 
prerun_cistopic_correlations<-function(cistopic_objs,dat,col_in=cividis(3),prefix="epi_cistopic_cor"){
  #run correlations on pairwise gene lists
  #add variable feature filter (run on dat and then filter cistopic_ps)
  cistopic_ps<-lapply(cistopic_objs,cistopic_peakscores)
  cistopic_ps<-purrr::reduce(map(cistopic_ps, ~ as_tibble(.x, rownames = "rn")), full_join, by = "rn")
  #filter to variable peaks
  DefaultAssay(dat)<-"peaks"
  dat <- RunTFIDF(dat)
  dat <- FindVariableFeatures(dat,nfeatures=10000)
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
    col=colorRamp2(quantile(unlist(cor_topics),c(0.1,0.9,0.99),na.rm=TRUE),cividis(3)))
  print(ph)
  dev.off()

  #then make gene list per topic cluster and output (with weights (i.e. how often they occur))
  member_split_df<-data.frame(
    sample=unlist(lapply(strsplit(as.character(names(member_split)),"_Scores_"),"[",1)),
    topic=paste0("Topic_",unlist(lapply(strsplit(as.character(names(member_split)),"_Scores_Topic"),"[",2))),
    cluster=member_split)

  peaks_out<-lapply(
    unique(member_split_df$cluster),
    function(i) report_peaks(i=i,member_split_df=member_split_df,cistopic_objs=cistopic_objs))

  peaks_out<-do.call("rbind",peaks_out)
  write.csv(peaks_out,file=paste0(prefix,"metaprogram.toppeaks.csv"))

  #find cancer hallmark signatures
    top_p_H <- do.call("rbind",
    lapply(unique(peaks_out$metaprogram_cluster), 
    function(i) {
    program_name=i
    program_genes=unlist(peaks_out[peaks_out$metaprogram_cluster==program_name,]$SYMBOL)
    out<-runGSEA(program_genes, universe=row.names(dat@assays$RNA), category = "H")
    out$program<-paste0(program_name)
    return(out)
    }
    ))
  pltdat<-top_p_H %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)
  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_H","_","cistopic",".pdf"),width=20)

  #Add module scores
  peaks_out$peaks<-gsub(":","-",peaks_out$peaks)

  metaprogram_peaks<-list()
  for(i in unique(peaks_out$metaprogram_cluster)){
    metaprogram_peaks[[paste0("MetaTopic_",i)]]<-unlist(peaks_out[peaks_out$metaprogram_cluster==i,]$peaks)
  }
  dat <- AddModuleScore_UCell(dat, features = metaprogram_peaks, assay="peaks", ncores=1, name = paste0("_","cistopic_metaprogram"),maxRank=2000)

  dat[["cistopic_metaprograms"]]<-CreateAssayObject(data=t(dat@meta.data[colnames(dat@meta.data) %in% paste0(names(metaprogram_peaks),"_",assay)]))
  #plot per cell type
  assay="cistopic_metaprogram"
  plt<-VlnPlot(dat, features=paste0(names(metaprogram_peaks),"_",assay), group.by = "HBCA_predicted.id",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bycelltype","_",assay,".pdf"),width=20)

  #plot per diagnosis
  dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
  plt<-VlnPlot(dat, features=paste0(names(metaprogram_peaks),"_",assay), group.by = "Diag_MolDiag",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bydiagnosis","_",assay,".pdf"),width=20)

  Idents(dat)<-dat$Diag_MolDiag
  diff<-FindAllMarkers(dat,assay="cistopic_metaprograms",logfc.threshold=0)
  write.table(diff,file=paste0(prefix,"_metaprograms","_bydiagnosis","_",assay,".differential_cistopics.tsv"),sep="\t",col.names=T,row.names=T)
  plt<-DotPlot(dat,assay="cistopic_metaprograms",features=row.names(dat@assays$cistopic_metaprograms))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bydiagnosis","_",assay,".dotplot.pdf"))
  return(dat)
}


cistopic_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/cistopic_objects"
#cistopic_objs<-list.files(path=cistopic_path,pattern="*cistopic.cistopicObject.rds$",full.names=TRUE)
#process_cistopic_topics(cistopic_objs=cistopic_objs,met=met,out_prefix="cistopic_allcells",cistopic_path=cistopic_path)


cistopic_epi_objs<-list.files(path=cistopic_path,pattern="*cistopic_epithelial.cistopicObject.rds$",full.names=TRUE)
out<-prerun_cistopic_correlations(cistopic_objs=cistopic_epi_objs,prefix="epi_cistopic_cor",col_in=cividis(3))
saveRDS(dat,file="epi_cistopic_SeuratObject.Rds")



###### ChromVAR on DMR sites ####
library(BiocParallel)
register(MulticoreParam(5)) # Use 50 cores
library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38) #locally installed
library(TFBSTools)
library(SummarizedExperiment)
library(Signac)
library(optparse)
library(circlize)
library(viridis)

setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3")


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
opt$object_input="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/epi_cistopic_SeuratObject.Rds"
dat=readRDS(opt$object_input)
metaprogram_peaks<-read.table(file="epi_cistopic_cormetaprogram.toppeaks.csv",sep=",",header=T)

motif_ix <- dat@assays$peaks@motifs@data
motif_names<-dat@assays$peaks@motifs@motif.names

motif_counts_per_metaprogram<-function(cluster){
  peaks<-metaprogram_peaks[metaprogram_peaks$metaprogram_cluster==cluster,]$peaks
  peaks<-gsub(":","-",peaks)
  motif_counts<-colSums(motif_ix[peaks,])
  return(motif_counts)
}

hypergeo_motif<-function(motif,metaprogram){
  pval<-phyper(
  q=motif_enrichment[metaprogram,motif],
  m=unname(motif_ix_background[motif]),
  n=nrow(motif_ix)-motif_ix_background[motif],
  k=nrow(metaprogram_peaks[metaprogram_peaks$metaprogram_cluster==metaprogram,]),
  lower.tail=FALSE
  )
return(c(motif,metaprogram,pval))
}

motif_enrichment<-as.data.frame(do.call("rbind",lapply(unique(metaprogram_peaks$metaprogram_cluster),motif_counts_per_metaprogram)))
row.names(motif_enrichment)<-unique(metaprogram_peaks$metaprogram_cluster)
motif_ix_background<-colSums(motif_enrichment) #just limiting to peaks in metaprograms
hypergeo_out<-list()

for(i in colnames(motif_enrichment)){
  for (j in row.names(motif_enrichment)){
    hypergeo_out<-c(hypergeo_out,list(hypergeo_motif(i,j)))
  }
}

motif_enrichment<-as.data.frame(do.call("rbind",hypergeo_out))
colnames(motif_enrichment)<-c("motif_ix","metaprogram_cluster","pval")
motif_enrichment$padj<-p.adjust(motif_enrichment$pval, method = "bonferroni", n = length(motif_enrichment$pval))
motif_enrichment$tf_name<-motif_names[motif_enrichment$motif_ix]
motif_enrichment$log10_padj<-log10(motif_enrichment$padj)
#get list of TFs that are sig enriched at least once
tf_subset<-unique(motif_enrichment[motif_enrichment$padj<0.01,]$motif_ix)
library(ComplexHeatmap)
library(reshape2)
motif_enrichment<-motif_enrichment[motif_enrichment$motif_ix %in% tf_subset,]

out_heatmap<-dcast(motif_enrichment, motif_ix ~ metaprogram_cluster, value.var="log10_padj")
row.names(out_heatmap)<-out_heatmap$motif_ix
out_heatmap<-out_heatmap[2:ncol(out_heatmap)]

col=colorRamp2(rev(seq(0.01,1,length.out=10)),cividis(10))

out_heatmap<-log10(out_heatmap)
pdf("test_tf_enrichment.pdf",width=5,height=30)
row_labels=motif_enrichment[match(row.names(out_heatmap), unique(motif_enrichment$motif_ix)),]$tf_name
Heatmap(out_heatmap,row_labels=row_labels,row_names_gp = gpar(fontsize = 5),col=col)
dev.off()