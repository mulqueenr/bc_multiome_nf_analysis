```bash
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(GeneNMF)
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

option_list = list(
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

inmf<-function(dat,assay="SCT",ndim=15,prefix="allcells"){
  DefaultAssay(dat)<-assay
  seu.list <- SplitObject(dat, split.by = "sample")
  geneNMF.programs <- multiNMF(seu.list, assay=assay,  k=10:25, slot="data", min.cells.per.sample = 100)
  geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                          nprograms=20,
                                          max.genes=100,
                                          hclust.method="ward.D2",
                                          min.confidence=0.3)
  t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, function(x) head(x,n=10))))

  ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
  ggsave(ph,file=paste0(prefix,"_metaprograms",assay,".pdf"))

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

  #Add module scores
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  dat <- AddModuleScore_UCell(dat, features = mp.genes, assay=assay, ncores=1, name = paste0("_",assay))

  #plot per cell type
  plt<-VlnPlot(dat, features=paste0(names(mp.genes),"_",assay), group.by = "HBCA_predicted.id",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bycelltype","_",assay,".pdf"),width=20)

  #plot per diagnosis
  dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
  plt<-VlnPlot(dat, features=paste0(names(mp.genes),"_",assay), group.by = "Diag_MolDiag",pt.size = 0, stack=TRUE)
  ggsave(plt,file=paste0(prefix,"_metaprograms","_bydiagnosis","_",assay,".pdf"),width=20)
  return(dat)
}

inmf(dat,assay="SCT")
inmf(dat,assay="GeneActivity")
#inmf(dat,assay="chromvar")

dat<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
dat<-inmf(dat,assay="SCT",prefix="epi")
dat<-inmf(dat,assay="GeneActivity",prefix="epi")
#dat<-inmf(dat,assay="chromvar",prefix="epi")

```

Run on cisTopic and TITAN methods
```R
library(GeneNMF)
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
library(cisTopic)


titan_top_genes<-function(x){
  titan<-readRDS(x)
  sample_name=strsplit(basename(x),"[.]")[[1]][1]
  TitanTopicGenes <- TopTopicGenes(titan, ngenes = 100)
  return(TitanTopicGenes)
}

compare_topics<-function(x_name,y_name){
  x=titan_topic_genes[x_name]
  y=titan_topic_genes[y_name]
  lapply(colnames(titan_topic_genes[[x_name]]),function(xi){
    lapply(colnames(titan_topic_genes[[y_name]]),function(yj){
      c(x_name,y_name,xi,yj,length(intersect(unlist(titan_topic_genes[[x_name]][,xi]),unlist(titan_topic_genes[[y_name]][,yj]))))
      })
      }) 
}

report_genes<-function(i,member_split_df){
  tmp=member_split_df[member_split_df$cluster==i,]
  genes<-unlist(lapply(1:nrow(tmp), function(x){
      sample=tmp[x,"sample"]
      topic=tmp[x,"topic"]
      titan_topic_genes[[sample]][,topic]
      }))
  genes_out<-as.data.frame(table(genes))
  genes_out$cluster<-i
  return(genes_out)
}
 

process_titan_topics<-function(titan_objs,met,out_prefix,titan_path){
  titan_topic_genes<-lapply(titan_objs,function(x) titan_top_genes(x))
  names(titan_topic_genes)<-unlist(lapply(titan_objs,function(x) strsplit(basename(x),"[.]")[[1]][1]))
  out<-lapply(names(titan_topic_genes),function(x){lapply(names(titan_topic_genes),function(y){compare_topics(x,y)})})
  df<-as.data.frame(matrix(unlist(out),ncol=5,byrow=T))
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
  saveRDS(sum_da_dend,file=paste0(titan_path,"/",out_prefix,".dend.rds"))
  member_split<-cutree(sum_da_dend,k_in$k)
  saveRDS(member_split,file=paste0(titan_path,"/",out_prefix,".members.rds"))

  top_ha = columnAnnotation(topic_cluster = as.character(member_split),
   col = list(topic_cluster=setNames(
      nm=unique(as.character(member_split)), colorRampPalette(brewer.pal(8, "Accent"))(length(unique(as.character(member_split))))
    )))

  pdf(paste0(titan_path,"/",out_prefix,".heatmap.pdf"))
  plt<-Heatmap(out_mat, name = "TITAN Topics", 
      col=colorRamp2(c(0,25,50),c("black","purple","white")),
      row_names_gp = gpar(fontsize = 0.5),
      column_names_gp = gpar(fontsize =0.5),
      width=20,height=20, right_annotation = ha,
      cluster_columns=sum_da_dend,cluster_rows=sum_da_dend,top_annotation=top_ha)
  print(plt)
  dev.off()

  #then make gene list per topic cluster and output (with weights (i.e. how often they occur))
  member_split_df<-data.frame(
    sample=unlist(lapply(strsplit(as.character(names(member_split))," "),"[",1)),
    topic=unlist(lapply(strsplit(as.character(names(member_split))," "),"[",2)),
    cluster=member_split)

  genes_out<-lapply(
    unique(member_split_df$cluster),
    function(i) report_genes(i=i,member_split_df=member_split_df))

  genes_out<-do.call("rbind",genes_out)
  saveRDS(genes_out,file=paste0(titan_path,"/",out_prefix,".genes_per_cluster.rds"))
  genes_out_filt<-genes_out[genes_out$Freq>1,]
  write.table(genes_out_filt,sep="\t",col.names=T,row.names=F,quote=F,file=paste0(titan_path,"/",out_prefix,".genes_per_cluster.filt.tsv"))
 
  return(genes_out)
}

titan_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/titan_objects"

#titan all cells
titan_objs<-list.files(path=titan_path, pattern="*titan.titanObject.rds$",full.names=TRUE)
titan_cluster_genes<-process_titan_topics(titan_objs=titan_objs,met=met,out_prefix="titan_allcells",titan_path=titan_path)

#titan epithelial
titan_objs<-list.files(path=titan_path, pattern="*titan_epithelial.titanObject.rds$",full.names=TRUE)
titan_epi_cluster_genes<-process_titan_topics(titan_objs=titan_objs,met=met,out_prefix="titan_epithelial",titan_path=titan_path)


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



met<-dat@meta.data[!duplicated(dat@meta.data$sample),]
row.names(met)<-met$sample
met$full_diag<-paste(met$Diagnosis,met$Mol_Diagnosis)

