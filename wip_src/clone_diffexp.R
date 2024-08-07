```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(Signac)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(grid)
library(dplyr) 
library(ggplot2)
library(ggrepel)
library(patchwork)
library(presto)
library(seriation)
library(org.Hs.eg.db)
library(ggtern)
library(optparse)


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

hist_col=c("NAT"="#99CCFF","DCIS"="#CCCCCC","IDC"="#FF9966","ILC"="#006633")
clin_col=c("IDC ER+/PR-/HER2-"="#f9bdbd",
                      "IDC ER+/PR+/HER2-"="#fe549d",
                      "NAT NA"="#c5eae7",          
                      "DCIS DCIS"="#707b90",
                      "ILC ER+/PR+/HER2-"="#fdd503",
                      "IDC ER+/PR-/HER2+"="#328983",
                      "ILC ER+/PR-/HER2-"="#123524")
sampled_col=c("Primary"="#8A4C80","Metastasis"="#4c9173","NAT"="#99CCFF")


#Grab top overlapping TFs
topTFs <- function(markers_list,celltype, padj.cutoff = 1e-2,rna=NA,ga=NA,motifs=NA) {
  ctmarkers_rna <- dplyr::filter(rna, SoupXRNA.group == celltype) %>% 
    arrange(-SoupXRNA.auc)

    if(is.data.frame(motifs)) {
    ctmarkers_motif <- dplyr::filter(motifs, chromvar.group == celltype) %>% 
      arrange(-chromvar.auc)
    }

    if(is.data.frame(ga)) {
    ctmarkers_ga<- dplyr::filter(ga, GeneActivity.group == celltype) %>% 
      arrange(-GeneActivity.auc)
    }

    if(is.data.frame(motifs) && is.data.frame(ga)){    
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
      top_tfs <- inner_join(
        x = top_tfs ,
        y = ctmarkers_ga [,c(2, 11, 6, 7)], by = "gene"
      )
    }else if(is.data.frame(motifs)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
    } else if (is.data.frame(ga)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_ga[,c(2, 11, 6, 7)], by = "gene"
      )
    } 
  auc_colnames<-grep(".auc$",colnames(top_tfs))
  top_tfs$avg_auc <-  rowMeans(top_tfs[auc_colnames])
  top_tfs <- arrange(top_tfs, -avg_auc)
  top_tfs$celltype<-celltype
  return(top_tfs)
}


#Identify top markers
Identify_Marker_TFs<-function(x,group_by.="predicted.id",assay.="SoupXRNA",prefix.,pval_filt=1){
    #x[[assay.]]<-as(object = x[[assay.]], Class = "Assay")
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by., 
      groups_use=unname(unlist(unique(x@meta.data[group_by.]))),
      y=unname(unlist(unique(x@meta.data[group_by.]))), 
      assay = 'data', seurat_assay = assay.)
    markers<-markers[markers$padj<=pval_filt,]
    write.table(markers,file=paste0(prefix.,"_",assay.,"_DE_table.tsv"),sep="\t",row.names=F,col.names=T,quote=F)
    colnames(markers) <- paste(assay., colnames(markers),sep=".")
    if (assay. == "chromvar") {
      motif.names <- markers[,paste0(assay.,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="peaks")
    } else {
    markers$gene <- markers[,paste0(assay.,".feature")]
    }
    return(markers) 
}

#Average markers across groups
average_features<-function(x=hg38_atac,features=da_tf_markers$motif.feature,assay="chromvar",group_by.="predicted.id"){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    dat_motif<-x[[assay]]@data[features,]
    dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
    sum_motif<-split(dat_motif,x@meta.data[,group_by.]) #group by rows to seurat clusters
    sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean,na.rm=T)) #take average across group
    sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame

    sum_motif<-t(scale(sum_motif))
    sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    sum_motif<-sum_motif[complete.cases(sum_motif),]
    return(sum_motif)
}

#Make a heatmap of aligned multiple modalities
plot_top_TFs<-function(x=stromal,tf_markers=da_tf_markers,prefix="stromal",group_by.="predicted.id",CHROMVAR=TRUE,GA=TRUE,height.){
    tf_rna<-average_features(x=x,features=tf_markers$gene,assay="SoupXRNA",group_by.=group_by.)
    tf_rna<-tf_rna[row.names(tf_rna) %in% tf_markers$gene,]

    tf_motif<-average_features(x=x,features=tf_markers$chromvar.feature,assay="chromvar",group_by.=group_by.)
    tf_motif<-tf_motif[row.names(tf_motif) %in% tf_markers$chromvar.feature,]
    row.names(tf_motif)<-tf_markers[tf_markers$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]

    tf_ga<-average_features(x=x,features=tf_markers$gene,assay="GeneActivity",group_by.=group_by.)
    tf_ga<-tf_ga[row.names(tf_ga) %in% tf_markers$gene,]
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_ga<-tf_ga[markers_list,]

    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]

    #set up heatmap seriation and order by RNA
    o = seriate(max(tf_rna) - tf_rna, method = "BEA_TSP")
    saveRDS(o,file=paste0(prefix,".geneactivity.dend.rds")) 
    side_ha_rna<-data.frame(ga_motif=tf_markers[get_order(o,1),]$SoupXRNA.auc)
    colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.80,0.95)),plasma(3))

    side_ha_motif<-data.frame(chromvar_motif=tf_markers[get_order(o,1),]$chromvar.auc)
    colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.80,0.95)),cividis(3))
    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-tf_markers[tf_markers$gene %in% row.names(tf_motif),]$chromvar.feature
    plt<-MotifPlot(object = x,assay="peaks",motifs = motif_list[get_order(o,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

    side_ha_ga<-data.frame(ga_auc=tf_markers[get_order(o,1),]$GeneActivity.auc)
    colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.80,0.95)),magma(3))


    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), labels =row.names(tf_rna),labels_gp=gpar(fontsize=6)))


    rna_auc<-Heatmap(side_ha_rna,
        row_order = get_order(o,1),
        col=side_ha_col,
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=7))

    rna_plot<-Heatmap(tf_rna,
        row_order = get_order(o,1),
        column_order = get_order(o,2),
        name="SoupX",
        column_title="SoupX",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 8),
        show_row_names=FALSE,
        column_names_rot=90)

      ga_auc<-Heatmap(side_ha_ga,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      ga_plot<-Heatmap(tf_ga,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90)

      motif_auc<-Heatmap(side_ha_motif,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_row_names=FALSE,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      motif_plot<-Heatmap(tf_motif,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          #top_annotation=top_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
          right_annotation=gene_ha)
      
      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot,row_title=prefix)


    pdf(paste0(prefix,".tf.heatmap.pdf"),height=height.)
    print(plt1)
    dev.off()
}

#Final wrapper function
run_top_TFs<-function(obj=stromal,prefix="stromal",i="predicted.id",n_markers=5,plot_height=10){
  markers<-lapply(c("SoupXRNA","GeneActivity","chromvar"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay,prefix.=prefix))
  names(markers)<-c("SoupXRNA","GeneActivity","chromvar")
  markers_out<-do.call("rbind",lapply(unique(obj@meta.data[,i]),function(x) head(topTFs(markers_list=markers,celltype=x,rna=markers$SoupXRNA,ga=markers$GeneActivity,motifs=markers$chromvar),n=n_markers))) #grab top 5 TF markers per celltype
  dim(markers_out)
  markers_out<-markers_out[!duplicated(markers_out$gene),]
  dim(markers_out)
  saveRDS(markers_out,file=paste0(prefix,"_celltype_TF_markers.RDS"))
  da_tf_markers<-readRDS(paste0(prefix,"_celltype_TF_markers.RDS"))
  plot_top_TFs(x=obj,tf_markers=da_tf_markers,prefix=prefix,group_by.=i)
 
}

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
    get_average <- function(v) sum(v * row.w)/sum(row.w)
    average <- apply(x, 2, get_average)
    sweep(x, 2, average)
  }

run_scsubtype_perclone<-function(obj){
  #Read in the single cell RDS object as 'Mydata'
  DefaultAssay(obj)<-"RNA"
  obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")
  Mydata<-AggregateExpression(obj,assay="RNA",return.seurat=T,features=temp_allgenes,group.by="cluster")
  Mydata <- ScaleData(Mydata, features=temp_allgenes,assay="RNA") #running only on aneuploid epithelial cells
  Mydata[["RNA"]] <- as(object = Mydata[["RNA"]], Class = "Assay")
  tocalc<-as.data.frame(Mydata@assays$RNA@scale.data) #using RNA (not SoupXRNA corrected read counts)

  #calculate mean scsubtype scores
  outdat <- matrix(0,
                  nrow=ncol(sigdat),
                  ncol=ncol(tocalc),
                  dimnames=list(colnames(sigdat),
                                colnames(tocalc)))
  for(i in 1:ncol(sigdat)){
    row <- as.character(sigdat[,i])
    row<-unique(row[row != ""])
    genes<-which(rownames(tocalc) %in% row)
    temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
    outdat[i,]<-as.numeric(temp)
  }

  #final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
  final<-outdat
  final<-as.data.frame(final)
  is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
  finalm<-as.matrix(final)

  ##Obtaining the highest call
  finalmt<-as.data.frame(t(finalm))
  finalm.sweep.t<-center_sweep(finalmt)
  Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
  finalm.sweep.t$SCSubtypeCall <- Finalnames
  finalm.sweep.t$sample <- paste(unlist(lapply(strsplit(row.names(finalm.sweep.t),"_"),"[",1)),unlist(lapply(strsplit(row.names(finalm.sweep.t),"_"),"[",2)),sep="_")
  finalm.sweep.t$sample<-unlist(lapply(finalm.sweep.t$sample, function(x) gsub("-","_",x)))
  plt<-ggtern(finalm.sweep.t,aes(x=LumA_SC,z=LumB_SC,y=Basal_SC,color=SCSubtypeCall))+geom_point()+ggtitle(obj$sample[1])
  ggsave(plt,file=paste0("./clonal_de_tfs/","clonal_",sample_in,"_SCSubtype_tern.pdf"))
return(c(finalm.sweep.t,Finalnames))
}

plot_ga_clones_covplot<-function(obj=dat_epi,samp=sample_in){
  ga_markers<-Identify_Marker_TFs(obj,group_by.="cluster",assay.="GeneActivity",prefix.=paste0("./clonal_de_tfs/","clonal_",samp))
  ga_markers<-ga_markers %>% group_by(GeneActivity.group) %>% filter(GeneActivity.logFC>0.2) %>% arrange(desc(GeneActivity.auc)) %>% slice_head(n=5)
  ga_regions<-unlist(lapply(ga_markers$GeneActivity.feature,function(x) {
    gene_coord<-GRangesToString(LookupGeneCoords(obj,x))
    }))
  plt<-CoveragePlot(
    object = obj,
    region=ga_regions,
    #features = ga_markers$GeneActivity.feature,
    annotation = TRUE,
    extend.upstream = 1000, extend.downstream = 1000,
    peaks = TRUE,
    assay="peaks",
    #expression.assay="GeneActivity",
    group.by="cluster",nrow=1
  )
  ggsave(plt,file=paste0("./clonal_de_tfs/","clonal_",samp,"_GA_covplot.pdf"),height=length(unique(dat_epi$cluster)),width=length(ga_regions)*4,limitsize=F)
}

plot_de_clones_dotplot<-function(obj=dat_epi,samp=sample_in){
  de_markers<-Identify_Marker_TFs(obj,group_by.="cluster",assay.="SoupXRNA",prefix.=paste0("./clonal_de_tfs/","clonal_",samp))
  de_markers<-de_markers %>% group_by(SoupXRNA.group) %>% filter(SoupXRNA.logFC>0.2) %>% arrange(desc(SoupXRNA.auc)) %>% slice_head(n=10)
  plt<-DotPlot(
    object = obj,
    features = unique(de_markers$SoupXRNA.feature),
    group.by="cluster",
    assay="SoupXRNA", cluster.idents = TRUE) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_color_gradient2(low=plasma(3)[1],mid=plasma(3)[2],high=plasma(3)[3])
    
  ggsave(plt,file=paste0("./clonal_de_tfs/","clonal_",samp,"_DE_dotplot.pdf"),height=length(unique(dat_epi$cluster))/2+1,width=length(unique(de_markers$SoupXRNA.feature))/2,limitsize=F)
}

#make a directory for output
system("mkdir -p ./clonal_de_tfs")

#downloading PAM50 gene list from SCSubtype git repo
if(!file.exists("NatGen_Supplementary_table_S4.csv")){
  system("wget https://raw.githubusercontent.com/Swarbricklab-code/BrCa_cell_atlas/main/scSubtype/NatGen_Supplementary_table_S4.csv")
}
# read in scsubtype gene signatures
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv",col.names=c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

module_feats<-list()
module_feats[["Basal_SC"]]=as.vector(sigdat[,"Basal_SC"])
module_feats[["Her2E_SC"]]=as.vector(sigdat[,"Her2E_SC"])
module_feats[["LumA_SC"]]=as.vector(sigdat[,"LumA_SC"])
module_feats[["LumB_SC"]]=as.vector(sigdat[,"LumB_SC"])
module_feats<-lapply(module_feats,function(x) {x[x!=""]})


#run DE multimodal per sample clones
for (sample_in in unique(dat@meta.data[!is.na(dat@meta.data$cnv_defined_ploidy),]$sample)){
    print(paste("Running sample",sample_in))
    obj=subset(dat,sample==sample_in & cnv_defined_ploidy %in% c("diploid","aneuploid"))
    if(any(table(obj$cluster)<30)){
      print("Not enough cells per cluster. Skipping")
    } else {
        print((paste(sample_in,"has",as.character(ncol(obj)),"cells.")))
          print(paste("Calculating DE for sample",sample_in))
          run_top_TFs(obj,prefix=paste0("./clonal_de_tfs/","clonal_",sample_in),i="cluster",n_markers=10) #generate top TF markers
          if(obj$Diagnosis[1]=="IDC"){
          dat_epi<-subset(obj,cnv_defined_ploidy=="aneuploid" & HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")) #limit to IDC
          scsubtype<-run_scsubtype_perclone(dat_epi)
          #plot pam50 genes per cluster
          plt<-DotPlot(dat_epi,cluster.idents=T,group.by="cluster",dot.scale=5,scale.max=25,features=module_feats,assay="SCT") + 
                RotatedAxis()+ 
                scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-2,2)) + 
                theme(axis.text.x = element_text(size=6, angle=45)) + ggtitle(paste(sample_in,"RNA PAM50"))
          plt2<-DotPlot(dat_epi,cluster.idents=T,group.by="cluster",features=module_feats,assay="GeneActivity") + 
                RotatedAxis()+ 
                scale_color_gradient2(low="white",mid="lightblue",high="purple",limits=c(-1,1)) + 
                theme(axis.text.x = element_text(size=6, angle=45)) + ggtitle(paste(sample_in,"Gene Activity PAM50"))
          ggsave(plt/plt2,file=paste0("./clonal_de_tfs/","clonal_",sample_in,"_PAM50.pdf"),width=30,limitsize=F)
          plot_ga_clones_covplot(obj=dat_epi,samp=sample_in)
          plot_de_clones_dotplot(obj=dat_epi,samp=sample_in)
        }      
    }
}

#merge 
#gzip -f *motif.pdf
#convert `ls *pdf` cluster_summary.pdf
#slack -F cluster_summary.pdf ryan_todo

