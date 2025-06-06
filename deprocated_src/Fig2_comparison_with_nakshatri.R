```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(viridis)
library(optparse)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ComplexHeatmap)
library(dendextend)
library(ggdendro)
library(circlize)


#Identify top markers
Identify_Marker_TFs<-function(x,group_by="group",assay="shared_rna",pval_filt=1,assay_name="RNA"){
    x[[assay]]<-as(object = x[[assay]], Class = "Assay")
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by, 
      groups_use=unname(unlist(unique(x@meta.data[group_by]))),
      y=unname(unlist(unique(x@meta.data[group_by]))), 
      assay = 'data', seurat_assay = assay)
    markers<-markers[markers$padj<=pval_filt,]
    colnames(markers) <- paste(assay_name, colnames(markers),sep=".")
    if (assay == "chromvar") {
      motif.names <- markers[,paste0(assay_name,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="peaks")
    } else {
    markers$gene <- markers[,paste0(assay_name,".feature")]
    }
    return(markers) 
}

#Grab top overlapping TFs
topTFs <- function(markers_list,group_by, padj.cutoff = 1e-2,rna=NA,ga=NA,motifs=NA) {
  ctmarkers_rna <- dplyr::filter(rna, RNA.group == group_by) %>% 
    arrange(-RNA.auc)
    if(is.data.frame(motifs)) {
    ctmarkers_motif <- dplyr::filter(motifs, chromvar.group == group_by) %>% 
      arrange(-chromvar.auc)
    }
    if(is.data.frame(ga)) {
    ctmarkers_ga<- dplyr::filter(ga, GeneActivity.group == group_by) %>% 
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
  top_tfs$group<-group_by
  return(top_tfs)
}

#Average markers across groups
average_features<-function(x=out_subset,features=tf_$motif.feature,assay="chromvar",group_by){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    dat_motif<-x[[assay]]@data[features,]
    dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
    sum_motif<-split(dat_motif,x@meta.data[,group_by]) #group by rows to seurat clusters
    sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean,na.rm=T)) #take average across group
    sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame
    sum_motif<-t(scale(sum_motif))
    sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    sum_motif<-sum_motif[complete.cases(sum_motif),]
    return(sum_motif)
}

plot_top_tf_markers<-function(x=out_subset,group_by="group",prefix="lum_epi",n_markers=20,order_by_idents=TRUE){
    #define markers
    markers<-list(
        Identify_Marker_TFs(x=x,group_by=group_by,assay="shared_rna",assay_name="RNA"),
        Identify_Marker_TFs(x=x,group_by=group_by,assay="GeneActivity",assay_name="GeneActivity"),
        Identify_Marker_TFs(x=x,group_by=group_by,assay="chromvar",assay_name="chromvar"))
    names(markers)<-c("RNA","GeneActivity","chromvar")
    markers_out<-do.call("rbind",lapply(unique(x@meta.data[,group_by]),
        function(group) head(topTFs(markers_list=markers,group_by=group,
                        rna=markers$RNA,ga=markers$GeneActivity,motifs=markers$chromvar),
                        n=n_markers))) #grab top 5 TF markers per celltype
    markers_out<-markers_out[!duplicated(markers_out$gene),]
    dim(markers_out)
    #summarize markers
    tf_rna<-average_features(x=x,features=markers_out$gene,assay="RNA",group_by=group_by)
    tf_rna<-tf_rna[row.names(tf_rna) %in% tf_markers$gene,]
    tf_ga<-average_features(x=out_subset,features=markers_out$gene,assay="GeneActivity",group_by=group_by)
    tf_ga<-tf_ga[row.names(tf_ga) %in% markers_out$gene,]
    tf_motif<-average_features(x=out_subset,features=markers_out$chromvar.feature,assay="chromvar",group_by=group_by)
    tf_motif<-tf_motif[row.names(tf_motif) %in% markers_out$chromvar.feature,]
    row.names(tf_motif)<-markers_out[markers_out$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]
    average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering

    #set up heatmap seriation and order by RNA
    o_rows =dist(average_matrix) %>%
                          hclust() %>%
                          as.dendrogram()  #%>%
                          #ladderize()
    o_col =dist(t(average_matrix)) %>%
                      hclust() %>%
                      as.dendrogram()  %>%
                      ladderize()
    side_ha_rna<-data.frame(ga_motif=tf_markers[get_order(o_rows,1),]$RNA.auc)
    #colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.90,0.95)),plasma(3))
    colfun_rna=colorRamp2(c(0,1,2),plasma(3))

    side_ha_motif<-data.frame(chromvar_motif=tf_markers[get_order(o_rows,1),]$chromvar.auc)
    #colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.90,0.95)),cividis(3))
    colfun_motif=colorRamp2(c(0,1,2),cividis(3))

    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-tf_markers[tf_markers$gene %in% row.names(tf_motif),]$chromvar.feature
    plt<-MotifPlot(object = x,assay="peaks",motifs = motif_list[get_order(o_rows,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

    side_ha_ga<-data.frame(ga_auc=tf_markers[get_order(o_rows,1),]$GeneActivity.auc)
    #colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.90,0.95)),magma(3))
    colfun_ga=colorRamp2(c(0,1,2),magma(3))

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), labels =row.names(tf_rna),labels_gp=gpar(fontsize=6)))
    o_col<-if(order_by_idents){
        levels(Idents(out_subset))
    }else{
        as.factor(labels(o_col))}

    rna_auc<-Heatmap(side_ha_rna,
        cluster_rows = o_rows,
        col=side_ha_col,
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=7))

    rna_plot<-Heatmap(tf_rna,
        cluster_rows = o_rows,
        column_order=o_col,
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 8),
        show_row_names=FALSE,
        column_names_rot=90)

      ga_auc<-Heatmap(side_ha_ga,
          cluster_rows = o_rows,         
          col=side_ha_col,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      ga_plot<-Heatmap(tf_ga,
          cluster_rows = o_rows,                 
        column_order=o_col,
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90)

      motif_auc<-Heatmap(side_ha_motif,
          cluster_rows = o_rows,          
          col=side_ha_col,
          show_row_names=FALSE,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      motif_plot<-Heatmap(tf_motif,
          cluster_rows = o_rows,                 
        column_order=o_col,
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          #top_annotation=top_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
          right_annotation=gene_ha)
      
      #motif_image<-anno_image(paste0(prefix,".tf.heatmap.motif.svg"))

      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot,row_title=prefix)
    pdf(paste0(prefix,".tf.heatmap.pdf"),height=height)
    print(plt1)
    dev.off()
}


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#cistopic=readRDS(opt$cistopic)
#titan=readRDS(opt$titan)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.clone_annot.passqc.SeuratObject.rds"
dat=readRDS(opt$object_input)
ref<-readRDS(file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri/nakshatri_multiome.rds")
table(Idents(dat))
table(Idents(ref))
ref[["peaks"]]<-ref[["our_peaks"]] #rename for assay merging

#Filter RNA to common features between data sets
#Peaks, GeneActivity and Chromvar already only common features
length(row.names(dat[["GeneActivity"]])) == sum(row.names(dat[["GeneActivity"]]) %in% row.names(ref[["GeneActivity"]]))
length(row.names(dat[["chromvar"]])) == sum(row.names(dat[["chromvar"]]) %in% row.names(ref[["chromvar"]]))
length(row.names(dat[["peaks"]])) == sum(row.names(dat[["peaks"]]) %in% row.names(ref[["peaks"]]))


dat$data<-"primary"
ref$data<-"nakshatri"
shared_rna<-row.names(dat[["RNA"]])[row.names(dat[["RNA"]]) %in% row.names(ref[["RNA"]])]
shared_rna<-row.names(ref[["RNA"]])[row.names(ref[["RNA"]]) %in% shared_rna]

out<-merge(dat,ref)
out[["shared_rna"]]<-CreateAssayObject(counts=out[["RNA"]]@counts[shared_rna,])
DefaultAssay(out)<-"shared_rna"
out<-NormalizeData(out)
out<-FindVariableFeatures(out)
out<-ScaleData(out)

out$celltype<-out$assigned_celltype
out@meta.data[is.na(out@meta.data$celltype),]$celltype<-out@meta.data[is.na(out@meta.data$celltype),]$author_cell_type
out@meta.data[is.na(out@meta.data$Mol_Diagnosis),]$Mol_Diagnosis<-"nakshatri"
out$group<-paste(out$Diagnosis,out$Mol_Diagnosis)
Idents(out)<-factor(out$group,
levels=c("NA nakshatri","NAT nakshatri","DCIS DCIS",
"ILC ER+/PR+/HER2-","ILC ER+/PR-/HER2-",
"IDC ER+/PR-/HER2+","IDC ER+/PR+/HER2-","IDC ER+/PR-/HER2-"))

#saveRDS(out,file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri/nakshatri_ours_merged.SeuratObject.rds")


#lum epithelial (by molecular type)
out_subset<-subset(out,cells=row.names(out@meta.data[out@meta.data$celltype %in% c("cancer_luminal_epithelial","luminal_epithelial","LHS","LASP"),]))
plot_top_tf_markers(x=out_subset,group_by="group",prefix="lum_epi",n_markers=20)

#basal epithelial (by molecular type)
out_subset<-subset(out,cells=row.names(out@meta.data[out@meta.data$celltype %in% c("BM","basal_epithelial"),]))
plot_top_tf_markers(x=out_subset,group_by="group",prefix="basal_epi",n_markers=20)

#tcell (by molecular type)
out_subset<-subset(out,cells=row.names(out@meta.data[out@meta.data$celltype %in% c("tcell","T-cells"),]))
plot_top_tf_markers(x=out_subset,group_by="group",prefix="tcells",n_markers=20)

#fibroblast (by molecular type)
out_subset<-subset(out,cells=row.names(out@meta.data[out@meta.data$celltype %in% c("fibroblast","Fibroblast"),]))
plot_top_tf_markers(x=out_subset,group_by="group",prefix="fibroblast",n_markers=20)

#cancer lum epi (by clone)
out_subset<-subset(out,cells=row.names(out@meta.data[!is.na(out@meta.data$merge_cluster),]))
out_subset<-subset(out_subset,
    cells=row.names(out_subset@meta.data[out_subset@meta.data$merge_cluster %in% names(table(out_subset$merge_cluster))[table(out_subset$merge_cluster)>50],]))
Idents(out_subset)<-out_subset$merge_cluster
plot_top_tf_markers(x=out_subset,group_by="merge_cluster",prefix="clones",n_markers=20,order_by_idents=FALSE)

#our data cell type
out_subset<-subset(out,cells=row.names(out@meta.data[!is.na(out@meta.data$assigned_celltype),]))
Idents(out_subset)<-factor(out_subset$assigned_celltype,
    levels=c("cancer_luminal_epithelial","luminal_epithelial","basal_epithelial",
    "adipocyte","endothelial_lymphatic","endothelial_vascular",
    "pericyte","fibroblast","mast","myeloid","bcell","plasma","pDC","tcell"))
plot_top_tf_markers(x=out_subset,group_by="assigned_celltype",prefix="our_celltypes",n_markers=20,order_by_idents=TRUE)

#ref data cell type
out_subset<-subset(out,cells=row.names(out@meta.data[!is.na(out@meta.data$author_cell_type),]))
Idents(out_subset)<-factor(out_subset@meta.data$author_cell_type,
    levels=c("LASP","LHS","BM",
    "Adi-1","Adi-2","Endo-1","Endo-2","Fibroblasts","Macrophages","T-cells"))
plot_top_tf_markers(x=out_subset,group_by="author_cell_type",prefix="ref_celltypes",n_markers=20,order_by_idents=TRUE)
