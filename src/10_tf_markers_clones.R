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
library(presto,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(seriation)
library(org.Hs.eg.db)
library(dendextend)
library(msigdbr,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(fgsea,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(BSgenome.Hsapiens.UCSC.hg38)
library(dendextend)
library(ggdendro)
library(circlize)
library(ggtern)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);


####################################################
#           Fig 2 Molecular type Tornado Plot      #
###################################################
#similarity matrix for diagnoses on chromvar motifs
#define markers
#define markers
dat$cellstate<-paste(dat$assigned_celltype,dat$Diagnosis)
dat_sub<-subset(dat,assigned_celltype %in% c("cancer","luminal_hs")) #epi/cancer only
dat_sub<-subset(dat_sub, cells=row.names(dat_sub@meta.data[!(dat_sub@meta.data$merged_assay_clones %in% c("contamination")),]) )
Idents(dat_sub)<-dat_sub$Diag_MolDiag
dat_sub<-subset(dat_sub,downsample=1000)

table(Idents(dat_sub))

#use wilcox
da_peaks <- presto:::wilcoxauc.Seurat(X = dat_sub, group_by = "Diag_MolDiag", 
  y=Idents(dat_sub), 
  assay = 'data', seurat_assay = "ATAC")

#use LR for peaks with counts as latent vars
#da_peaks<-FindAllMarkers(dat_sub,assay="ATAC",test.use="LR",latent.vars="nCount_ATAC",min.pct=0.1)

#add chromvar motif scan of motifs per DA peaks
da_peaks<-merge(da_peaks,as.data.frame(dat@assays$ATAC@motifs@data),by.x="feature",by.y= 'row.names')
write.table(da_peaks,file="epithelial_diagnosis_da_peaks.csv",col.names=T,row.names=T,sep=",")

clin_col=c(
"DCIS DCIS"="#cccccb", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"IDC ER+/PR-/HER2+"="#f37872", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"NAT NA"="#c2d9ea")

#setup outnames to be more nice
outnames=gsub(x=names(clin_col),pattern="[+]",replacement="plus")
outnames=gsub(x=outnames,pattern="[-]",replacement="neg")
outnames=gsub(x=outnames,pattern="[/]",replacement="_")
outnames=gsub(x=outnames,pattern="[ ]",replacement="_")
outnames=setNames(nm=names(clin_col),outnames)


#all DA peaks
tornado_plot<-function(obj=dat_sub,da_peak_set=da_peaks,i,peak_count=100){
    print(i)
    top_peaks <- da_peak_set %>% filter(group==i) %>%  filter(logFC>0) %>% arrange(padj) %>% slice_head(n=peak_count)
    print(head(top_peaks))

    obj<-RegionMatrix(obj,key="DA_mat",
    regions=StringToGRanges(top_peaks$feature),
    upstream=5000,downstream=5000,
    assay="ATAC")

    plt<-RegionHeatmap(obj,key="DA_mat",
    upstream=5000,downstream=5000,
    order=TRUE, window=(10000)/100,normalize=FALSE,
    assay="ATAC", idents=levels(Idents(obj)),
    nrow=length(unique(da_peak_set$group)))+ 
    #scale_fill_gradient2(low="black",mid=unname(clin_col[i]),midpoint=0.04,high="white",na.value="black",breaks=seq(0,0.05,0.005))+
    ggtitle(i)
    #cols=setNames(nm=c("ATAC"),c(unname(clin_col[i]))),

    print("Returning plot...")
    return(plt)
}

#all peaks
plt_list<-lapply(names(clin_col),function(x) {
    tornado_plot(obj=dat_sub,da_peak_set=da_peaks,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1,guides='collect')
ggsave(plt,file="tornado.all_da_peaks.epi_diag.pdf",width=20,height=20)

#da peaks that overlap with ESR1 motif only
motif_name="ESR1"
motif<-names(dat@assays$ATAC@motifs@motif.names[which(dat@assays$ATAC@motifs@motif.names==motif_name)])
da_peaks_motif_filt<-da_peaks[da_peaks$gene %in% names(which(dat@assays$ATAC@motifs@data[,motif])),]

plt_list<-lapply(unique(da_peaks_motif_filt$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks_motif_filt,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1)
ggsave(plt,file="tornado.ESR1_overlap_da_peaks.diag.pdf",width=20,height=20)



####################################################
#           Fig 2  TF motifs IDC v ILC           #
###################################################
#define markers
dat$cellstate<-paste(dat$assigned_celltype,dat$Diagnosis)
dat_sub<-subset(dat,assigned_celltype %in% c("cancer","luminal_hs")) #epi/cancer only
dat_sub<-subset(dat_sub, cells=row.names(dat_sub@meta.data[!(dat_sub@meta.data$merged_assay_clones %in% c("contamination")),]) )
dat_sub<-subset(dat_sub,cellstate %in% c("cancer IDC","cancer ILC") )
dat_sub<-subset(dat_sub,Diagnosis %in% c("ILC","IDC"))
Idents(dat_sub)<-dat_sub$sample

dat_sub<-subset(dat_sub,downsample=100)
Idents(dat_sub)<-dat_sub$Diagnosis
dat_sub<-subset(dat_sub,downsample=310) #min ILC

#following https://stuartlab.org/signac/articles/motif_vignette
markers <- FindMarkers(
  object =dat_sub,
  assay="chromvar",
  ident.1 = "IDC",
  ident.2 = 'ILC',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  logfc.threshold=0,
  fc.name = "avg_diff"
)
motif.names <- row.names(markers)
markers$gene <- ConvertMotifID(dat_sub, id = motif.names,assay="ATAC") #or ATAC as assay

#filter motifs by those differentially expressed as well
rna_markers <- presto:::wilcoxauc.Seurat(X = dat_sub, group_by = "Diagnosis", 
  y=c("IDC","ILC"), 
  assay = 'data', seurat_assay = "SCT")
rna_markers_idc <- rna_markers %>% filter(padj<0.05) %>% filter(group=="IDC") %>% filter(logFC>0)
rna_markers_ilc <- rna_markers %>% filter(padj<0.05) %>% filter(group=="ILC") %>% filter(logFC>0)

chromvar_markers_idc <- markers %>% filter(avg_diff > 0)  %>% filter(gene %in% rna_markers_idc$feature)
chromvar_markers_ilc <- markers %>% filter(avg_diff < 0)  %>% filter(gene %in% rna_markers_ilc$feature)

markers <- rbind(chromvar_markers_idc,chromvar_markers_ilc)
markers$fill_col<-"gray"
markers[markers$p_val_adj<0.05 & markers$avg_diff>0.25,]$fill_col<-"idc"
markers[markers$p_val_adj<0.05 & markers$avg_diff<as.numeric(-0.25),]$fill_col<-"ilc"

markers$label_name<-NA

min_label<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_diff<0)%>% slice_min(avg_diff,n=20)
max_label<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_diff>0) %>% slice_max(avg_diff,n=20)

markers[row.names(min_label),]$label_name<-min_label$gene
markers[row.names(max_label),]$label_name<-max_label$gene
plt<-ggplot(markers,aes(x=avg_diff,y=-log10(p_val_adj),color=fill_col,label=label_name))+geom_point()+geom_text_repel(max.overlaps=Inf)+theme_minimal()
ggsave(plt,file="ILC_IDC.chromvar_DE.pdf")

markers<-Identify_Marker_TFs(x=dat_sub,group_by="cellstate",assay="chromvar",assay_name="chromvar")
markers_out <- markers %>% group_by(chromvar.group) %>% filter(abs(chromvar.logFC)>2)

length(unique(markers_out$chromvar.feature))
markers_out<-markers_out[!duplicated(markers_out$chromvar.feature),]
tf_motif<-average_features(x=dat_sub,features=markers_out$chromvar.feature,assay="chromvar",group_by="cellstate")
row.names(tf_motif)<-markers_out$gene
colfun_motif=colorRamp2(c(0,0.5,1),cividis(3))
cell_counts=as.data.frame(table(dat_sub$cellstate))
row.names(cell_counts)<-cell_counts$Var1
row_ha = rowAnnotation(cell_count = anno_barplot(log10(cell_counts[colnames(tf_motif),]$Freq)))

motif_plot<-Heatmap(cor(tf_motif,method="spearman"),
    name="TF Motif",
    column_title="TF Motif",
    col=colfun_motif,
    column_names_gp = gpar(fontsize = 5),
    row_names_gp = gpar(fontsize = 5),
    show_row_names=TRUE,
    column_names_rot=90, right_annotation = row_ha)
pdf("test_motif.distance.pdf")
print(motif_plot)
dev.off()
#69 motifs

####################################################
#           Fig 3 Heatmap By Clones                #
###################################################
#same functions from bc_multiome_nf_analysis/src/8_tf_markers_coverage_plots.R

#Identify top markers
Identify_Marker_TFs<-function(x,group_by,assay,pval_filt=1,assay_name){
      if (assay != "chromvar") {
        x[[assay]]<-as(object = x[[assay]], Class = "Assay")
        }
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by, 
      groups_use=unname(unlist(unique(x@meta.data[group_by]))),
      y=unname(unlist(unique(x@meta.data[group_by]))), 
      assay = 'data', seurat_assay = assay)
    markers<-markers[markers$padj<=pval_filt,]
    colnames(markers) <- paste(assay_name, colnames(markers),sep=".")
    if (assay == "chromvar") {
      motif.names <- markers[,paste0(assay_name,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="ATAC") #or ATAC as assay
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
average_features<-function(x=out_subset,features=tf_$motif.feature,assay,group_by){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    x[[assay]]<-as(object = x[[assay]], Class = "Assay")
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

plot_top_tf_markers<-function(x=out_subset,group_by,prefix,n_markers=20,order_by_idents=TRUE){
    #define markers
    markers<-list(
        Identify_Marker_TFs(x=x,group_by=group_by,assay="RNA",assay_name="RNA"),
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
    tf_rna<-tf_rna[row.names(tf_rna) %in% markers_out$gene,]
    tf_ga<-average_features(x=x,features=markers_out$gene,assay="GeneActivity",group_by=group_by)
    tf_ga<-tf_ga[row.names(tf_ga) %in% markers_out$gene,]
    tf_motif<-average_features(x=x,features=markers_out$chromvar.feature,assay="chromvar",group_by=group_by)
    tf_motif<-tf_motif[row.names(tf_motif) %in% markers_out$chromvar.feature,]
    row.names(tf_motif)<-markers_out[markers_out$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]
    average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering
    #average_matrix=tf_ga #just cluster by GENE activity expression
    #set up heatmap seriation and order by GA
    o_rows =dist(average_matrix) %>%
                          hclust() %>%
                          as.dendrogram()  #%>%
                          #ladderize()
    o_col =dist(t(average_matrix),method="maximum") %>%
                      hclust() %>%
                      as.dendrogram()  %>%
                      ladderize()
    side_ha_rna<-data.frame(ga_motif=markers_out[get_order(o_rows,1),]$RNA.auc)
    #colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.90,0.95)),plasma(3))
    colfun_rna=colorRamp2(c(0,1,2),plasma(3))

    side_ha_motif<-data.frame(chromvar_motif=markers_out[get_order(o_rows,1),]$chromvar.auc)
    #colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.90,0.95)),cividis(3))
    colfun_motif=colorRamp2(c(0,1,2),cividis(3))

    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-markers_out[markers_out$gene %in% markers_list,]$chromvar.feature
    plt<-MotifPlot(object = x,assay="ATAC",motifs = motif_list[get_order(o_rows,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

    side_ha_ga<-data.frame(ga_auc=markers_out[get_order(o_rows,1),]$GeneActivity.auc)
    #colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.90,0.95)),magma(3))
    colfun_ga=colorRamp2(c(0,1,2),magma(3))

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), labels =row.names(tf_rna),labels_gp=gpar(fontsize=6)))
    o_col<-if(order_by_idents){
        levels(Idents(x))
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
    
    pdf(paste0(prefix,".tf.heatmap.pdf"))
    print(draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot,row_title=prefix))
    dev.off()
}

#merged clone calls
underpower_clones<-names(which(table(dat$merged_assay_clones)<30))
dat@meta.data[dat@meta.data$merged_assay_clones %in% underpower_clones,]$merged_assay_clones<-NA
table(dat$merged_assay_clones)

dat_cnv<-subset(dat,cells=names(dat$merged_assay_clones[!is.na(dat$merged_assay_clones)]))
plot_top_tf_markers(x=dat_cnv,group_by="merged_assay_clones",prefix="clones",n_markers=3,order_by_idents=FALSE)


####################################################
#           Fig 2 DA TF motifs IDC v ILC           #
###################################################

#similarity matrix for diagnoses on chromvar motifs
#define markers
dat_sub<-subset(dat_epi,cellstate %in% c("cancer IDC","cancer ILC") )
dat_sub<-subset(dat_sub,Diagnosis %in% c("ILC","IDC"))
Idents(dat_sub)<-dat_sub$sample
dat_sub<-subset(dat_sub,downsample=100)
Idents(dat_sub)<-dat_sub$cellstate

#following https://stuartlab.org/signac/articles/motif_vignette
markers <- FindMarkers(
  object =dat_sub,
  assay="chromvar",
  ident.1 = "cancer IDC",
  ident.2 = 'cancer ILC',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  logfc.threshold=0,
  fc.name = "avg_diff"
)
motif.names <- row.names(markers)
markers$gene <- ConvertMotifID(dat_epi, id = motif.names,assay="ATAC") #or ATAC as assay

#filter motifs by those differentially expressed as well
rna_markers <- presto:::wilcoxauc.Seurat(X = dat_sub, group_by = "cellstate", 
  groups_use=c("cancer IDC","cancer ILC"),
  y=c("cancer IDC","cancer ILC"), 
  assay = 'data', seurat_assay = "SCT")
rna_markers_idc <- rna_markers %>% filter(padj<0.05) %>% filter(group=="cancer IDC") %>% filter(logFC>0)
rna_markers_ilc <- rna_markers %>% filter(padj<0.05) %>% filter(group=="cancer ILC") %>% filter(logFC>0)

chromvar_markers_idc <- markers %>% filter(avg_diff > 0)  %>% filter(gene %in% rna_markers_idc$feature)
chromvar_markers_ilc <- markers %>% filter(avg_diff < 0)  %>% filter(gene %in% rna_markers_ilc$feature)

markers <- rbind(chromvar_markers_idc,chromvar_markers_ilc)
markers$fill_col<-"gray"
markers[markers$p_val_adj<0.05 & markers$avg_diff>0.25,]$fill_col<-"idc"
markers[markers$p_val_adj<0.05 & markers$avg_diff<as.numeric(-0.25),]$fill_col<-"ilc"

markers$label_name<-NA

min_label<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_diff<0)%>% slice_min(avg_diff,n=10)
max_label<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_diff>0) %>% slice_max(avg_diff,n=10)

markers[row.names(min_label),]$label_name<-min_label$gene
markers[row.names(max_label),]$label_name<-max_label$gene
plt<-ggplot(markers,aes(x=avg_diff,y=-log10(p_val_adj),color=fill_col,label=label_name))+geom_point()+geom_text_repel(max.overlaps=Inf)+theme_minimal()
ggsave(plt,file="ILC_IDC.chromvar_DE.pdf")

markers<-Identify_Marker_TFs(x=dat_sub,group_by="cellstate",assay="chromvar",assay_name="chromvar")
markers_out <- markers %>% group_by(chromvar.group) %>% filter(abs(chromvar.logFC)>2)


length(unique(markers_out$chromvar.feature))
markers_out<-markers_out[!duplicated(markers_out$chromvar.feature),]
tf_motif<-average_features(x=dat_sub,features=markers_out$chromvar.feature,assay="chromvar",group_by="cellstate")
row.names(tf_motif)<-markers_out$gene
colfun_motif=colorRamp2(c(0,0.5,1),cividis(3))
cell_counts=as.data.frame(table(dat_sub$cellstate))
row.names(cell_counts)<-cell_counts$Var1
row_ha = rowAnnotation(cell_count = anno_barplot(log10(cell_counts[colnames(tf_motif),]$Freq)))

motif_plot<-Heatmap(cor(tf_motif,method="spearman"),
    name="TF Motif",
    column_title="TF Motif",
    col=colfun_motif,
    column_names_gp = gpar(fontsize = 5),
    row_names_gp = gpar(fontsize = 5),
    show_row_names=TRUE,
    column_names_rot=90, right_annotation = row_ha)
pdf("test_motif.distance.pdf")
print(motif_plot)
dev.off()
#69 motifs


####################################################
#           Extra Plots                         #
###################################################

#epi only
dat_cancer<-subset(dat,assigned_celltype %in% c("cancer","luminal_hs","luminal_asp","basal_myoepithelial"))
dat_cancer$cellstate<-paste(dat_cancer$assigned_celltype,dat_cancer$Diagnosis,dat_cancer$Mol_Diagnosis)
plot_top_tf_markers(x=dat_cancer,group_by="cellstate",prefix="epi_celltypes",n_markers=20,order_by_idents=FALSE)

#cancer only by diag
dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer$cellstate<-paste(dat_cancer$assigned_celltype,dat_cancer$Diagnosis,dat_cancer$Mol_Diagnosis)
plot_top_tf_markers(x=dat_cancer,group_by="cellstate",prefix="cancer_diag",n_markers=20,order_by_idents=FALSE)

#epi only with ref
bhat_nakshatri<-readRDS(opt$ref_object)
table(bhat_nakshatri$cell_type)
bhat_nakshatri$assigned_celltype<-bhat_nakshatri$cell_type
bhat_nakshatri[["ATAC"]]<-bhat_nakshatri[["our_peaks"]]
DefaultAssay(bhat_nakshatri)<-"ATAC"
bhat_nakshatri<-subset(bhat_nakshatri,assigned_celltype %in% c("luminal adaptive secretory precursor cell of mammary gland",
"luminal hormone-sensing cell of mammary gland","basal-myoepithelial cell of mammary gland"))
bhat_nakshatri$Diagnosis<-"normal"
bhat_nakshatri$Mol_Diagnosis<-"normal"

dat_epi<-merge(
  subset(dat,assigned_celltype %in% c("cancer","luminal_hs","luminal_asp","basal_myoepithelial")),
  bhat_nakshatri)

dat_epi$cellstate<-paste(dat_epi$assigned_celltype,dat_epi$Diagnosis)

dat_epi<-subset(dat_epi,cellstate %in% names(table(dat_epi$cellstate))[which(table(dat_epi$cellstate)>200)])
dat_epi[["RNA"]]<-JoinLayers(dat_epi[["RNA"]])
dat_epi<-ScaleData(dat_epi,assay="RNA")
dat_epi[["ATAC"]]@motifs<-dat_epi@assays$our_peaks@motifs
plot_top_tf_markers(x=dat_epi,group_by="cellstate",prefix="nakshatri_epi_diag",n_markers=5,order_by_idents=FALSE)

Idents(dat_epi)<-dat_epi$cellstate
markers <- presto:::wilcoxauc.Seurat(X = dat_epi, group_by = "cellstate", 
  groups_use=unname(unlist(unique(dat_epi$cellstate))),
  y=unname(unlist(unique(dat_epi$cellstate))), 
  assay = 'data', seurat_assay = "ATAC")

da_peaks <- markers %>% filter(padj<0.01) %>% filter(logFC>0.25)
#[1] 39525    10

#add chromvar motif scan of motifs per DA peaks
da_peaks<-merge(da_peaks,as.data.frame(dat@assays$ATAC@motifs@data),by.x="feature",by.y= 'row.names')
#add closest gene
annot<-ClosestFeature(
  object = dat_epi,
  regions = da_peaks$feature
)
da_peaks<-merge(da_peaks,annot,by.x="feature",by.y= 'query_region')
da_peaks <- da_peaks %>% filter(distance<2000)

write.table(da_peaks,file="epithelial_diagnosis_da_peaks.csv",col.names=T,row.names=T,sep=",")

library(GeneNMF)
library(rGREAT)

# #GSEA of clone DMRs
gsea_enrichment<-function(prefix,dmrs,
gene_universe,
category="C3",subcategory="TFT:GTRD",
out_setname="TFT"){
  top_p_gsea <- do.call("rbind",
    lapply(unique(dmrs$group), 
    function(i) {
    #gene set
    gene_list<-dmrs |> dplyr::filter(group==i) |> pull(gene_name)
    out<-runGSEA(gene_list, universe=gene_universe, category = category,subcategory=subcategory)
    out$celltype<-i
    return(out)
    }
    ))
 pltdat<-top_p_gsea %>% group_by(celltype) %>% slice_max(order_by = -padj, n = 5)
 plt<-ggplot(pltdat,aes(x=celltype,y=pathway))+
 geom_point(aes(size = -log10(padj), fill = overlap/size), shape=21)+
 theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 ggsave(plt,file=paste0(prefix,"_GSEA","_",out_setname,".pdf"),width=20)
 saveRDS(top_p_gsea,file=paste0(prefix,"_GSEA","_",out_setname,".rds"))
 }


  #set up gene universe on atac peaks
  prefix="epi_diag"
  dmrs=da_peaks
  gene_universe<-dmrs %>% pull(gene_name)
  gene_universe<-gene_universe[!duplicated(gene_universe)]

  #run gsea enrichment at gene level on different sets
   gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
   category="C3",subcategory="TFT:GTRD",out_setname="TFT") #find enrichment in tft (transcription factor targets)

   gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
   category="C1",subcategory=NULL,out_setname="position") #find enrichment in c1 signatures (positional)

   gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
   category="H",subcategory=NULL,out_setname="hallmark") #find cancer hallmark signatures

#run on RNA
rna<-FindAllMarkers(dat_epi,assay="RNA")
rna$group<-rna$cluster
rna$gene_name<-rna$gene

#set up gene universe on RNA
prefix="epi_diag_rna"
dmrs=rna
gene_universe<-dmrs %>% pull(gene_name)
gene_universe<-gene_universe[!duplicated(gene_universe)]

#run gsea enrichment at gene level on different sets
gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
category="C3",subcategory="TFT:GTRD",out_setname="TFT") #find enrichment in tft (transcription factor targets)

gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
category="C1",subcategory=NULL,out_setname="position") #find enrichment in c1 signatures (positional)

gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
category="H",subcategory=NULL,out_setname="hallmark") #find cancer hallmark signatures





#EPI
Idents(dat_epi)<-dat_epi$cellstate
clin_col=c(
"cancer ILC ER+/PR+/HER2-"="#f6bea1", 
"cancer ILC ER+/PR-/HER2-"="#b9db98", 
"cancer IDC ER+/PR-/HER2+"="#f37872", 
"cancer IDC ER+/PR+/HER2-"="#8d86c0", 
"cancer IDC ER+/PR-/HER2-"="#7fd0df", 
"luminal_hs"="#4c3c97",
"luminal_asp"="#7161ab",
"basal_myoepithelial"="#ee6fa0",
"basal-myoepithelial cell of mammary gland"="#c497a9",
"luminal adaptive secretory precursor cell of mammary gland"="#6c6585",
"luminal hormone-sensing cell of mammary gland"="#6f6799")


marker_list<-FindAllMarkers(dat_epi,assay="GeneActivity")
#dat_epi<-PrepSCTFindMarkers(dat_epi)
#rna_marker_list<-FindAllMarkers(dat_epi,assay="SCT")

marker_subset<-marker_list %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1) %>% 
group_by(cluster) %>% 
slice_max(n=20,order_by=avg_log2FC,with_ties=FALSE)
#filter by RNA sig
#marker_subset[marker_subset$gene %in% rna_marker_list$gene,]

marker<-list()
#iglesia et al
marker[["her2_tumors"]]<-c("ERBB2","GRB7")
marker[["lum_tumor"]]<-c("ESR1","GRHL2")
marker[["basal-like"]]<-c("SFRP1","KRT17")
#nakshatri et al
marker[["tcells"]]<-c("IL7R","IFNG","GZMK","FCGR3A")
marker[["macrophage"]]<-c("LYVE1","ACKR1")
marker[["fibroblast"]]<-c("CXCL12","GLI2")
marker[["lum_hs"]]<-c("ESR1","FOXA1","GATA3")
marker[["lum_hs_asp"]]<-c("EHF","ELF5","KIT")
marker[["basal"]]<-c("TP63","KRT14")


# first compute the GC content for each peak
DefaultAssay(dat_epi)<-"ATAC"
dat_epi <- RegionStats(dat_epi, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
dat_epi <- LinkPeaks(
  object = dat_epi,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = unlist(marker_subset))

for(i in 1:nrow(marker_subset)){
    x<-marker_subset$gene[i]
    if (x in row.names(dat_epi@assays$SCT@data)){
      celltype<-marker_subset$cluster[i]
      col<-clin_col[celltype]
      celltype<-gsub(celltype,pattern="[+]",replacement="plus")
      celltype<-gsub(celltype,pattern="[-]",replacement="minus")
      celltype<-gsub(celltype,pattern="[ ]",replacement="_")
      celltype<-gsub(celltype,pattern="[/]",replacement="_")
      print(paste("Plotting...",x))
      print(celltype)

      cov_plot <- CoveragePlot(
        object = dat_epi,
        region = x,
        annotation = TRUE,
        peaks = TRUE,links=FALSE)+ scale_color_manual(values=col)

      link_plot <- LinkPlot(
        object = dat_epi,
        region = x)+scale_color_gradient2(limits=c(0,0.3),low="white",high="magenta")

      expr_plot <- ExpressionPlot(
        object = dat_epi,
        features = x,
        assay = "SCT") + scale_color_manual(values=col)

      plt<-CombineTracks(
        plotlist = list(cov_plot, link_plot),
        expression.plot = expr_plot,
        heights = c( 10, 3),
        widths = c(10, 1))
      ggsave(plt,file=paste0("coverage.",celltype,".",x,".pdf"))
    }
  }

da_peaks<-FindAllMarkers(dat_epi,assay="ATAC")


####################################################
#           Fig 2 Tornado Plots                  #
###################################################
#write out DA sites for tornado plots
Idents(dat)<-dat$Diag_MolDiag
dat_cancer<-subset(dat,assigned_celltype %in% c("cancer","luminal_hs","luminal_asp","basal_myoepithelial"))

#split out cancer cells by diagnosis, keep all others
dat_cancer$cellstate<-unlist(lapply(1:nrow(dat_cancer@meta.data),function(x){
  cell_tmp=dat_cancer@meta.data[x,]
  if(cell_tmp$assigned_celltype!="cancer"){
    return(cell_tmp$assigned_celltype)
  }else{
    return(paste(cell_tmp$Diagnosis,cell_tmp$assigned_celltype))
  }
}))

Idents(dat_cancer)<-factor(dat_cancer$cellstate,
    levels=c("luminal_hs","luminal_asp","basal_myoepithelial","ILC cancer","IDC cancer"))

dat_cancer<-subset(x = dat_cancer, downsample = 1500) #downsample to ~ lowest cell count
da_peaks<-FindAllMarkers(dat_cancer,assay="ATAC")

#add chromvar motif scan of motifs per DA peaks
da_peaks<-merge(da_peaks,as.data.frame(dat@assays$ATAC@motifs@data),by.x="gene",by.y= 'row.names')
write.table(da_peaks,file="epithelial_diagnosis_da_peaks.csv",col.names=T,row.names=T,sep=",")

#setup outnames to be more nice
outnames=gsub(x=names(clin_col),pattern="[+]",replacement="plus")
outnames=gsub(x=outnames,pattern="[-]",replacement="neg")
outnames=gsub(x=outnames,pattern="[/]",replacement="_")
outnames=gsub(x=outnames,pattern="[ ]",replacement="_")
outnames=setNames(nm=names(clin_col),outnames)

#making it smaller by removing other assays
DefaultAssay(dat_cancer)<-"ATAC"

#subsetting object 

table(Idents(dat_cancer))

volcano_col=c(
"luminal_hs"="#4c3c97",
"luminal_asp"="#7161ab",
"basal_myoepithelial"="#ee6fa0",
"NAT cancer"="#99CCFF",
"DCIS cancer"="#CCCCCC",
"IDC cancer"="#FF9966",
"ILC cancer"="#006633")

clin_col=c(
"DCIS DCIS"="#cccccb", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"IDC ER+/PR-/HER2+"="#f37872", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"NAT NA"="#c2d9ea")

#all DA peaks
tornado_plot<-function(obj=dat_cancer,da_peak_set=da_peaks,i=da_peaks$cluster,peak_count=100){
    print(i)
    top_peaks <- da_peak_set %>% filter(cluster==i) %>%  filter(avg_log2FC>0) %>% arrange(p_val_adj) %>% slice_head(n=peak_count)
    print(head(top_peaks))

    obj<-RegionMatrix(obj,key="DA_mat",
    regions=StringToGRanges(top_peaks$gene),
    upstream=5000,downstream=5000,
    assay="ATAC")

    plt<-RegionHeatmap(obj,key="DA_mat",
    upstream=5000,downstream=5000,
    order=TRUE, window=(10000)/100,normalize=TRUE,
    assay="ATAC", idents=levels(Idents(obj)),
    nrow=length(unique(da_peak_set$cluster)))+ 
    scale_fill_gradient2(low="black",mid=unname(volcano_col[i]),midpoint=0.03,high="white",na.value="black",breaks=seq(0,0.05,0.01))+
    ggtitle(i)
    #cols=setNames(nm=c("ATAC"),c(unname(clin_col[i]))),

    print("Returning plot...")
    return(plt)
}

#all peaks
plt_list<-lapply(unique(da_peaks$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1,guides='collect')
ggsave(plt,file="tornado.all_da_peaks.diag.pdf",width=20,height=20)

#da peaks that overlap with ESR1 motif only
motif_name="ESR1"
motif<-names(dat@assays$ATAC@motifs@motif.names[which(dat@assays$ATAC@motifs@motif.names==motif_name)])
da_peaks_motif_filt<-da_peaks[da_peaks$gene %in% names(which(dat@assays$ATAC@motifs@data[,motif])),]

plt_list<-lapply(unique(da_peaks_motif_filt$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks_motif_filt,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1)
ggsave(plt,file="tornado.ESR1_overlap_da_peaks.diag.pdf",width=20,height=20)


#da peaks that overlap with ESR1 motif only
motif_name="GRHL1"
motif<-names(dat@assays$ATAC@motifs@motif.names[which(dat@assays$ATAC@motifs@motif.names==motif_name)])
da_peaks_motif_filt<-da_peaks[da_peaks$gene %in% names(which(dat@assays$ATAC@motifs@data[,motif])),]

plt_list<-lapply(unique(da_peaks_motif_filt$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks_motif_filt,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1)
ggsave(plt,file="tornado.GRHL1_overlap_da_peaks.diag.pdf",width=20,height=20)


#da peaks that overlap with ESR1 motif only
motif_name="FOXA1"
motif<-names(dat@assays$ATAC@motifs@motif.names[which(dat@assays$ATAC@motifs@motif.names==motif_name)])
da_peaks_motif_filt<-da_peaks[da_peaks$gene %in% names(which(dat@assays$ATAC@motifs@data[,motif])),]

plt_list<-lapply(unique(da_peaks_motif_filt$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks_motif_filt,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1)
ggsave(plt,file="tornado.FOXA1_overlap_da_peaks.diag.pdf",width=20,height=20)


