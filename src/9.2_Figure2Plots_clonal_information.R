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
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character"),
  make_option(c("-r", "--ref_object"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri/nakshatri_multiome.geneactivity.rds", 
              help="Nakshatri reference object for epithelial comparisons", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
dat<-subset(dat,merged_assay_clones != "contamination")

clin_col=c("IDC ER+/PR-/HER2+"="#f37872", 
"DCIS"="#cccccb", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"NAT NA"="#c2d9ea")

write.table(dat@meta.data,row.names=T,col.names=T,sep="\t",file="cell_metadata.tsv")

#####Plot of heatmap of all clones ####

cnv_col<-c("0"="#002C3E", "0.5"="#78BCC4", "1"="#F7F8F3", "1.5"="#F7444E", "2"="#aa1407", "3"="#440803")
#from Curtis et al.

dat_sub<-subset(dat,cells=names(which(!is.na(dat$merged_assay_clones))))
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


####################################################
#           Fig Clone Heatmap                      #
###################################################


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

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Average markers across groups
average_features<-function(x=out_subset,features=tf_$motif.feature,assay,group_by,slot_name="data"){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    x[[assay]]<-as(object = x[[assay]], Class = "Assay")
    if(slot_name=="data"){
      dat_motif<-x[[assay]]@data[features,]
      dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
      sum_motif<-split(dat_motif,x@meta.data[,group_by]) #group by rows to seurat clusters
      sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean,na.rm=T)) #take mean across group
      sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame
      sum_motif<-t(scale(sum_motif))
      sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    } else {
      dat_motif<-x[[assay]]@counts[features,]
      dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
      sum_motif<-split(dat_motif,x@meta.data[,group_by]) #group by rows to seurat clusters
      sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean)) #take mean across group
      sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame
      sum_motif<-t(sum_motif)
      sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    }

    #sum_motif<-sum_motif[complete.cases(sum_motif),]
    return(sum_motif)
}

#modified to plot by sample and group by pairwise (also added a column annotation)
plot_top_tf_markers<-function(x=out_subset,group_by,prefix,n_markers=20,order_by_idents=FALSE,plot_by,outdir){
    gtf<-rtracklayer::readGFF(file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
    gtf<-makeGRangesFromDataFrame(gtf,keep.extra.columns=TRUE)
    gtf<-gtf[gtf$type=="gene",]
    gtf<-gtf[gtf@seqnames %in% paste0("chr",1:22),]

    #define markers
    markers<-list(
        Identify_Marker_TFs(x=x,group_by=group_by,assay="RNA",assay_name="RNA"),
        Identify_Marker_TFs(x=x,group_by=group_by,assay="GeneActivity",assay_name="GeneActivity"),
        Identify_Marker_TFs(x=x,group_by=group_by,assay="chromvar",assay_name="chromvar"))
    names(markers)<-c("RNA","GeneActivity","chromvar")
    markers_out<-do.call("rbind",lapply(unique(x@meta.data[,group_by]),
        function(group) head(topTFs(markers_list=markers,group_by=group,
                        rna=markers$RNA,ga=markers$GeneActivity,motifs=markers$chromvar),
                        n=3))) #grab top N TF markers per celltype
    markers_out<-markers_out[!duplicated(markers_out$gene),]
    dim(markers_out)

    #summarize markers over groups
    tf_rna<-average_features(x=x,features=markers_out$gene,assay="RNA",group_by=plot_by)
    tf_rna<-tf_rna[row.names(tf_rna) %in% markers_out$gene,]
    tf_ga<-average_features(x=x,features=markers_out$gene,assay="GeneActivity",group_by=plot_by)
    tf_ga<-tf_ga[row.names(tf_ga) %in% markers_out$gene,]
    tf_motif<-average_features(x=x,features=markers_out$chromvar.feature,assay="chromvar",group_by=plot_by)
    tf_motif<-tf_motif[row.names(tf_motif) %in% markers_out$chromvar.feature,]
    row.names(tf_motif)<-markers_out[markers_out$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_rna),row.names(tf_ga),gtf$gene_name))

    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]

    #assign cnv windows to genes
    cnv_granges<-makeGRangesFromDataFrame(data.frame(
        seqnames=unlist(lapply(strsplit(row.names(x@assays$cnv$data),"-"),"[",1)),
        start=unlist(lapply(strsplit(row.names(x@assays$cnv$data),"-"),"[",2)),
        end=unlist(lapply(strsplit(row.names(x@assays$cnv$data),"-"),"[",3))))
    gtf<-gtf[gtf$gene_name %in% markers_list,]
    overlaps<-findOverlaps(gtf,cnv_granges)
    overlaps<-overlaps[!duplicated(overlaps@to),]
    gtf$cnv_windows<-NA
    names(gtf)<-1:length(gtf)
    gtf[overlaps@from,]$cnv_windows<-overlaps@to
    gtf<-gtf[!is.na(gtf$cnv_windows),]
    gtf<-gtf[!duplicated(gtf$cnv_windows),]

    tf_rna<-tf_rna[gtf$gene_name,]
    tf_motif<-tf_motif[gtf$gene_name,]
    tf_ga<-tf_ga[gtf$gene_name,]
    tf_cnv<-average_features(x=x,features=row.names(dat@assays$cnv@counts)[gtf$cnv_windows],assay="cnv",slot_name="counts",group_by=plot_by)

    #average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering
    average_matrix=tf_motif #just cluster only on TF motifs

    #set up heatmap seriation and order by GA
    o_rows =dist(average_matrix) %>%
                          hclust() %>%
                          as.dendrogram()  #%>%
                          #ladderize()
    o_col =dist(t(average_matrix),method="euclidean") %>%
                      hclust() %>%
                      as.dendrogram()  #%>%
                      #ladderize()
    
    cnv_col<-colorRamp2(c(0,0.5,1,1.5,2,3),c("#002C3E","#78BCC4","#F7F8F3","#F7444E","#aa1407","#440803"))

    side_ha_rna<-data.frame(ga_motif=markers_out[get_order(o_rows,1),]$RNA.auc)
    colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.90,0.95)),plasma(3))
    #colfun_rna=colorRamp2(c(0,1,2),plasma(3))

    side_ha_motif<-data.frame(chromvar_motif=markers_out[get_order(o_rows,1),]$chromvar.auc)
    colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.90,0.95)),cividis(3))
    #colfun_motif=colorRamp2(c(0,1,2),cividis(3))

    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-markers_out[markers_out$gene %in% row.names(tf_motif),]$chromvar.feature
    
    #plot into tmp_motif folder
    #note anno_bar reorders, so just supply in motif list order here
    system(paste0("rm -rf ",outdir,"/tmp_motifs"))
    system(paste0("mkdir -p ",outdir,"/tmp_motifs"))
    lapply(1:length(motif_list),function(i) {
      plt<-MotifPlot(
                    object = x,
                    assay="ATAC",
                    motifs = motif_list[i],ncol=1)+
                    theme_void()+
                    theme(strip.text = element_blank())
      if(nchar(i)==1){i<-paste0("0",i)}
      ggsave(plt,
            file=paste0(i,"_",prefix,".tf.heatmap.motif.png"),
            path=paste0(outdir,"/tmp_motifs/"),
            height=2,
            width=2,
            limitsize=F)
        })
    motif_plots<-list.files(
                path=paste0(outdir,"/tmp_motifs"),
                pattern="*motif.png",
                full.names=TRUE)

    side_ha_ga<-data.frame(ga_auc=markers_out[get_order(o_rows,1),]$GeneActivity.auc)
    colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.90,0.95)),magma(3))
    #colfun_ga=colorRamp2(c(0,1,2),magma(3))

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), 
                                            labels =row.names(tf_rna),
                                            labels_gp=gpar(fontsize=6)),
                            motifs = anno_image(motif_plots))

    o_col_split<-unlist(lapply(strsplit(colnames(tf_rna),"_"),'[',1))

    cnv_plot<-Heatmap(tf_cnv,
        cluster_rows = o_rows,
        name="CNV",
        col=cnv_col,
        column_title="CNV",
        column_names_gp = gpar(fontsize = 6),
        show_row_names=FALSE,
        column_names_rot=90,
        #column_split = o_col_split,
        cluster_column_slices = FALSE)

    rna_auc<-Heatmap(side_ha_rna,
        cluster_rows = o_rows,
        col=side_ha_col,
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=7))

    rna_plot<-Heatmap(tf_rna,
        cluster_rows = o_rows,
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 6),
        show_row_names=FALSE,
        column_names_rot=90,
        #column_split = o_col_split,
        cluster_column_slices = FALSE)

      rna_order<-draw(rna_plot)

      ga_auc<-Heatmap(side_ha_ga,
          cluster_rows = o_rows,         
          col=side_ha_col,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      ga_plot<-Heatmap(tf_ga,
          cluster_rows = o_rows,                 
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          column_names_gp = gpar(fontsize = 6),
          show_row_names=FALSE,
          column_names_rot=90,
        #column_split = o_col_split,
          column_order = unlist(column_order(rna_order)),
          cluster_column_slices = FALSE)

      motif_auc<-Heatmap(side_ha_motif,
          cluster_rows = o_rows,          
          col=side_ha_col,
          show_row_names=FALSE,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      motif_plot<-Heatmap(tf_motif,
          cluster_rows = o_rows,                 
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          #top_annotation=top_ha,
          column_names_gp = gpar(fontsize = 6),
          show_row_names=FALSE,
          column_names_rot=90,
          #column_split = o_col_split,
          column_order = unlist(column_order(rna_order)),
          cluster_column_slices = FALSE,
          right_annotation=gene_ha)
      
      #motif_image<-anno_image(paste0(prefix,".tf.heatmap.motif.svg"))
    
    pdf(paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf")),width=30,height=30)
    print(draw(cnv_plot+ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot,row_title=prefix))
    dev.off()
    print(paste("Plotted... ",paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf"))))
}

run_correlations<-function(x_data=tf_cnv,y_data=tf_rna,x_name="cnv",y_name="rna"){
  cor_out<-as.data.frame(do.call("rbind",
            lapply(1:nrow(x_data), function(i){
              if(sum(is.finite(x_data[i,]))>10 & sum(is.finite(y_data[i,]))>10){
                cor<-cor.test(c(as.numeric(x_data[i,])),c(as.numeric(y_data[i,])),method="spearman",use="pairwise.complete")
                return(
                  c(x_name,y_name,row.names(x_data)[i],row.names(y_data)[i],
                    median(as.numeric(x_data[i,]),na.rm=T),median(as.numeric(y_data[i,]),na.rm=T),
                    mean(as.numeric(x_data[i,]),na.rm=T),mean(as.numeric(y_data[i,]),na.rm=T),
                    var(as.numeric(x_data[i,]),na.rm=T),var(as.numeric(y_data[i,]),na.rm=T),
                    cor$statistic,as.numeric(cor$estimate),as.numeric(cor$p.value),-log10(as.numeric(cor$p.value))))}
              })))
  colnames(cor_out)<-c("x","y","x_rowname","y_rowname","x_median","y_median","x_mean","y_mean","x_var","y_var","cor_statistic","cor_rho","cor_pval","cor_neglog10_pval")
  cor_out <- cor_out %>% mutate(across(c(x_median, y_median, x_mean, y_mean,x_var,y_var,cor_rho,cor_pval,cor_neglog10_pval), as.numeric))
  return(cor_out)
}



correlate_tf_modalities<-function(x=out_subset,prefix,plot_by,outdir){
    gtf<-rtracklayer::readGFF(file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
    gtf<-makeGRangesFromDataFrame(gtf,keep.extra.columns=TRUE)
    gtf<-gtf[gtf$type=="gene",]
    gtf<-gtf[gtf@seqnames %in% paste0("chr",1:22),]

    #summarize markers over groups
    #grab all chromvar motif names
    tf_list<- ConvertMotifID(x, id = Features(x@assays$chromvar),assay="ATAC") #or ATAC as assay
    #clean up tf list names
    tf_list<-unlist(lapply(strsplit(tf_list,"\\("),"[",1))
    #shared motifs across assays
    shared_motif_idx<-Reduce(intersect, c(list(which(tf_list %in% Features(x@assays$RNA))), 
                                        list(which(tf_list %in% Features(x@assays$GeneActivity)))))
    #passing motifs
    tf_list_motifs<-Features(x@assays$chromvar)[shared_motif_idx]
    
    tf_list_genes<-tf_list[shared_motif_idx]
    length(tf_list_motifs)==length(tf_list_genes)

    print("Genes not overlapping all assays:")
    tf_list[!(tf_list %in% Features(x@assays$RNA))] #note that we are still excluding the combination motifs e.g. FOSL2::JUNB
    
    tf_rna<-average_features(x=x,features=tf_list_genes,assay="RNA",group_by=plot_by)

    tf_ga<-average_features(x=x,features=tf_list_genes,assay="GeneActivity",group_by=plot_by)

    tf_motif<-average_features(x=x,features=tf_list_motifs,assay="chromvar",group_by=plot_by)

    row.names(tf_motif)<- ConvertMotifID(x, id = tf_list_motifs,assay="ATAC") #or ATAC as assay
    row.names(tf_motif)<-unlist(lapply(strsplit(row.names(tf_motif),"\\("),"[",1))
    tf_motif<-tf_motif[row.names(tf_rna),]

    #assign cnv windows to genes
    cnv_granges<-makeGRangesFromDataFrame(data.frame(
        seqnames=unlist(lapply(strsplit(row.names(x@assays$cnv$data),"-"),"[",1)),
        start=unlist(lapply(strsplit(row.names(x@assays$cnv$data),"-"),"[",2)),
        end=unlist(lapply(strsplit(row.names(x@assays$cnv$data),"-"),"[",3))))
    gtf<-gtf[gtf$gene_name %in% row.names(tf_rna),]
    overlaps<-findOverlaps(gtf,cnv_granges)
    overlaps<-overlaps[!duplicated(overlaps@to),]
    gtf$cnv_windows<-NA
    names(gtf)<-1:length(gtf)
    gtf[overlaps@from,]$cnv_windows<-overlaps@to
    gtf<-gtf[!is.na(gtf$cnv_windows),]
    gtf<-gtf[!duplicated(gtf$cnv_windows),]

    tf_rna<-tf_rna[gtf$gene_name,]
    tf_motif<-tf_motif[gtf$gene_name,]
    tf_ga<-tf_ga[gtf$gene_name,]
    tf_cnv<-average_features(x=x,features=row.names(dat@assays$cnv@counts)[gtf$cnv_windows],assay="cnv",slot_name="counts",group_by=plot_by)

    #correlation plots between 
    #1. cnv and gene activity
    cnv_by_ga<-run_correlations(x_data=tf_cnv,y_data=tf_ga,x_name="cnv",y_name="ga")
    #2. cnv and rna
    cnv_by_rna<-run_correlations(x_data=tf_cnv,y_data=tf_rna,x_name="cnv",y_name="rna")
    #3. cnv and tf
    cnv_by_motif<-run_correlations(x_data=tf_cnv,y_data=tf_motif,x_name="cnv",y_name="motif")
    #4. ga and rna
    ga_by_rna<-run_correlations(x_data=tf_ga,y_data=tf_rna,x_name="ga",y_name="rna")
    #5. ga and tf
    ga_by_tf<-run_correlations(x_data=tf_ga,y_data=tf_motif,x_name="ga",y_name="motif")
    #6. rna and tf
    rna_by_tf<-run_correlations(x_data=tf_rna,y_data=tf_motif,x_name="rna",y_name="motif")

  plt1<-ggplot(cnv_by_ga,aes(x=x_median,y=y_median,color=as.numeric(cor_rho),size=as.numeric(cor_neglog10_pval)))+geom_point()+scale_color_gradient2(limits = c(-1, 1),low="blue",mid="white",high="red")+  scale_size(limits = c(0, 10))
  plt2<-ggplot(cnv_by_rna,aes(x=x_median,y=y_median,color=as.numeric(cor_rho),size=as.numeric(cor_neglog10_pval)))+geom_point()+scale_color_gradient2(limits = c(-1, 1),low="blue",mid="white",high="red")+  scale_size(limits = c(0, 10))
  plt3<-ggplot(cnv_by_motif,aes(x=x_median,y=y_median,color=as.numeric(cor_rho),size=as.numeric(cor_neglog10_pval)))+geom_point()+scale_color_gradient2(limits = c(-1, 1),low="blue",mid="white",high="red")+  scale_size(limits = c(0, 10))
  plt4<-ggplot(ga_by_rna,aes(x=x_median,y=y_median,color=as.numeric(cor_rho),size=as.numeric(cor_neglog10_pval)))+geom_point()+scale_color_gradient2(limits = c(-1, 1),low="blue",mid="white",high="red")+  scale_size(limits = c(0, 10))
  plt5<-ggplot(ga_by_tf,aes(x=x_median,y=y_median,color=as.numeric(cor_rho),size=as.numeric(cor_neglog10_pval)))+geom_point()+scale_color_gradient2(limits = c(-1, 1),low="blue",mid="white",high="red")+  scale_size(limits = c(0, 10))
  plt6<-ggplot(rna_by_tf,aes(x=x_median,y=y_median,color=as.numeric(cor_rho),size=as.numeric(cor_neglog10_pval)))+geom_point()+scale_color_gradient2(limits = c(-1, 1),low="blue",mid="white",high="red")+  scale_size(limits = c(0, 10))

  plt<-(plt1/plt2/plt3)|(plt4/plt5/plt6)+plot_layout(guides = "collect")
  ggsave(plt,file=paste0(outdir,"/",paste0(prefix,".tf.cor.pdf")),width=20,height=20)
    print(paste("Plotted... ",paste0(outdir,"/",paste0(prefix,".tf.cor.pdf"))))
}

clone_filter<-names(which(table(dat$merged_assay_clones)>=30))
dat_sub<-subset(dat,merged_assay_clones %in% clone_filter)
dat_sub<-subset(dat_sub,cells=row.names(dat_sub@meta.data)[!endsWith(dat_sub$merged_assay_clones,suffix="_normal")])
#dat_sub<-subset(dat_sub,Diag_MolDiag %in% c("IDC ER+/PR+/HER2-","IDC ER+/PR-/HER2-"))

plot_top_tf_markers(x=dat_sub,
                    group_by="merged_assay_clones",
                    plot_by="merged_assay_clones",
                    prefix="clone_tf",
                    n_markers=3,
                    order_by_idents=FALSE,
                    outdir=paste0(output_directory,"/cnvclones_by_diagnosis"))

#clonal barplots

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

#use only clones with at least 10 cell 
clone_filter<-names(which(table(dat$merged_assay_clones)>=10))

cellcycle_freq<-dat@meta.data %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,Phase,merged_assay_clones) %>% 
  count(Phase,.drop=TRUE)

plt<-ggplot(cellcycle_freq, aes(fill=factor(Phase,levels=names(Phase_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=Phase_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="phase_perclone_percbarplot.pdf",width=10,height=10)

#box plot of wu scores
emt_and_d_scores<-dat@meta.data %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,merged_assay_clones) 

plt<-ggplot(emt_and_d_scores, aes(y=Wu_EMTScores, x=merged_assay_clones, fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="emtscore_perclone_boxplot.pdf",width=10,height=10)


plt<-ggplot(emt_and_d_scores, aes(y=Wu_DScores, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="dscore_perclone_boxplot.pdf",width=10,height=10)

#check esr1 expression on clones
dat<-SCTransform(dat)

gene_expression<-FetchData(object = dat, vars = c("ESR1","KRT5","PGR"), assay="SCT")
chromatin_expression_esr1<-FetchData(object = dat, vars = c("MA0112.3"), assay="chromvar")
#chromatin_expression_pgr<-FetchData(object = dat, vars = c("MA2327.1"), assay="chromvar")
dat_tmp<-AddMetaData(dat,gene_expression)
dat_tmp<-AddMetaData(dat_tmp,chromatin_expression_esr1,col.name="ESR1_TF")
#dat_tmp<-AddMetaData(dat_tmp,chromatin_expression_pgr,col.name="PGR_TF")

esr1_clones<-dat_tmp@meta.data %>% 
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


plt3<-ggplot(esr1_clones, aes(y=PGR, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

plt4<-ggplot(esr1_clones, aes(y=ESR1_TF, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +theme(axis.text.x = element_text(angle = 90))+
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt1/plt2/plt3/plt4,file="ESR1_perclone_boxplot.pdf",width=10,height=10)


#per clone, 30 cells per clone minimum
scsubtype_freq<-dat@meta.data %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) %>% 
  count(scsubtype,.drop=TRUE)

plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")+ theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="scsubtype_perclone_percbarplot.pdf",width=10,height=10)
