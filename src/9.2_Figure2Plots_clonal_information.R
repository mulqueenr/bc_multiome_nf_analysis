```R
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
library(matrixStats)

library(msigdbr,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(fgsea,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(presto,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
setwd("/home/groups/MohammedLab/bc_multiome/seurat_object")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)


outdir="/home/groups/MohammedLab/bc_multiome/fig2"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}


clin_col=c("IDC ER+/PR-/HER2+"="#f37872", 
"DCIS DCIS"="#cccccb", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"NAT NA"="#c2d9ea")

write.table(dat@meta.data,row.names=T,col.names=T,sep="\t",file="cell_metadata.tsv")

#####Plot of heatmap of all clones ####

cnv_col<-c("0"="#002C3E", "0.5"="#78BCC4", "1"="#F7F8F3", "1.5"="#F7444E", "2"="#aa1407", "3"="#440803")
#colors from Curtis et al.

dat<-JoinLayers(dat,assay="RNA")
dat_sub<-subset(dat,cells=names(which(!is.na(dat$merged_assay_clones))))
dat_sub<-subset(dat_sub,cells=names(which(dat_sub$merged_assay_clones != "contamination")))
dat_sub<-subset(dat_sub,cells=names(which(dat_sub$assigned_celltype %in% c("cancer"))))
dat_sub<-subset(dat_sub,cells=row.names(dat_sub@meta.data)[which(!endsWith(dat_sub$merged_assay_clones,suffix="_normal"))])
table(dat_sub$merged_assay_clones,dat_sub$sample)

table(dat_sub$Diag_MolDiag,dat_sub$sample)

windows<-data.frame(chr=unlist(lapply(strsplit(row.names(dat_sub@assays$cnv@counts),"-"),"[",1)),
                    start=unlist(lapply(strsplit(row.names(dat_sub@assays$cnv@counts),"-"),"[",2)),
                    end=unlist(lapply(strsplit(row.names(dat_sub@assays$cnv@counts),"-"),"[",3)))
          
windows<-makeGRangesFromDataFrame(windows)

#relevant CNV genes from curtis work
#from https://www.nature.com/articles/s41416-024-02804-6#Sec20
#change RAB7L1 to RAB29
#lost RAB7L1

cnv_genes<-c('ESR1','PGR','DLEU2L', 'TRIM46', 'FASLG', 'KDM5B', 'RAB7L1', 'PFN2', 'PIK3CA', 'EREG', 'AIM1', 'EGFR', 'ZNF703', 'MYC', 'SEPHS1', 'ZMIZ1', 'EHF', 'POLD4', 'CCND1', 'P2RY2', 'NDUFC2-KCTD14', 'FOXM1', 'MDM2', 'STOML3', 'NEMF', 'IGF1R', 'TP53I13', 'ERBB2', 'SGCA', 'RPS6KB1', 'BIRC5', 'NOTCH3', 'CCNE1', 'RCN3', 'SEMG1', 'ZNF217', 'TPD52L2', 'PCNT', 'CDKN2AIP', 'LZTS1', 'PPP2R2A', 'CDKN2A', 'PTEN', 'RB1', 'CAPN3', 'CDH1', 'MAP2K4', 'GJC2', 'TERT', 'RAD21', 'ST3GAL1', 'SOCS1')
cnv_genes_class<-c('amp','amp','amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'amp', 'amp', 'amp', 'amp', 'amp')
cnv_genes<-setNames(cnv_genes_class,cnv_genes)
cnv_genes<-cnv_genes[names(cnv_genes) %in% dat_sub@assays$ATAC@annotation$gene_name]
cnv_genes_windows<-dat_sub@assays$ATAC@annotation[dat_sub@assays$ATAC@annotation$gene_name %in% names(cnv_genes),] #filter annotation to genes we want
cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
windows<-findOverlaps(windows,cnv_genes_windows)

#filter to genes that actually show changes

annot<-data.frame(
  window_loc=queryHits(windows),
  gene=cnv_genes_windows$gene_name,
  cnv_class=unname(cnv_genes[cnv_genes_windows$gene_name]))

annot$col<-ifelse(annot$cnv_class=="amp","red","blue")

cell_cnv<-t(as.data.frame(dat_sub@assays$cnv@counts))

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
output_directory=paste0(dirname(getwd()),"/cnv_calls")

row_title_color<-unique(dat_sub@meta.data[row.names(dat_sub@meta.data) %in% colnames(dat_sub@assays$cnv@counts),]$merged_assay_clones)
title_color<-data.frame(Diag_MolDiag=dat_sub@meta.data[!duplicated(dat_sub@meta.data$merged_assay_clones),]$Diag_MolDiag,merged_assay_clones=dat_sub@meta.data[!duplicated(dat_sub@meta.data$merged_assay_clones),]$merged_assay_clones)

title_color$col<-clin_col[title_color$Diag_MolDiag]

pdf(paste0(outdir,"/","FIG2_all_samples.cnv.heatmap.pdf"),height=90,width=40)
Heatmap(cell_cnv,
  col=cnv_col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = FALSE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  row_title_gp = gpar(col = title_color$col),
  bottom_annotation=hc,
  row_split=dat_sub@meta.data[row.names(dat_sub@meta.data) %in% colnames(dat_sub@assays$cnv@counts),]$merged_assay_clones,
  column_split=factor(unlist(lapply(strsplit(row.names(dat_sub@assays$cnv@counts),"-"),"[",1)),levels=paste0("chr",1:22)),
  border = TRUE)
dev.off()
print(paste0(outdir,"/","FIG2_all_samples.cnv.heatmap.pdf"))

```

####################################################
#           Fig Clone Heatmap                      #
###################################################

```R
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
plot_top_tf_markers<-function(x=out_subset,colfun,group_by,prefix,n_markers=20,order_by_idents=FALSE,plot_by,outdir){
    #read in track to assign CNV window locations
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
                        n=n_markers))) #grab top N TF markers per celltype
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
    row.names(tf_cnv)<-row.names(tf_rna)
    average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering
    #average_matrix=tf_motif #just cluster only on cnvs

    #set up heatmap seriation and order by average z score
    #first cluster columns, then sort row orders by average

    o_cols =t(average_matrix) %>% dist(method="maximum")  %>% 
                          hclust(method="ward.D2") %>%
                          as.dendrogram() %>%
                          ladderize(decreasing=FALSE) %>% labels()

    o_rows = average_matrix %>% dist(method="maximum")  %>% 
                          hclust(method="ward.D2") %>%
                          as.dendrogram() %>%
                          ladderize(decreasing=FALSE) %>% labels()

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

    #colfun_ga=colorRamp2(c(-2,0,2),c("#053061","#ffffff","#67001f"))
    #colfun_motif=colorRamp2(c(-2,0,2),c("#4d4d4d","#ffffff","#e08214"))
    #colfun_rna=colorRamp2(c(-2,0,2),c("#313695","#ffffff","#a50026"))
    cnv_col<-colorRamp2(c(0,0.5,1,1.5,2,3),c("#002C3E","#78BCC4","#F7F8F3","#F7444E","#aa1407","#440803"))
    colfun_ga<-colfun
    colfun_rna<-colfun
    colfun_motif<-colfun

    tf_rna<-tf_rna[o_rows,o_cols]
    tf_cnv<-tf_cnv[o_rows,o_cols]
    tf_ga<-tf_ga[o_rows,o_cols]
    tf_motif<-tf_motif[o_rows,o_cols]

    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), 
                                            labels =row.names(tf_rna),
                                            labels_gp=gpar(fontsize=6)),
                            motifs = anno_image(motif_plots))

    cnv_plot<-Heatmap(tf_cnv,
        cluster_rows = FALSE,
        cluster_columns= FALSE,
        name="CNV",
        col=cnv_col,
        column_title="CNV",
        column_names_gp = gpar(fontsize = 6),
        show_row_names=FALSE,
        column_names_rot=90,
        cluster_column_slices = FALSE)

    rna_plot<-Heatmap(tf_rna,
        cluster_rows = FALSE,
        cluster_columns= FALSE,
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 6),
        show_row_names=FALSE,
        column_names_rot=90,
        cluster_column_slices = FALSE)
    
      ga_plot<-Heatmap(tf_ga,
        cluster_rows = FALSE,
        cluster_columns= FALSE,
          column_title="Gene Activity",
          col=colfun_ga,
          column_names_gp = gpar(fontsize = 6),
          show_row_names=FALSE,
          column_names_rot=90,
          cluster_column_slices = FALSE)


      motif_plot<-Heatmap(tf_motif,
        cluster_rows = FALSE,
        cluster_columns= FALSE,
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          column_names_gp = gpar(fontsize = 6),
          show_row_names=FALSE,
          column_names_rot=90,
          cluster_column_slices = FALSE,
          right_annotation=gene_ha)
          
    pdf(paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf")),width=30,height=30)
    print(draw(cnv_plot+ga_plot+rna_plot+motif_plot,row_title=prefix))
    dev.off()
    print(paste("Plotted... ",paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf"))))
}


clone_filter<-names(which(table(dat_sub$merged_assay_clones)>=30))
dat_sub<-subset(dat_sub,cells=row.names(dat_sub@meta.data)[dat_sub$merged_assay_clones %in% clone_filter]) #limit to cells passing cnv
table(dat_sub$merged_assay_clones,dat_sub$sample)

colfun <- colorRamp2(
  #breaks = c(-3,-1,-0.5,0,0.5,1,3),
  breaks = c(-3,-2,-1,0,1,2,3),
  colors = c("#053061","#487590","#d1e5f0","#f7f7f7","#fddbc7","#a76146","#67001f"))

plot_top_tf_markers(x=dat_sub,
                    group_by="Diag_MolDiag", #groups to find marker genes, Diag_MolDiag?
                    plot_by="merged_assay_clones", #groups to plot by
                    prefix="FIG2_clone_tf",
                    colfun=colfun,
                    n_markers=5,
                    order_by_idents=FALSE,
                    outdir=outdir)

```

Plotting correlations across all genes
```R

#read in track to assign CNV window locations
gtf<-rtracklayer::readGFF(file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
gtf<-makeGRangesFromDataFrame(gtf,keep.extra.columns=TRUE)
gtf<-gtf[gtf$type=="gene",]
gtf<-gtf[gtf@seqnames %in% paste0("chr",1:22),]

#assign cnv windows to genes
cnv_granges<-makeGRangesFromDataFrame(data.frame(
    seqnames=unlist(lapply(strsplit(row.names(dat_sub@assays$cnv$data),"-"),"[",1)),
    start=unlist(lapply(strsplit(row.names(dat_sub@assays$cnv$data),"-"),"[",2)),
    end=unlist(lapply(strsplit(row.names(dat_sub@assays$cnv$data),"-"),"[",3))))

overlaps<-findOverlaps(gtf,cnv_granges)
gtf$cnv_windows<-NA
gtf[overlaps@from,]$cnv_windows<-overlaps@to
cnv_mat=GetAssayData(dat_sub,assay="cnv",layer="data")
cnv_mat=cnv_mat[gtf$cnv_windows,]
row.names(cnv_mat)<-gtf$gene_name

#correlations across modalities
ga_mat=t(scale(t(GetAssayData(dat_sub,assay="GeneActivity",layer="data"))))
ga_mat<-ga_mat[is.finite(rowSums(ga_mat)) & rowVars(ga_mat) > 0, ]

dat_sub<-JoinLayers(dat_sub,assay="RNA")
rna_mat=t(scale(t(LayerData(dat_sub,assay="RNA",layer="data"))))
rna_mat<-rna_mat[is.finite(rowSums(rna_mat)) & rowVars(rna_mat) > 0, ]

tf_mat=LayerData(dat_sub,assay="chromvar",layer="data") #pre scaled
row.names(tf_mat) <- ConvertMotifID(object = dat_sub, assay="ATAC",id = row.names(tf_mat))

pairwise_correlations<-function(x_name="chromvar",y_name="cnv",x_mat=tf_mat,y_mat=cnv_mat,cores=1){
  common_genes<-intersect(row.names(x_mat),row.names(y_mat))
  x_mat<-x_mat[row.names(x_mat) %in% common_genes,]
  y_mat<-y_mat[row.names(y_mat) %in% common_genes,]

  x_y_cor_mat<-mclapply(1:nrow(x_mat), function(x){
    cor_out<-cor.test(x_mat[x,],y_mat[x,],method="spearman",)
    return(c(row.names(x_mat)[x],
    mean(x_mat[x,]),var(x_mat[x,]),
    mean(y_mat[x,]),var(y_mat[x,]),
    cor_out$p.value,
    cor_out$estimate))
  },mc.cores=cores)

  x_y_cor_mat<-as.data.frame(do.call("rbind",x_y_cor_mat)) 
  x_y_cor_mat$V8 <- p.adjust(x_y_cor_mat$V6, method = "bonferroni")
  colnames(x_y_cor_mat)<-c("gene",
  paste0(x_name,"_","mean"),paste0(x_name,"_","var"),
  paste0(y_name,"_","mean"),paste0(y_name,"_","var"),
  paste0(x_name,"_",y_name,"_pval"),paste0(x_name,"_",y_name,"_spearmanrho"),
  paste0(x_name,"_",y_name,"_qval"))
  saveRDS(x_y_cor_mat,file=paste0(outdir,"/",x_name,"_",y_name,".spearman_corr.rds"))
}

#cnv correlations
pairwise_correlations(x_name="tf",y_name="cnv",x_mat=tf_mat,y_mat=cnv_mat)
pairwise_correlations(x_name="ga",y_name="cnv",x_mat=ga_mat,y_mat=cnv_mat,cores=4)
pairwise_correlations(x_name="rna",y_name="cnv",x_mat=rna_mat,y_mat=cnv_mat,cores=4)

#modality correlations
pairwise_correlations(x_name="ga",y_name="rna",x_mat=ga_mat,y_mat=rna_mat,cores=4)
pairwise_correlations(x_name="ga",y_name="tf",x_mat=ga_mat,y_mat=tf_mat)
pairwise_correlations(x_name="rna",y_name="tf",x_mat=rna_mat,y_mat=tf_mat)

corr_files<-list.files(outdir,pattern=".spearman_corr.rds",full.names=T)
corr_list<-lapply(corr_files,readRDS)

# Merge all data frames together by the common column
merged_df <- corr_list %>% purrr::reduce(dplyr::full_join, by = "gene")
row.names(merged_df)<-merged_df$gene
merged_df <- merged_df %>%
  mutate(across(where(is.character), as.numeric))
merged_df$gene<-row.names(merged_df)

#plot dotplot of correlations, scale by rho, color if significant
#label top 10 per quadrant that are significant
#color each quadrant different

# Assuming your data frame is called 'plot_data'
#quadrant1 and label1 are for plot1
#setting point colors by sig
#setting labels by quadrant top hits (distance from center)

merged_df <- merged_df %>%
  mutate(
    col1= case_when(
      rna_cnv_qval>=0.05 & ga_cnv_qval>=0.05 ~ "#333333",
      rna_cnv_qval<0.05 & ga_cnv_qval>=0.05 ~ "#FF0000",
      rna_cnv_qval>=0.05 & ga_cnv_qval<0.05 ~ "#000080",
      rna_cnv_qval<0.05 & ga_cnv_qval<0.05 ~ "#800080"),
    quadrant1 = case_when(
      rna_cnv_spearmanrho > 0 & ga_cnv_spearmanrho > 0 ~ "Q1_TopRight",
      rna_cnv_spearmanrho < 0 & ga_cnv_spearmanrho > 0 ~ "Q2_TopLeft",
      rna_cnv_spearmanrho < 0 & ga_cnv_spearmanrho < 0 ~ "Q3_BottomLeft",
      rna_cnv_spearmanrho > 0 & ga_cnv_spearmanrho < 0 ~ "Q4_BottomRight"),
    distance1 = sqrt(rna_cnv_spearmanrho^2 + ga_cnv_spearmanrho^2)) %>%
  group_by(quadrant1) %>%
  mutate(rank1 = rank(-distance1, ties.method = "first"),
    label1 = ifelse(rank1 <= 10, gene, "") # Only keep top 10
  ) %>% ungroup()


#do same for tf motifs

merged_df <- merged_df %>%
  mutate(
    col2= case_when(
      rna_tf_qval>=0.05 & tf_cnv_qval>=0.05 ~ "#333333",
      rna_tf_qval<0.05 & tf_cnv_qval>=0.05 ~ "#FF0000",
      rna_tf_qval>=0.05 & tf_cnv_qval<0.05 ~ "#000080",
      rna_tf_qval<0.05 & tf_cnv_qval<0.05 ~ "#800080"),
    quadrant2 = case_when(
      rna_tf_spearmanrho > 0 & tf_cnv_spearmanrho > 0 ~ "Q1_TopRight",
      rna_tf_spearmanrho < 0 & tf_cnv_spearmanrho > 0 ~ "Q2_TopLeft",
      rna_tf_spearmanrho < 0 & tf_cnv_spearmanrho < 0 ~ "Q3_BottomLeft",
      rna_tf_spearmanrho > 0 & tf_cnv_spearmanrho < 0 ~ "Q4_BottomRight"),
    distance2 = sqrt(rna_tf_spearmanrho^2 + tf_cnv_spearmanrho^2)) %>%
  group_by(quadrant2) %>%
  mutate(rank2 = rank(-distance2, ties.method = "first"),
    label2 = ifelse(rank2 <= 10, gene, "") # Only keep top 10
  ) %>% ungroup()


plt1<-ggplot(merged_df,aes(x=rna_cnv_spearmanrho,y=ga_cnv_spearmanrho,color=col1,label=label1))+
      geom_point(size=0.5,alpha=1)+
      scale_color_identity()+
      theme_minimal()+
      geom_text_repel(size=1,max.overlaps=50,min.segment.length=0)


plt2<-ggplot(merged_df,
  aes(x=rna_tf_spearmanrho,y=tf_cnv_spearmanrho,color=col2,label=label2))+
  geom_point(size=0.5,alpha=1)+
  scale_color_identity()+
  theme_minimal()+
  geom_text_repel(size=1,max.overlaps=50,min.segment.length=0)


ggsave(plt1/plt2,file=paste0(outdir,"/","spearman_corr.pdf"),height=15,width=10)
print(paste0(outdir,"/","spearman_corr.pdf"))

saveRDS(merged_df,file="gene_modality_correlations.rds")
