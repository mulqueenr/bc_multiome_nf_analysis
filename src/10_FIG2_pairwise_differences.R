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
library(grImport,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(BSgenome.Hsapiens.UCSC.hg38)
library(dendextend)
library(ggdendro)
library(circlize)
library(ggtern)
library(GeneNMF)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
dat<-subset(dat,cells=row.names(dat@meta.data)[isNA(dat@meta.data$merged_assay_clones) | dat@meta.data$merged_assay_clones != "contamination"])
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])
hist_col=c(
  "IDC"="#FF9966",
  "ILC"="#006633")

clin_col=c(
  "ER+/PR+/HER2-"="#8d86c0", 
  "ER+/PR-/HER2-"="#7fd0df")

scsubtype_col=c(
  "SC_Subtype_LumA_SC"="#2b2c76",
  "SC_Subtype_LumB_SC"="#86cada")

dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
DefaultAssay(dat)<-"ATAC"
dat <- RegionStats(dat, genome = BSgenome.Hsapiens.UCSC.hg38)

#### Pairwise comparisons
system(paste0("mkdir -p ",paste0(dirname(getwd()),"/pairwise_comparisons"))) #paste pairwise comparisons into directory one folder up
output_directory=paste0(dirname(getwd()),"/pairwise_comparisons")

#tornado plot of top DA peaks
tornado_plot<-function(obj=obj,da_peak_set=markers,i="IDC",peak_count=500,col=col,col_lim=0.03){
    print(i)
    top_peaks <- da_peak_set %>% 
                  filter(group==i) %>%  
                  filter(logFC>0) %>% 
                  slice_max(logFC,n=peak_count)
    print(head(top_peaks))

    obj_mat<-RegionMatrix(obj,key="DA_mat",
    regions=StringToGRanges(top_peaks$feature),
    upstream=5000,downstream=5000,
    assay="ATAC")

    plt<-RegionHeatmap(obj_mat,key="DA_mat",
      upstream=5000,downstream=5000,
      order=TRUE, 
      window=(10000)/100, normalize=TRUE,
      assay="ATAC", 
      idents=levels(Idents(obj)),
      cols=col[i],max.cutoff=col_lim,
      nrow=length(unique(da_peak_set$group)))+ 
    ggtitle(i)

    print("Returning plot...")
    return(plt)
}

#volcano plot of top DE genes and chromvar TFs
volcano_plot<-function(obj=obj,de_features_set=markers,prefix,feature_count=10,outname,assay,group1,group2,col,outdir){
  
  feature_markers_group1<-de_features_set %>% filter(padj<0.05) %>% slice_max(n=feature_count,order_by=logFC)
  feature_markers_group2<-de_features_set %>% filter(padj<0.05) %>% slice_min(n=feature_count,order_by=logFC)
  group_1_count=de_features_set %>% filter(padj<0.05) %>% filter(logFC>0) %>% nrow()
  group_2_count=de_features_set %>% filter(padj<0.05) %>% filter(logFC<0) %>% nrow()

  de_features_set$fill_col<-"#808080"
  de_features_set[de_features_set$padj<0.05 & de_features_set$logFC>0,]$fill_col<-col[group1]
  de_features_set[de_features_set$padj<0.05 & de_features_set$logFC<0,]$fill_col<-col[group2]

  de_features_set$label<-NA
  de_features_set[de_features_set$feature %in% feature_markers_group1$feature,]$label<-feature_markers_group1$feature
  de_features_set[de_features_set$feature %in% feature_markers_group2$feature,]$label<-feature_markers_group2$feature

  plt<-ggplot(de_features_set,
              aes(x=logFC,
                  y=-log10(padj),
                  color=fill_col,
                  label=label))+
        geom_point(size=0.5,aes(alpha=0.5,stroke=0))+geom_hline(yintercept=-log10(0.05))+
        scale_color_identity()+ coord_cartesian(clip = "off") + theme_minimal() +
        ggrepel::geom_label_repel(size=1,box.padding=0.1,label.padding = 0.1,max.overlaps=Inf,min.segment.length = 0,vjust = "inward")+
        ggtitle(paste(group1,":",as.character(group_1_count),
        "\n",group2,":",as.character(group_2_count)))+ theme(legend.position="none")

  ggsave(plt,
          file=paste("pairwise",outname,"volcano",assay,"pdf",sep="."),
          path=outdir,
          units="in",width=2,height=2)
}

#gsea of top DE genes
gsea_enrichment<-function(annot,species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          outname=outname,
                          de_features_set,
                          col,group1,group2,
                          assay){
  pathwaysDF <- msigdbr(species=species, 
                        category=category, 
                        subcategory = subcategory)

  #limit pathways to genes in our data
  pathwaysDF<-pathwaysDF[pathwaysDF$ensembl_gene %in% unique(annot[annot$gene_biotype=="protein_coding",]$gene_id),]
  
  pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

  group1_features<-de_features_set %>%
    dplyr::filter(group == group1) %>%
    #dplyr::filter(padj<0.05) %>% 
    dplyr::arrange(logFC) %>%
    dplyr::select(feature, logFC)

  ranks<-setNames(nm=group1_features$feature,group1_features$logFC)

  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 10,
                    nproc = 1)

  topPathwaysUp <- fgseaRes %>% filter(ES > 0) %>% slice_max(NES,n=10) %>% dplyr::select(pathway)
  topPathwaysDown <- fgseaRes %>% filter(ES < 0) %>% slice_max(abs(NES),n=10) %>% dplyr::select(pathway)
  topPathways <- unlist(c(topPathwaysUp, rev(topPathwaysDown)))

  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))
  #plot gsea ranking of genes by top pathways
  #pdf(paste0("pairwise.",outname,".",assay,".",out_setname,".gsea.pdf"),width=20,height=10)
  #print(plt)
  #dev.off()

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% slice_max(NES,n= 10),
    fgseaRes  %>% slice_min(NES,n= 10))
  
  nes_plt_dat$col<-"#808080"
  nes_plt_dat[nes_plt_dat$NES>0 & nes_plt_dat$pval<0.05,]$col<-col[group1]
  nes_plt_dat[nes_plt_dat$NES<0 & nes_plt_dat$pval<0.05,]$col<-col[group2]

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= col)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(out_setname)+ylim(c(-4,4))
  return(patchwork::wrap_plots(list(plt1,plt2),ncol=2))
}

plot_gsea<-function(obj,annot,dmrs,
                    outname=outname,
                    assay=assay,col=col,group1,group2,outdir){

  #run gsea enrichment on different sets
  tft_plt<-gsea_enrichment(species="human",
              category="C3",
              subcategory="TFT:GTRD",
              out_setname="TFT",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  position_plt<-gsea_enrichment(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  hallmark_plt<-gsea_enrichment(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)
  plt<-patchwork::wrap_plots(list(tft_plt,position_plt,hallmark_plt),nrow=3,axes="collect_x")
  ggsave(plt,
          file=paste("pairwise",outname,"GSEA",assay,"pdf",sep="."),
          path=outdir,
          width=20,height=10)


}

cov_plot_per_gene<-function(obj,col,i,group1,group2,obj_group1,obj_group2){
  annot_plot<-AnnotationPlot(object=obj,region=i)

  cov_plot <- CoveragePlot(
    object = obj,
    region = i,
    annotation = FALSE,
    peaks = TRUE,links=FALSE)+
    scale_fill_manual(values=col)

  link_plot_1 <- LinkPlot(
    object = obj_group1,
    region = i)+
    scale_color_gradient2(limits=c(0,0.3),low="white",high=col[group1])

  link_plot_2 <-LinkPlot(
    object = obj_group2,
    region = i)+
    scale_color_gradient2(limits=c(0,0.3),low="white",high=col[group2])

  expr_plot <- ExpressionPlot(
    object = obj,
    features = i,
    assay = "SCT") + scale_fill_manual(values=col)

  plt<-CombineTracks(
    plotlist = list(cov_plot, annot_plot, link_plot_1,link_plot_2),
    expression.plot = expr_plot,
    heights = c(10, 2, 3, 3),
    widths = c(10, 3))
  return(plt)
}

#coverage of top GA+RNA genes
coverage_plot<-function(obj,markers_rna,markers_ga,col,outname,group1,group2,group_by,outdir){
  DefaultAssay(obj)<-"ATAC"

  da_combined<-merge(markers_rna,markers_ga,by="feature")
  da_combined$avg_logFC<-rowMeans(da_combined[,c('logFC.x', 'logFC.y')], na.rm=TRUE) #dont actually need this
  da_combined$avg_AUC<-rowMeans(da_combined[,c('auc.x', 'auc.y')], na.rm=TRUE) #dont actually need this

  group1_enriched<-da_combined %>% filter(padj.x<0.05) %>% filter(padj.y<0.05) %>% filter(logFC.x>0 & logFC.y>0) %>% slice_max(avg_AUC,n=10)
  group2_enriched<-da_combined %>% filter(padj.x<0.05) %>% filter(padj.y<0.05) %>% filter(logFC.x<0 & logFC.y<0) %>% slice_max(avg_AUC,n=10)
  genes=c(group1_enriched$feature,group2_enriched$feature)
  
  obj_group1<-subset(obj,cells=row.names(obj@meta.data)[obj@meta.data[,group_by] %in% c(group1)])
  obj_group2<-subset(obj,cells=row.names(obj@meta.data)[obj@meta.data[,group_by] %in% c(group2)])

  # link peaks to genes
  obj_group1<- LinkPeaks(
    object = obj_group1,
    peak.assay = "ATAC",
    expression.assay = "SCT",
    genes.use = genes)
  
  obj_group2<- LinkPeaks(
    object = obj_group2,
    peak.assay = "ATAC",
    expression.assay = "SCT",
    genes.use = genes)

  plt_list<-lapply(genes,function(i) {cov_plot_per_gene(i=i,obj=obj,col=col,group1,group2,obj_group1,obj_group2)})
  plt<-wrap_plots(plt_list,nrow=2,ncol=10)
  ggsave(plt,
          file=paste("pairwise",outname,"coverage","pdf",sep="."),
          path=outdir,
          width=50,height=20,limitsize=FALSE)

}

#wrapper of pairwise functions
pairwise_comparison<-function(obj=dat_cancer,
                              group_by="Diagnosis",
                              group1="IDC",
                              group2="ILC",
                              outname="diagnosis",
                              motif_name="ESR1",
                              downsample_cells_per_sample=50,
                              plot_tornado=TRUE,
                              col,
                              outdir){

  #make output directory
  system(paste0("mkdir -p ",outdir))
  #subset to relevent groups
  obj <-subset(obj, cells=row.names(obj@meta.data)[obj@meta.data[,group_by] %in% c(group1,group2)])
  DefaultAssay(obj)<-"ATAC"
  annot<-Annotation(obj)

  #downsample to ~ equal cell counts per sample
  obj_full<-obj
  Idents(obj)<-obj$sample
  downsample_cell_table=table(Idents(obj))
  obj <-subset(obj,downsample=downsample_cells_per_sample)
  Idents(obj)<-obj@meta.data[,group_by]

  ##########peaks and tornado plots##############
  print("Generating tornado plots...")
  assay="ATAC"

  #all peaks
  markers <- presto:::wilcoxauc.Seurat(
    X = obj, 
    group_by = group_by, 
    groups_use=c(group1,group2),
    y=c(group1,group2), 
    seurat_assay = assay)
  
  markers<- markers %>% filter(padj<0.05)
  write.table(markers,
              col.names=T,
              row.names=F,
              file=paste0(outdir,"/",paste("pairwise.all_da_peaks",outname,"tsv",sep=".")),
              sep="\t")

  if(plot_tornado){
  plt_list<-lapply(unique(markers$group),function(j) {
      tornado_plot(obj=obj,da_peak_set=markers,i=j,col=col)})
  plt<-wrap_plots(plt_list,nrow=1,guides='collect')
  ggsave(plt,
        file=paste("pairwise",outname,"tornado",assay,"all_da_peaks","pdf",sep="."),
        path=outdir,
        width=20,height=20)
  }

  # #ESR1 only
  # #da peaks that overlap with ESR1 motif only
  # motif<-names(obj@assays$ATAC@motifs@motif.names[which(obj@assays$ATAC@motifs@motif.names==motif_name)])
  # da_peaks_motif_filt<-markers[markers$feature %in% names(which(obj@assays$ATAC@motifs@data[,motif])),]

  # if(plot_tornado){
  # plt_list<-lapply(unique(da_peaks_motif_filt$group),function(j) {
  #     tornado_plot(obj=obj,da_peak_set=da_peaks_motif_filt,i=j,col=col)
  #     })
  # plt<-wrap_plots(plt_list,nrow=1,guides='collect')
  # ggsave(plt,file=paste("pairwise",outname,"tornado",assay,"ESR1_peaks","pdf",sep="."),width=20,height=20)
  # }

  ##########SCT plots##############
  print("Running SCT pairwise comparisons...")
  assay="SCT"
  obj<-SCTransform(obj,vars.to.regress = "nCount_RNA")
  markers_rna <- presto:::wilcoxauc.Seurat(X = obj, 
                                     group_by = group_by,
                                     groups_use=c(group1,group2),
                                     seurat_assay = assay,
                                     assay="data",
                                     y=group1)

  markers_rna<- markers_rna %>% filter(group==group1)
 
  write.table(markers_rna,
              col.names=T,
              row.names=F,
              file=paste0(outdir,"/",paste("pairwise",outname,assay,"tsv",sep=".")),
              sep="\t")

  # #all genes
  # volcano_plot(obj=obj,
  #             de_features_set=markers_rna,
  #             outname=paste0(outname,".allgenes"),
  #             assay=assay,group1=group1,group2=group2,col=col)

  # plot_gsea(obj=obj,
  #           dmrs=markers_rna,
  #           outname=paste0(outname,".allgenes"),
  #           assay=assay,group1=group1,group2=group2,col=col,annot=annot)
  
  #protein coding only
  print("Running SCT pairwise comparisons (protein coding only)...")
  markers_rna<-markers_rna[markers_rna$feature %in% annot[annot$gene_biotype=="protein_coding",]$gene_name,]
  volcano_plot(obj=obj,
              de_features_set=markers_rna,
              outname=paste0(outname,".proteincoding"),
              assay=assay,group1=group1,group2=group2,col=col,outdir=outdir)

  plot_gsea(obj=obj,
            dmrs=markers_rna,
            outname=paste0(outname,".proteincoding"),
            assay=assay,group1=group1,group2=group2,col=col,annot=annot,outdir=outdir)

  # ##########GENE ACTIVITY plots##############
  print("Running Gene Activity pairwise comparisons...")
  assay="GeneActivity"

  markers_ga <- presto:::wilcoxauc.Seurat(X = obj, 
                                      group_by = group_by,
                                      groups_use=c(group1,group2),
                                      seurat_assay = assay,
                                      assay="data")
   markers_ga<- markers_ga %>% filter(group==group1)
  write.table(markers_ga,
              col.names=T,
              row.names=F,
              file=paste0(outdir,"/",paste("pairwise",outname,assay,"tsv",sep=".")),
              sep="\t")

  # #  volcano_plot(obj=obj,
  # #              de_features_set=markers_ga,
  # #              outname=outname,
  # #              assay=assay,group1=group1,group2=group2,col=col)

  # #  plot_gsea(obj=obj,
  # #            dmrs=markers_ga,
  # #            outname=outname,
  # #            assay=assay,group1=group1,group2=group2,col=col,annot=annot)

  print("Running Gene Activity pairwise comparisons (protein coding only)...")
  markers_ga<-markers_ga[markers_ga$feature %in% annot[annot$gene_biotype=="protein_coding",]$gene_name,]

  volcano_plot(obj=obj,
              de_features_set=markers_ga,
              outname=paste0(outname,".proteincoding"),
              assay=assay,group1=group1,group2=group2,col=col,outdir=outdir)

  plot_gsea(obj=obj,
            dmrs=markers_ga,
            outname=paste0(outname,".proteincoding"),
            assay=assay,group1=group1,group2=group2,col=col,annot=annot,outdir=outdir)

  ########Coverage plots of GA and RNA################
  coverage_plot(obj=obj,
                markers_rna=markers_rna,
                markers_ga=markers_ga,
                col=col,
                outname=outname,
                group1=group1,
                group2=group2,
                group_by=group_by,
                outdir=outdir)

  ##########CHROMVAR plots##############
  print("Running CHROMVAR pairwise comparisons...")
  assay="chromvar"
  markers_tf  <- FindMarkers(
                            object = obj,
                            ident.1 = group1,
                            ident.2 = group2,
                            only.pos = FALSE,
                            mean.fxn = rowMeans,
                            fc.name = "avg_diff",
                            assay=assay,
                            test="LR")
  markers_tf$group<-group1
  markers_tf$logFC<-markers_tf$avg_diff
  markers_tf$padj<-markers_tf$p_val
  markers_tf$pval<-markers_tf$p_val_adj
  markers_tf$feature<-row.names(markers_tf)
  markers_tf$feature<-ConvertMotifID(object=obj,
                                  assay="ATAC",
                                  id=markers_tf$feature)
  write.table(markers_tf,
              col.names=T,
              row.names=F,
              file=paste0(outdir,"/",paste("pairwise",outname,assay,"tsv",sep=".")),
              sep="\t")

  volcano_plot(obj=obj,
              de_features_set=markers_tf,
              outname=outname,
              assay=assay,group1=group1,group2=group2,col=col,outdir=outdir)

}

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

#modified to plot by sample and group by pairwise (also added a column annotation)
plot_top_tf_markers<-function(x=out_subset,group_by,prefix,n_markers=20,order_by_idents=TRUE,plot_by,outdir){
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
    #summarize markers over samples (rather than by groups)
    tf_rna<-average_features(x=x,features=markers_out$gene,assay="RNA",group_by=plot_by)
    tf_rna<-tf_rna[row.names(tf_rna) %in% markers_out$gene,]
    tf_ga<-average_features(x=x,features=markers_out$gene,assay="GeneActivity",group_by=plot_by)
    tf_ga<-tf_ga[row.names(tf_ga) %in% markers_out$gene,]
    tf_motif<-average_features(x=x,features=markers_out$chromvar.feature,assay="chromvar",group_by=plot_by)
    tf_motif<-tf_motif[row.names(tf_motif) %in% markers_out$chromvar.feature,]
    row.names(tf_motif)<-markers_out[markers_out$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]
    #average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering
    average_matrix=tf_motif #just cluster only on TF motifs

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
    #colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.90,0.95)),magma(3))
    colfun_ga=colorRamp2(c(0,1,2),magma(3))

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), 
                                            labels =row.names(tf_rna),
                                            labels_gp=gpar(fontsize=6)),
                            motifs = anno_image(motif_plots))

    o_col_split<-unlist(lapply(strsplit(colnames(tf_rna)," "),'[',1))

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
        column_split = o_col_split,
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
        column_split = o_col_split,
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
        column_split = o_col_split,
          column_order = unlist(column_order(rna_order)),
          cluster_column_slices = FALSE,
          right_annotation=gene_ha)
      
      #motif_image<-anno_image(paste0(prefix,".tf.heatmap.motif.svg"))
    
    pdf(paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf")))
    print(draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot,row_title=prefix))
    dev.off()
    print(paste("Plotted... ",paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf"))))
}


#modified to plot by sample and group by pairwise (also added a column annotation)
plot_top_tf_markers_tfonly<-function(x=out_subset,group_by,prefix,n_markers=20,order_by_idents=TRUE,plot_by,outdir){
    #define markers
    markers<-Identify_Marker_TFs(x=x,group_by=group_by,assay="chromvar",assay_name="chromvar")
    markers_out<-as.data.frame(
      markers %>% group_by(chromvar.group) %>% slice_max(order_by=chromvar.auc,n=n_markers)
      )

    markers_out<-markers_out[!duplicated(markers_out$gene),]
    dim(markers_out)
    tf_motif<-average_features(x=x,features=markers_out$chromvar.feature,assay="chromvar",group_by=plot_by)
    tf_motif<-tf_motif[row.names(tf_motif) %in% markers_out$chromvar.feature,]
    #average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering
    average_matrix=tf_motif #just cluster only on TF motifs

    #set up heatmap seriation and order by GA
    o_rows =dist(average_matrix) %>%
                          hclust() %>%
                          as.dendrogram()  #%>%
                          #ladderize()
    o_col =dist(t(average_matrix),method="maximum") %>%
                      hclust() %>%
                      as.dendrogram()  %>%
                      ladderize()

    side_ha_motif<-data.frame(chromvar_motif=markers_out[get_order(o_rows,1),]$chromvar.auc)
    #colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.90,0.95)),cividis(3))
    colfun_motif=colorRamp2(c(-2,-1,0,1,2),cividis(5))

    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-markers_out$chromvar.feature
    
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

    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(markers_out)), 
                                            labels =markers_out$gene,
                                            labels_gp=gpar(fontsize=6)),
                            motifs = anno_image(motif_plots))

    o_col_split<-unlist(lapply(strsplit(colnames(tf_motif)," "),'[',1))

    motif_auc<-Heatmap(side_ha_motif,
        cluster_rows = o_rows,
        col=side_ha_col,
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
        column_split = o_col_split,
          cluster_column_slices = FALSE,
          right_annotation=gene_ha)
      
      #motif_image<-anno_image(paste0(prefix,".tf.heatmap.motif.svg"))
    
    pdf(paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf")))
    print(draw(motif_auc+motif_plot,row_title=prefix))
    dev.off()
    print(paste("Plotted... ",paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf"))))
}

#######################
#cancer only idc and ilc
#######################

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
#15 IDC vs 5 ILC, 100 cells each
pairwise_comparison(obj=dat_cancer,
                    group_by="Diagnosis",
                    group1="IDC",
                    group2="ILC",
                    outname="diagnosis",
                    motif_name="ESR1",
                    col=hist_col,
                    downsample_cells_per_sample=100,
                    outdir=paste0(output_directory,"/pairwise_by_diagnosis")
                    )

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))        
dat_cancer<-subset(dat_cancer, Diagnosis %in% c("IDC","ILC"))
dat_cancer$sample_diag<-paste(dat_cancer$Diagnosis,dat_cancer$sample)
Idents(dat_cancer)<-dat_cancer$sample_diag
plot_top_tf_markers_tfonly(x=dat_cancer,
                    group_by="Diagnosis",
                    plot_by="sample_diag",
                    prefix="pairwise_by_diagnosis",
                    n_markers=20,
                    order_by_idents=TRUE,
                    outdir=paste0(output_directory,"/pairwise_by_diagnosis"))

#######################
#PR+/- of cancer only IDC
#######################

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
dat_cancer<-subset(dat_cancer, Mol_Diagnosis %in% c("ER+/PR+/HER2-","ER+/PR-/HER2-"))
#5 IDC ER+/PR-/HER2- vs 8 IDC ER+/PR+/HER2- 95 cells each
pairwise_comparison(obj=dat_cancer,
                    group_by="Mol_Diagnosis",
                    group1="ER+/PR+/HER2-",
                    group2="ER+/PR-/HER2-",
                    outname="PR.subtype",
                    motif_name="ESR1",
                    col=clin_col,
                    downsample_cells_per_sample=95,
                    outdir=paste0(output_directory,"/pairwise_by_moleculardiag"))

rna<-read.table(paste0(output_directory,"/pairwise_by_moleculardiag/","pairwise.PR.subtype.SCT.tsv"),header=T)
rna<-rna %>% filter(padj<0.05) %>% filter(feature %in% c("RANKL","TNFRSF11A","CCND1","CDKN1A","DUSP1","EGFR","PGR","PGRMC1","IGF1R","AR","MKI67","FGFR4","LCK","FRK","MST1R"))

atac<-read.table(paste0(output_directory,"/pairwise_by_moleculardiag/","pairwise.PR.subtype.GeneActivity.tsv"),header=T)
datac<-atac %>% filter(padj<0.05) %>% filter(feature %in% c("RANKL","TNFRSF11A","CCND1","CDKN1A","DUSP1","EGFR","PGR","PGRMC1","IGF1R","AR","MKI67","FGFR4","LCK","FRK","MST1R"))

Idents(dat_cancer)<-paste(dat_cancer$sample,dat_cancer$Mol_Diagnosis)
coverage_plot(obj=dat_cancer,markers_rna=rna,markers_ga=atac,col=clin_col,outname="apriori_genes",group1="ER+/PR+/HER2-",group2="ER+/PR-/HER2-",group_by="Mol_Diagnosis",outdir=paste0(output_directory,"/pairwise_by_moleculardiag"))

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
dat_cancer<-subset(dat_cancer, Mol_Diagnosis %in% c("ER+/PR+/HER2-","ER+/PR-/HER2-"))
dat_cancer$moldiag<-paste(dat_cancer$Mol_Diagnosis,dat_cancer$sample)
Idents(dat_cancer)<-dat_cancer$moldiag
plot_top_tf_markers_tfonly(x=dat_cancer,
                    group_by="Mol_Diagnosis",
                    plot_by="moldiag",
                    prefix="pairwise_moleculardiag",
                    n_markers=20,
                    order_by_idents=FALSE,
                    outdir=paste0(output_directory,"/pairwise_by_moleculardiag"))

#######################
#scsubtype of cancer only IDC
#######################

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
#9 IDC SC_Subtype_LumA_SC vs 13 SC_Subtype_LumB_SC 50 cells each
pairwise_comparison(obj=dat_cancer,
                    group_by="scsubtype",
                    group1="SC_Subtype_LumA_SC",
                    group2="SC_Subtype_LumB_SC",
                    outname="IDC.scsubtype",
                    motif_name="ESR1",
                    col=scsubtype_col,
                    downsample_cells_per_sample=50,
                    outdir=paste0(output_directory,"/pairwise_by_scsubtype"))

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
dat_cancer<-subset(dat_cancer, scsubtype %in% c("SC_Subtype_LumA_SC","SC_Subtype_LumB_SC"))
dat_cancer$sample_scsubtype<-paste(dat_cancer$scsubtype,dat_cancer$sample)
Idents(dat_cancer)<-dat_cancer$sample_scsubtype
plot_top_tf_markers_tfonly(x=dat_cancer,
                    group_by="scsubtype",
                    plot_by="sample_scsubtype",
                    prefix="pairwise_scsubtype_by_sample",
                    n_markers=20,
                    order_by_idents=FALSE,
                    outdir=paste0(output_directory,"/pairwise_by_scsubtype"))

# ####################################################
# #           Fig 3 Heatmap By Clones                #
# ###################################################

# dat_cnv<-subset(dat,cells=names(dat$merged_assay_clones[!is.na(dat$merged_assay_clones)]))
# plot_top_tf_markers(x=dat_cnv,group_by="merged_assay_clones",prefix="clones",n_markers=3,order_by_idents=FALSE)

# ################################################
# ####Correlation of CNV count to SCT value####
# ################################################

# clone_filter<-names(which(table(dat$merged_assay_clones)>=30))
# dat_cnv<-subset(dat,cells=names(dat$merged_assay_clones[!is.na(dat$merged_assay_clones)]))

# windows<-data.frame(chr=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",1)),
#                     start=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",2)),
#                     end=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",3)))
# windows<-makeGRangesFromDataFrame(windows)

# #fetch cnv relevant to each gene and correlate cnv profile of cells with RNA expression
# gene_cnv_cor<-function(i,assay="SCT"){
#   cnv_name<-windows[queryHits(hits)[i],]
#   cnv_name<-paste(seqnames(cnv_name),start(cnv_name),end(cnv_name),sep="-")
#   cnv_val<-FetchData(dat_cnv[["cnv"]],vars=cnv_name,layer="data")

#   gene_name<-cnv_genes_windows[subjectHits(hits)[i],]$gene_name
#   gene_val<-FetchData(dat_cnv[[assay]],vars=gene_name,layer="data")

#   out<-cor(gene_val,cnv_val)
#   out<-c(row.names(out),colnames(out),unname(out))

#   return(out)
#   }


# #process RNA
# all_genes<-Features(dat_cnv,assay="SCT")
# annot<-dat@assays$ATAC@annotation
# all_genes<-all_genes[all_genes %in% annot$gene_name]
# cnv_genes_windows<-annot[annot$gene_name %in% all_genes,] #filter annotation to genes we want
# cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
# hits<-findOverlaps(query=windows,subject=cnv_genes_windows)

# rna_out<-mclapply(1:length(hits),gene_cnv_cor,mc.cores=10)
# rna_out<-as.data.frame(do.call("rbind",rna_out))

# colnames(rna_out)<-c("gene","cnv_window","rna_correlation")
# rna_out<-rna_out[complete.cases(rna_out),]
# rna_out$rna_correlation<-as.numeric(rna_out$rna_correlation)


# #process GeneActivity
# all_genes<-Features(dat_cnv,assay="GeneActivity")
# annot<-dat@assays$ATAC@annotation
# all_genes<-all_genes[all_genes %in% annot$gene_name]
# cnv_genes_windows<-annot[annot$gene_name %in% all_genes,] #filter annotation to genes we want
# cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
# hits<-findOverlaps(query=windows,subject=cnv_genes_windows)

# ga_out<-mclapply(1:length(hits),gene_cnv_cor,mc.cores=10,assay="GeneActivity")
# ga_out<-as.data.frame(do.call("rbind",ga_out))
# colnames(ga_out)<-c("gene","cnv_window","ga_correlation")
# ga_out<-ga_out[complete.cases(ga_out),]
# ga_out$ga_correlation<-as.numeric(ga_out$ga_correlation)

# combined_out<-merge(rna_out,ga_out,by=c("gene","cnv_window"))

# #add cnv variance per window
# cnv_var<-apply(dat_cnv[["cnv"]]@data, 1, var)
# combined_out$cnv_var<-cnv_var[combined_out$cnv_window]

# combined_out$label<-NA
# top_cor<-combined_out %>% arrange(desc(rna_correlation)) %>% head(n=20)
# bottom_cor<-combined_out %>% arrange(desc(rna_correlation)) %>% tail(n=10)

# combined_out[combined_out$gene %in% top_cor$gene,]$label<-top_cor$gene
# combined_out[combined_out$gene %in% bottom_cor$gene,]$label<-bottom_cor$gene

# plt<-ggplot(combined_out,
#             aes(
#               x=rna_correlation,
#               y=ga_correlation,
#               color=cnv_var,
#               label=label))+
#               geom_point()+
#               geom_text_repel(max.overlaps=Inf)+
#               theme_minimal()
# ggsave(plt,file="rna_ga_cnv_correlation.dotplot.pdf")


