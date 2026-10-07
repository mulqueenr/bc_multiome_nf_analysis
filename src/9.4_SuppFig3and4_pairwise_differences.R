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
setwd("/home/groups/MohammedLab/bc_multiome/seurat_object")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)


outdir="/home/groups/MohammedLab/bc_multiome/suppfig3"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}


setwd(outdir)

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

DefaultAssay(dat)<-"ATAC"
dat <- RegionStats(dat, genome = BSgenome.Hsapiens.UCSC.hg38)

#### Pairwise comparisons


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
  #https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
  
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

  group1_enriched<-da_combined %>% filter(padj.x<0.05) %>% filter(padj.y<0.05) %>% filter(logFC.x>0 & logFC.y>0) %>% slice_max(avg_logFC,n=10)
  group2_enriched<-da_combined %>% filter(padj.x<0.05) %>% filter(padj.y<0.05) %>% filter(logFC.x<0 & logFC.y<0) %>% slice_min(avg_logFC,n=10)
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
    tf_rna<-tf_rna[complete.cases(tf_rna),]
    tf_ga<-average_features(x=x,features=markers_out$gene,assay="GeneActivity",group_by=plot_by)
    tf_ga<-tf_ga[row.names(tf_ga) %in% markers_out$gene,]
    tf_ga<-tf_ga[complete.cases(tf_ga),]

    tf_motif<-average_features(x=x,features=markers_out$chromvar.feature,assay="chromvar",group_by=plot_by)
    tf_motif<-tf_motif[row.names(tf_motif) %in% markers_out$chromvar.feature,]
    tf_motif<-tf_motif[complete.cases(tf_motif),]

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

    o_cols =t(average_matrix) %>% dist()  %>% 
                          hclust() %>%
                          as.dendrogram() %>%
                          ladderize(decreasing=FALSE) %>% labels()

    o_rows = average_matrix %>% dist()  %>% 
                          hclust() %>%
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

```


```R
#######################
#PR+/- of cancer only IDC
#######################

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
dat_cancer<-subset(dat_cancer, Mol_Diagnosis %in% c("ER+/PR+/HER2-","ER+/PR-/HER2-"))
dat_cancer<-subset(dat_cancer, cells=row.names(dat_cancer@meta.data)[which(!endsWith(dat_cancer$merged_assay_clones,suffix="_normal"))]) 

#5 IDC ER+/PR-/HER2- vs 8 IDC ER+/PR+/HER2- 95 cells each
pairwise_comparison(obj=dat_cancer,
                    group_by="Mol_Diagnosis",
                    group1="ER+/PR+/HER2-",
                    group2="ER+/PR-/HER2-",
                    outname="PR.subtype",
                    motif_name="ESR1",
                    col=clin_col,
                    downsample_cells_per_sample=95,
                    outdir=outdir)

rna<-read.table(paste0(outdir,"/pairwise.PR.subtype.SCT.tsv"),header=T)
rna<-rna %>% filter(padj<0.05) %>% filter(feature %in% c("RANKL","TNFRSF11A","CCND1","CDKN1A","DUSP1","EGFR","PGR","PGRMC1","IGF1R","AR","MKI67","FGFR4","LCK","FRK","MST1R"))

atac<-read.table(paste0(output_directory,"/pairwise_by_moleculardiag/","pairwise.PR.subtype.GeneActivity.tsv"),header=T)
atac<-atac %>% filter(padj<0.05) %>% filter(feature %in% c("GRHL2","NRIP1","TRPS1","CO4A","TLE3","GATA3","FKBP4","FKBP5","HS90A","HS90B","RANKL","TNFRSF11A","CCND1","CDKN1A","DUSP1","EGFR","PGR","PGRMC1","IGF1R","AR","MKI67","FGFR4","LCK","FRK","MST1R"))

Idents(dat_cancer)<-paste(dat_cancer$merged_assay_clones)

coverage_plot(obj=dat_cancer,markers_rna=rna,markers_ga=atac,col=clin_col,outname="apriori_genes",group1="ER+/PR+/HER2-",group2="ER+/PR-/HER2-",group_by="Mol_Diagnosis",outdir=paste0(output_directory,"/pairwise_by_moleculardiag"))

col=clin_col
region_i="chr11-101120114-101134860"
  annot_plot<-AnnotationPlot(object=dat_cancer,region=region_i)
  cov_plot <- CoveragePlot(
    object = dat_cancer,
    region = region_i,
    group.by = "Mol_Diagnosis",
    split.by = "merged_assay_clones",
    annotation = FALSE,
    peaks = TRUE,links=FALSE)+
    scale_fill_manual(values=col)

  expr_plot <- ExpressionPlot(
    object = dat_cancer,
    group.by = "merged_assay_clones",
    features = "PGR",
    assay = "SCT") + scale_fill_manual(values=col)

 # link_plot_1 <- LinkPlot(
 #   object = "ER+/PR+/HER2-",
 #   region = region_i)+
 #   scale_color_gradient2(limits=c(0,0.3),low="white",high=col[group1])

  #link_plot_2 <-LinkPlot(
  #  object = "ER+/PR-/HER2-",
  #  region = region_i)+
  #  scale_color_gradient2(limits=c(0,0.3),low="white",high=col[group2])

  plt<-CombineTracks(
    plotlist = list(cov_plot, annot_plot ),#link_plot_1,link_plot_2
    expression.plot = expr_plot,
    heights = c(10, 2), #3, 3
    widths = c(10, 3))
ggsave(plt,file=paste0(outdir,"/","PGR_promoter.clones.coverage.pdf"),height=20)

colfun <- colorRamp2(
  #breaks = c(-3,-1,-0.5,0,0.5,1,3),
  breaks = c(-3,-2,-1,0,1,2,3),
  colors = c("#053061","#487590","#d1e5f0","#f7f7f7","#fddbc7","#a76146","#67001f"))


Idents(dat_cancer)<-dat_cancer$Diag_MolDiag
Idents(dat_cancer)<-factor(dat_cancer$merged_assay_clones,levels=levels(reorder(dat_cancer$merged_assay_clones,dat_cancer$Diag_MolDiag)))
plot_top_tf_markers(x=dat_cancer,
                    group_by="Diag_MolDiag",
                    plot_by="merged_assay_clones",
                    prefix="pairwise_moleculardiag_pr",
                    colfun=colfun,
                    n_markers=20,
                    order_by_idents=TRUE,
                    outdir=outdir)
```


IDC vs ILC PR+
```r
#######################
#IDC vs ILC
#######################

dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC","ILC"))
dat_cancer<-subset(dat_cancer, Mol_Diagnosis %in% c("ER+/PR+/HER2-"))
dat_cancer<-subset(dat_cancer, cells=row.names(dat_cancer@meta.data)[which(!endsWith(dat_cancer$merged_assay_clones,suffix="_normal"))]) 

#5 IDC ER+/PR-/HER2- vs 8 IDC ER+/PR+/HER2- 95 cells each
pairwise_comparison(obj=dat_cancer,
                    group_by="Diagnosis",
                    group1="IDC",
                    group2="ILC",
                    outname="IDC_ILC.subtype",
                    motif_name="ESR1",
                    col=hist_col,
                    downsample_cells_per_sample=95,
                    outdir=outdir)

rna<-read.table(paste0(outdir,"/pairwise.IDC_ILC.subtype.SCT.tsv"),header=T)
rna<-rna %>% filter(padj<0.05) %>% filter(feature %in% c("RANKL","TNFRSF11A","CCND1","CDKN1A","DUSP1","EGFR","PGR","PGRMC1","IGF1R","AR","MKI67","FGFR4","LCK","FRK","MST1R"))

atac<-read.table(paste0(outdir,"/pairwise.IDC_ILC.subtype.GeneActivity.tsv"),header=T)
atac<-atac %>% filter(padj<0.05) %>% filter(feature %in% c("GRHL2","NRIP1","TRPS1","CO4A","TLE3","GATA3","FKBP4","FKBP5","HS90A","HS90B","RANKL","TNFRSF11A","CCND1","CDKN1A","DUSP1","EGFR","PGR","PGRMC1","IGF1R","AR","MKI67","FGFR4","LCK","FRK","MST1R"))

Idents(dat_cancer)<-paste(dat_cancer$merged_assay_clones)

#coverage_plot(obj=dat_cancer,markers_rna=rna,markers_ga=atac,col=hist_col,outname="apriori_genes",group1="IDC",group2="ILC",group_by="Mol_Diagnosis",outdir=outdir)

#plotting CDH1 activity

col=hist_col
chr16:68,737,292 - 68,835,540
region_i="chr16-68735292-68837540"
  annot_plot<-AnnotationPlot(object=dat_cancer,region=region_i)
  cov_plot <- CoveragePlot(
    object = dat_cancer,
    region = region_i,
    group.by = "Diagnosis",
    split.by = "merged_assay_clones",
    annotation = FALSE,
    peaks = TRUE,links=FALSE)+
    scale_fill_manual(values=col)

  expr_plot <- ExpressionPlot(
    object = dat_cancer,
    group.by = "merged_assay_clones",
    features = "CDH1",
    assay = "SCT") + scale_fill_manual(values=col)

 # link_plot_1 <- LinkPlot(
 #   object = "ER+/PR+/HER2-",
 #   region = region_i)+
 #   scale_color_gradient2(limits=c(0,0.3),low="white",high=col[group1])

  #link_plot_2 <-LinkPlot(
  #  object = "ER+/PR-/HER2-",
  #  region = region_i)+
  #  scale_color_gradient2(limits=c(0,0.3),low="white",high=col[group2])

  plt<-CombineTracks(
    plotlist = list(cov_plot, annot_plot ),#link_plot_1,link_plot_2
    expression.plot = expr_plot,
    heights = c(10, 2), #3, 3
    widths = c(10, 3))
ggsave(plt,file=paste0(outdir,"/","CDH1.clones.coverage.pdf"),height=20)


colfun <- colorRamp2(
  #breaks = c(-3,-1,-0.5,0,0.5,1,3),
  breaks = c(-3,-2,-1,0,1,2,3),
  colors = c("#053061","#487590","#d1e5f0","#f7f7f7","#fddbc7","#a76146","#67001f"))


Idents(dat_cancer)<-factor(dat_cancer$merged_assay_clones,levels=levels(reorder(dat_cancer$merged_assay_clones,dat_cancer$Diagnosis)))
plot_top_tf_markers(x=dat_cancer,
                    group_by="Diagnosis",
                    plot_by="merged_assay_clones",
                    prefix="pairwise_idc_ilc",
                    colfun=colfun,
                    n_markers=10,
                    order_by_idents=TRUE,
                    outdir=outdir)

```