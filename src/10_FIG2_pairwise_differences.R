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
#tornado plot of top DA peaks
#volcano plot of top DE genes and chromvar TFs
#gsea of top DE genes

tornado_plot<-function(obj=obj,da_peak_set=markers,i="IDC",peak_count=200,col=col,col_lim=0.03){
    print(i)
    top_peaks <- da_peak_set %>% 
                  filter(group==i) %>%  
                  filter(logFC>0) %>% 
                  arrange(padj,desc(auc)) %>% 
                  slice_head(n=peak_count)
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

volcano_plot<-function(obj=obj,de_features_set=markers,prefix,feature_count=25,outname,assay,group1,group2,col){
  
  feature_markers_group1<-de_features_set %>% filter(padj<0.05) %>% arrange(desc(logFC)) %>% head(n=feature_count)
  feature_markers_group2<-de_features_set %>% filter(padj<0.05) %>% arrange(logFC) %>% head(n=feature_count)
  group_1_count=de_features_set %>% filter(padj<0.05) %>% filter(logFC>0) %>% nrow()
  group_2_count=de_features_set %>% filter(padj<0.05) %>% filter(logFC<0) %>% nrow()

  de_features_set$fill_col<-"#808080"
  de_features_set[de_features_set$padj<0.05 & de_features_set$logFC>0,]$fill_col<-col[group1]
  de_features_set[de_features_set$padj<0.05 & de_features_set$logFC<0,]$fill_col<-col[group2]

  de_features_set$label<-NA
  de_features_set[de_features_set$feature %in% feature_markers_group1$feature,]$label<-feature_markers_group1$feature
  de_features_set[de_features_set$feature %in% feature_markers_group2$feature,]$label<-feature_markers_group2$feature
  de_features_set$padj_pseudocount<-de_features_set$padj+1e-30

  plt<-ggplot(de_features_set,
              aes(x=logFC,
                  y=-log10(padj_pseudocount),
                  color=fill_col,
                  label=label))+
        geom_point()+geom_hline(yintercept=-log10(0.05))+
        scale_color_identity()+ coord_cartesian(clip = "off") + theme_minimal() +
        ggrepel::geom_text_repel(size=2,max.overlaps=Inf,min.segment.length = 0,vjust = "inward")+
        ggtitle(paste(group1,":",as.character(group_1_count),
        "\n",group2,":",as.character(group_2_count)))

  ggsave(plt,file=paste0("pairwise.volcano.",assay,".",outname,".pdf"))
}

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
    dplyr::filter(padj<0.05) %>% 
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(feature, logFC)

  ranks<-setNames(nm=group1_features$feature,group1_features$logFC)

  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 10,
                    nproc = 1)

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))
  #plot gsea ranking of genes by top pathways
  #pdf(paste0("pairwise.",outname,".",assay,".",out_setname,".gsea.pdf"),width=20,height=10)
  #print(plt)
  #dev.off()

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% arrange(desc(NES)) %>% head(n= 10),
    fgseaRes  %>% arrange(desc(NES)) %>% tail(n= 10))
  
  nes_plt_dat$col<-"#808080"
  nes_plt_dat[nes_plt_dat$NES>0 & nes_plt_dat$pval<0.05,]$col<-col[group1]
  nes_plt_dat[nes_plt_dat$NES<0 & nes_plt_dat$pval<0.05,]$col<-col[group2]

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= col)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(out_setname)
  return(patchwork::wrap_plots(list(plt1,plt2),ncol=2))
}

plot_gsea<-function(obj,annot,dmrs,
                    outname=outname,
                    assay=assay,col=col,group1,group2){

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
  plt<-patchwork::wrap_plots(list(tft_plt,position_plt,hallmark_plt),nrow=3)
  ggsave(plt,file=paste0("pairwise.",outname,".",assay,".gsea.NES.pdf"),width=20,height=10)
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

coverage_plot<-function(obj,marker_rna,marker_ga,col,outname,group1,group2,group_by){
  DefaultAssay(obj)<-"ATAC"

  #rank by average AUC
  da_combined<-merge(markers_rna,markers_ga,by="feature")
  da_combined$avg_logFC<-rowMeans(da_combined[,c('logFC.x', 'logFC.y')], na.rm=TRUE) #dont actually need this
  da_combined$avg_AUC<-rowMeans(da_combined[,c('auc.x', 'auc.y')], na.rm=TRUE) #dont actually need this
  group1_enriched<-da_combined %>% filter(padj.x<0.05) %>% filter(padj.y<0.05) %>%
                    arrange(desc(avg_AUC)) %>% slice_head(n=10)

  group2_enriched<-da_combined %>% filter(padj.x<0.05) %>% filter(padj.y<0.05) %>%
                     arrange(desc(avg_AUC)) %>% slice_tail(n=10)
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
  ggsave(plt,file=paste0("pairwise.coverage.",outname,".pdf"),width=50,height=20,limitsize=FALSE)
  
}

pairwise_comparison<-function(obj=dat_cancer,
                              group_by="Diagnosis",
                              group1="IDC",
                              group2="ILC",
                              outname="diagnosis",
                              motif_name="ESR1",
                              col){
  #subset to relevent groups
  obj <-subset(obj, cells=row.names(obj@meta.data)[obj@meta.data[,group_by] %in% c(group1,group2)])
  DefaultAssay(obj)<-"ATAC"
  annot<-Annotation(obj)

  #downsample to ~ equal cell counts per sample
  obj_full<-obj
  Idents(obj)<-obj$sample
  downsample_cell_table=table(Idents(obj))
  obj <-subset(obj,downsample=50)
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
  
  markers<- markers %>% filter(group==group1) %>% filter(padj<0.05)
  write.table(markers,col.names=T,row.names=F,file=paste0("pairwise.all_da_peaks.",outname,".tsv"),sep="\t")
  plt_list<-lapply(unique(markers$group),function(j) {
      tornado_plot(obj=obj,da_peak_set=markers,i=j,col=col)})
  plt<-wrap_plots(plt_list,nrow=1,guides='collect')
  ggsave(plt,file=paste0("pairwise.tornado.all_da_peaks.",outname,".pdf"),width=20,height=20)
  
  #ESR1 only
  #da peaks that overlap with ESR1 motif only
  motif<-names(obj@assays$ATAC@motifs@motif.names[which(obj@assays$ATAC@motifs@motif.names==motif_name)])
  da_peaks_motif_filt<-markers[markers$feature %in% names(which(obj@assays$ATAC@motifs@data[,motif])),]
  plt_list<-lapply(unique(da_peaks_motif_filt$group),function(j) {
      tornado_plot(obj=obj,da_peak_set=da_peaks_motif_filt,i=j,col=col)
      })
  plt<-wrap_plots(plt_list,nrow=1,guides='collect')
  ggsave(plt,file=paste0("pairwise.tornado.ESR1only_peaks.",outname,".pdf"),width=20,height=20)
  
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
 
  write.table(markers_rna,col.names=T,row.names=F,file=paste0("pairwise.",assay,".",outname,".tsv"),sep="\t")

  #all genes
  volcano_plot(obj=obj,
              de_features_set=markers_rna,
              feature_count=25,
              outname=paste0(outname,".allgenes"),
              assay=assay,group1=group1,group2=group2,col=col)

  plot_gsea(obj=obj,
            dmrs=markers_rna,
            outname=paste0(outname,".allgenes"),
            assay=assay,group1=group1,group2=group2,col=col,annot=annot)
  
  #protein coding only
  print("Running SCT pairwise comparisons (protein coding only)...")
  markers_rna<-markers_rna[markers_rna$feature %in% annot[annot$gene_biotype=="protein_coding",]$gene_name,]
  volcano_plot(obj=obj,
              de_features_set=markers_rna,
              feature_count=25,
              outname=paste0(outname,".proteincoding"),
              assay=assay,group1=group1,group2=group2,col=col)

  plot_gsea(obj=obj,
            dmrs=markers_rna,
            outname=paste0(outname,".proteincoding"),
            assay=assay,group1=group1,group2=group2,col=col,annot=annot)

  # ##########GENE ACTIVITY plots##############
  print("Running Gene Activity pairwise comparisons...")
  assay="GeneActivity"

  markers_ga <- presto:::wilcoxauc.Seurat(X = obj, 
                                      group_by = group_by,
                                      groups_use=c(group1,group2),
                                      seurat_assay = assay,
                                      assay="data")
   markers_ga<- markers_ga %>% filter(group==group1)
  write.table(markers_ga,col.names=T,row.names=F,file=paste0("pairwise.",assay,".",outname,".tsv"),sep="\t")

   volcano_plot(obj=obj,
               de_features_set=markers_ga,
               feature_count=25,
               outname=outname,
               assay=assay,group1=group1,group2=group2,col=col)

   plot_gsea(obj=obj,
             dmrs=markers_ga,
             outname=outname,
             assay=assay,group1=group1,group2=group2,col=col,annot=annot)

  print("Running Gene Activity pairwise comparisons (protein coding only)...")
  markers_ga<-markers_ga[markers_ga$feature %in% annot[annot$gene_biotype=="protein_coding",]$gene_name,]

  volcano_plot(obj=obj,
              de_features_set=markers_ga,
              feature_count=25,
              outname=paste0(outname,".proteincoding"),
              assay=assay,group1=group1,group2=group2,col=col)

  plot_gsea(obj=obj,
            dmrs=markers_ga,
            outname=paste0(outname,".proteincoding"),
            assay=assay,group1=group1,group2=group2,col=col,annot=annot)

  ########Coverage plots of GA and RNA################
  coverage_plot(obj=obj,
                marker_rna=markers_rna,
                marker_ga=markers_ga,
                col=col,
                outname=outname,
                group1=group1,
                group2=group2,
                group_by=group_by)

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
  write.table(markers_tf,col.names=T,row.names=F,file=paste0("pairwise.",assay,".",outname,".tsv"),sep="\t")

  volcano_plot(obj=obj,
              de_features_set=markers_tf,
              feature_count=25,
              outname=outname,
              assay=assay,group1=group1,group2=group2,col=col)

}

#cancer only idc and ilc
dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
pairwise_comparison(obj=dat_cancer,
                    group_by="Diagnosis",
                    group1="IDC",
                    group2="ILC",
                    outname="diagnosis",
                    motif_name="ESR1",
                    col=hist_col)

#PR+/- of cancer only IDC
dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
pairwise_comparison(obj=dat_cancer,
                    group_by="Mol_Diagnosis",
                    group1="ER+/PR+/HER2-",
                    group2="ER+/PR-/HER2-",
                    outname="PR.subtype",
                    motif_name="ESR1",
                    col=clin_col)

#scsubtype of cancer only IDC
dat_cancer<-subset(dat,assigned_celltype %in% c("cancer"))
dat_cancer<-subset(dat_cancer,Diagnosis %in% c("IDC"))
pairwise_comparison(obj=dat_cancer,
                    group_by="scsubtype",
                    group1="SC_Subtype_LumA_SC",
                    group2="SC_Subtype_LumB_SC",
                    outname="IDC.scsubtype",
                    motif_name="ESR1",
                    col=scsubtype_col)


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


# marker<-list()
# #iglesia et al
# marker[["her2_tumors"]]<-c("ERBB2","GRB7")
# marker[["lum_tumor"]]<-c("ESR1","GRHL2")
# marker[["basal-like"]]<-c("SFRP1","KRT17")
# #nakshatri et al
# marker[["tcells"]]<-c("IL7R","IFNG","GZMK","FCGR3A")
# marker[["macrophage"]]<-c("LYVE1","ACKR1")
# marker[["fibroblast"]]<-c("CXCL12","GLI2")
# marker[["lum_hs"]]<-c("ESR1","FOXA1","GATA3")
# marker[["lum_hs_asp"]]<-c("EHF","ELF5","KIT")
# marker[["basal"]]<-c("TP63","KRT14")