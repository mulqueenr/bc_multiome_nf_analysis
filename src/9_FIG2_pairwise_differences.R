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
library(fgsea)
library(msigdbr)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

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


#### Pairwise comparisons
#tornado plot of top DA peaks
#volcano plot of top DE genes, geneactivity, and chromvar TFs
#gsea of top DE genes, geneactivity

#add cols in regionheatmap call
tornado_plot<-function(obj=obj,da_peak_set=markers,i="IDC",peak_count=100,col=col){
    print(i)
    top_peaks <- da_peak_set %>% 
                  filter(group==i) %>%  
                  filter(logFC>0) %>% 
                  arrange(padj,desc(logFC)) %>% 
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
      cols=col[i],
      nrow=length(unique(da_peak_set$group)))+ 
    ggtitle(i)

    print("Returning plot...")
    return(plt)
}

volcano_plot<-function(obj=obj,de_features_set=markers,prefix,feature_count=25,outname,assay,group1,group2,col){
  
  de_features_set<-de_features_set %>% filter(group==group1)
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

  plt<-ggplot(de_features_set,
              aes(x=logFC,
                  y=-log10(padj),
                  color=fill_col,
                  label=label))+
        geom_point()+geom_hline(yintercept=-log10(0.05))+
        scale_color_identity()+
        geom_text_repel(max.overlaps=Inf)+
        theme_minimal()+ggtitle(paste(group1,":",as.character(group_1_count),
        "\n",group2,":",as.character(group_2_count)))

  ggsave(plt,file=paste0("pairwise.volcano.",assay,".",outname,".pdf"))
}

####ADD COVERAGE PLOT FOR HIGH GA DIFFERENCES ACROSS MARKERS
gsea_enrichment<-function(species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          de_features_set,
                          col){
  pathwaysDF <- msigdbr(species=species, 
                        category=category, 
                        subcategory = subcategory)
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
                    maxSize  = 500,
                    nproc = 1)

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  plt<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)
  #plot gsea ranking of genes by top pathways
  pdf(paste0("pairwise.",outname,".",assay,".",out_setname,"gsea.pdf"),width=20,height=10)
  print(plt)
  dev.off()

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% arrange(desc(NES)) %>% head(n= 10),
    fgseaRes  %>% arrange(desc(NES)) %>% tail(n= 10))
  
  nes_plt_dat$col<-"#808080"
  nes_plt_dat[nes_plt_dat$NES>0 & nes_plt_dat$padj<0.05,]$col<-col[group1]
  nes_plt_dat[nes_plt_dat$NES<0 & nes_plt_dat$padj<0.05,]$col<-col[group2]

  plt<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= col)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()
  ggsave(plt,file=paste0("pairwise.",outname,".",assay,".",out_setname,".gsea.NES.pdf"),width=20,height=10)
}

plot_gsea<-function(obj,dmrs,
                    outname=outname,
                    assay=assay,col=col){

  #run gsea enrichment on different sets
  gsea_enrichment(species="human",
              category="C3",
              subcategory="TFT:GTRD",
              out_setname="TFT",
              de_features_set=dmrs,
              col=col)

  gsea_enrichment(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",
              de_features_set=dmrs,
              col=col)

  gsea_enrichment(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",
              de_features_set=dmrs,
              col=col)
}

pairwise_comparison<-function(obj=dat_cancer,
                              group_by="Diagnosis",
                              group1="IDC",
                              group2="ILC",
                              outname="diagnosis",
                              motif_name="ESR1",
                              col){
  #subset to relavent groups
  obj <-subset(obj, cells=row.names(obj@meta.data)[obj@meta.data[,group_by] %in% c(group1,group2)])
  Idents(obj)<-obj@meta.data[,group_by]

  ##########peaks and tornado plots##############
  assay="ATAC"
  #all peaks
  markers <- presto:::wilcoxauc.Seurat(
    X = obj, 
    group_by = group_by, 
    groups_use=c(group1,group2),
    y=c(group1,group2), 
    seurat_assay = assay)
  
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
  assay="SCT"
  markers <- presto:::wilcoxauc.Seurat(X = obj, 
                                      group_by = group_by,
                                      groups_use=c(group1,group2),
                                      seurat_assay = assay,
                                      assay="scale.data")
  
  #all genes
  volcano_plot(obj=obj,
              de_features_set=markers,
              feature_count=25,
              outname=paste0(outname,".allgenes"),
              assay=assay,group1=group1,group2=group2,col=col)

  plot_gsea(obj=obj,
            dmrs=markers,
            outname=paste0(outname,".allgenes"),
            assay=assay,col=col)
  

  #protein coding only
  annot<-Annotation(obj)
  markers<-markers[markers$feature %in% annot[annot$gene_biotype=="protein_coding",]$gene_name,]
  volcano_plot(obj=obj,
              de_features_set=markers,
              feature_count=25,
              outname=paste0(outname,".proteincoding"),
              assay=assay,group1=group1,group2=group2,col=col)

  plot_gsea(obj=obj,
            dmrs=markers,
            outname=paste0(outname,".proteincoding"),
            assay=assay,col=col)

  # ##########GENE ACTIVITY plots##############
  # assay="GeneActivity"

  # obj<-ScaleData(obj,
  #   assay="GeneActivity",
  #   vars.to.regress="nCount_GeneActivity",
  #   features=Features(obj,assay="GeneActivity"))

  # markers <- presto:::wilcoxauc.Seurat(X = obj, 
  #                                     group_by = group_by,
  #                                     seurat_assay = assay,
  #                                     assay="scale.data")

  # volcano_plot(obj=obj,
  #             de_features_set=markers,
  #             feature_count=25,
  #             outname=outname,
  #             assay=assay,
  #             group1=group1,group2=group2,col=col)

  # plot_gsea(obj=obj,
  #           dmrs=markers,
  #           outname=outname,
  #           assay=assay,col=col)

  ##########CHROMVAR plots##############
  assay="chromvar"
  markers <- presto:::wilcoxauc.Seurat(X = obj, 
                                      group_by = group_by,
                                      seurat_assay = assay,
                                      assay="scale.data")

  markers$feature<-ConvertMotifID(object=obj,
                                  assay="ATAC",
                                  id=markers$feature)
  row.names(markers)<-markers$feature   
  volcano_plot(obj=obj,
              de_features_set=markers,
              feature_count=25,
              outname=outname,
              assay=assay,group1=group1,group2=group2)

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
#          Supplementary Plots                     #
###################################################














####################################################
#           SUPP Fig 2 Coverage Plots of Promoters     #
###################################################
# first compute the GC content for each peak
DefaultAssay(dat)<-"ATAC"
dat <- RegionStats(dat, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
dat<- LinkPeaks(
  object = dat,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = unlist(marker))

#set regions to print (promoter region + 10kb into gene)
markers<-c(
"FOXA1", #lum_hs
"EPCAM", #epi
"ELF5",#lum_asp
"KRT17", #basal
"ACACB", #adipo
"LAMA2",#fibro
"VWF", #endo
"CD163", #myeloid
"IL7R"#tcell
)


for(i in markers){
    x<-i
    col<-celltype_col
    gene_annot=dat@assays$ATAC@annotation[dat@assays$ATAC@annotation$gene_name==i,]
    gene_region<-GetTSSPositions(gene_annot, biotypes = "protein_coding")
    promoter_chr<-seqnames(gene_region)
    promoter_start=ifelse(strand(gene_region)=="-",start(gene_region)+2000,start(gene_region)-2000)
    promoter_end=ifelse(strand(gene_region)=="-",start(gene_region)-10000,start(gene_region)+10000)

    if(promoter_start<promoter_end){
      promoter_region=StringToGRanges(paste(promoter_chr,promoter_start,promoter_end,sep="-"))
    } else {
      promoter_region=StringToGRanges(paste(promoter_chr,promoter_end,promoter_start,sep="-"))
    }
    if(x %in% row.names(dat@assays$SCT@data)){
      print(paste("Plotting...",x,"promoter."))

      cov_plot <- CoveragePlot(
        object = dat,
        region = promoter_region,features = x,
        annotation = TRUE,
        peaks = TRUE,links=FALSE)+ scale_color_manual(values=col)
      ggsave(cov_plot,file=paste0("set_region_coverage.",x,".pdf"),width=3,height=10)
    }
  }


Idents(dat_epi)<-dat_epi$cellstate
markers <- presto:::wilcoxauc.Seurat(X = dat_epi, group_by = "cellstate", 
  groups_use=unname(unlist(unique(dat_epi$cellstate))),
  y=unname(unlist(unique(dat_epi$cellstate))), 
  assay = 'data', seurat_assay = "ATAC")

da_peaks <- markers %>% filter(padj<0.01) %>% filter(logFC>0.25)
#[1] 39525    10

#add chromvar motif scan of motifs per DA peaks
da_peaks<-merge(da_peaks,as.data.frame(dat@assays$ATAC@motifs@data),by.x="feature",by.y= 'row.names')

#run on RNA
rna<-FindAllMarkers(dat_epi,assay="RNA")
rna$group<-rna$cluster
rna$gene_name<-rna$gene




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
#           Fig 2  TF motifs IDC v ILC           #
###################################################

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
"NAT"="#99CCFF",
"DCIS"="#CCCCCC",
"IDC"="#FF9966",
"ILC"="#006633")

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

    obj_mat<-RegionMatrix(obj,key="DA_mat",
    regions=StringToGRanges(top_peaks$gene),
    upstream=5000,downstream=5000,
    assay="ATAC")

    plt<-RegionHeatmap(obj_mat,key="DA_mat",
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


