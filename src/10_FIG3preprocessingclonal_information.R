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
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character"),
  make_option(c("-r", "--ref_object"), type="character", default="/home/groups/CEDAR/mulqueen/bc_multiome/ref/nakshatri/nakshatri_multiome.geneactivity.rds", 
              help="Nakshatri reference object for epithelial comparisons", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])

############# Add clone and filtering info to object
read_clone_annot<-function(i){
clone_annot<-list.files(i,pattern="_ATAC_anno.csv",full.names=T)
  if(length(clone_annot)!= 0 ){
    sample_name<-gsub(basename(clone_annot),pattern="_ATAC_anno.csv",replace="")
    clones<-read.csv(clone_annot,row.names=1)
    row.names(clones)<-gsub(row.names(clones),pattern="[.]",replace="-")
    row.names(clones)<-paste(sample_name,row.names(clones),sep="_")
    clones$merge_cluster<-paste(sample_name,clones$merge_cluster,sep="_")
    clones$atac_cluster<-paste(sample_name,clones$atac_cluster,sep="_")
    clones$rna_cluster<-paste(sample_name,clones$rna_cluster,sep="_")
    return(clones)
  }
}



############# Add clone and filtering info to object
read_cnv_plots<-function(i){
  sample<-gsub(basename(i),pattern="_consensus_cnvs.bed",replacement="")
  cnv<-read.table(i)
  cnv<-cnv[5:ncol(cnv)]
  colnames(cnv)<-gsub(colnames(cnv),pattern="[.]",replacement="-")
  colnames(cnv)<-paste(sample,colnames(cnv),sep="_")
  return(cnv)
}

clone_annot_list<-list.files("/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/",pattern="*clone_val",full.name=T)
clone_annot<-do.call("rbind",lapply(clone_annot_list,read_clone_annot))
colnames(clone_annot)<-c("merged_assay_clones","atac_clones","rna_clones","clones_celltype","clones_log2reads")
dat<-AddMetaData(dat,clone_annot)

cnv_list<-list.files("/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/CNV_validate/",recursive=T,pattern="*_consensus_cnvs.bed",full.name=T)
cnv_list<-cnv_list[!grepl(cnv_list,pattern="unfiltered")]
cnv_windows<-read.table(cnv_list[1])
cnv_windows<-cnv_windows[1:4]
colnames(cnv_windows)<-c("chr","start","end","win")

cnv_out<-do.call("cbind",lapply(cnv_list,read_cnv_plots))
row.names(cnv_out)<-paste(cnv_windows$chr,cnv_windows$start,cnv_windows$end,sep="-")
cnv_out<-cnv_out[colnames(cnv_out) %in% row.names(dat@meta.data)]
cnv_assay <- CreateAssayObject(counts = cnv_out)
dat[["cnv"]]<-cnv_assay

#normal clones based on lack of CNV calls
normal_clones<-c(
"DCIS_03_1",
"IDC_01_4",
"IDC_02_1",
"IDC_05_1",
"IDC_06_3",
"IDC_08_2",
"IDC_09_2",
"IDC_10_3",
"IDC_11_2",
"IDC_12_4",
"ILC_02_1",
"ILC_04_2")

contamination_clones<-c(
"ILC_04_3",
"ILC_04_4"
)

dat@meta.data[dat@meta.data$merged_assay_clones %in% normal_clones,]$merged_assay_clones<-"normal"
dat@meta.data[dat@meta.data$merged_assay_clones %in% contamination_clones,]$merged_assay_clones<-"contamination"
dat@meta.data[dat@meta.data$merged_assay_clones %in% c("normal"),]$merged_assay_clones <-paste(dat@meta.data[dat@meta.data$merged_assay_clones %in% c("normal"),]$sample,"normal",sep="_")

saveRDS(dat,file="8_merged.cnv_clones.SeuratObject.rds")


#####Plot of heatmap of all clones ####
cnv_col<-c("0"="#002C3E", "0.5"="#78BCC4", "1"="#F7F8F3", "1.5"="#F7444E", "2"="#aa1407", "3"="#440803")
#from Curtis et al.

windows<-data.frame(chr=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",1)),
                    start=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",2)),
                    end=unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",3)))
          
windows<-makeGRangesFromDataFrame(windows)


#from https://www.nature.com/articles/s41416-024-02804-6#Sec20
#change RAB7L1 to RAB29
#lost RAB7L1

cnv_genes<-c('DLEU2L', 'TRIM46', 'FASLG', 'KDM5B', 'RAB7L1', 'PFN2', 'PIK3CA', 'EREG', 'AIM1', 'EGFR', 'ZNF703', 'MYC', 'SEPHS1', 'ZMIZ1', 'EHF', 'POLD4', 'CCND1', 'P2RY2', 'NDUFC2-KCTD14', 'FOXM1', 'MDM2', 'STOML3', 'NEMF', 'IGF1R', 'TP53I13', 'ERBB2', 'SGCA', 'RPS6KB1', 'BIRC5', 'NOTCH3', 'CCNE1', 'RCN3', 'SEMG1', 'ZNF217', 'TPD52L2', 'PCNT', 'CDKN2AIP', 'LZTS1', 'PPP2R2A', 'CDKN2A', 'PTEN', 'RB1', 'CAPN3', 'CDH1', 'MAP2K4', 'GJC2', 'TERT', 'RAD21', 'ST3GAL1', 'SOCS1')
cnv_genes_class<-c('amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'amp', 'amp', 'amp', 'amp', 'amp')
cnv_genes<-setNames(cnv_genes_class,cnv_genes)
cnv_genes<-cnv_genes[names(cnv_genes) %in% dat@assays$ATAC@annotation$gene_name]
cnv_genes_windows<-dat@assays$ATAC@annotation[dat@assays$ATAC@annotation$gene_name %in% names(cnv_genes),] #filter annotation to genes we want
cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
windows<-findOverlaps(windows,cnv_genes_windows)
annot<-data.frame(
  window_loc=queryHits(windows),
  gene=cnv_genes_windows$gene_name,
  cnv_class=unname(cnv_genes[cnv_genes_windows$gene_name]))

annot$col<-ifelse(annot$cnv_class=="amp","red","blue")

hc = columnAnnotation(common_cnv = anno_mark(at = annot$window_loc, 
                        labels = annot$gene,
                        which="column",side="bottom",
                        labels_gp=gpar(col=annot$col)))



pdf("all_samples.cnv.heatmap.pdf",height=50,width=50)
Heatmap(t(as.data.frame(dat@assays$cnv@counts)),
  col=cnv_col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = FALSE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  bottom_annotation=hc,
  row_split=dat@meta.data[row.names(dat@meta.data) %in% colnames(dat@assays$cnv@counts),]$merged_assay_clones,
  column_split=factor(unlist(lapply(strsplit(row.names(dat@assays$cnv@counts),"-"),"[",1)),levels=paste0("chr",1:22)),
  border = TRUE)
dev.off()



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
#use only clones with at least 30 cell 
clone_filter<-names(which(table(dat$merged_assay_clones)>=30))

cellcycle_freq<-dat@meta.data %>% 
  filter(assigned_celltype=="cancer") %>%
  filter(sample %in% cancercell_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,Phase,merged_assay_clones) %>% 
  count(Phase,.drop=TRUE)


cellcycle_freq<-dat@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) %>% 
  count(Phase,.drop=TRUE)

plt<-ggplot(cellcycle_freq, aes(fill=factor(Phase,levels=names(Phase_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=Phase_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="phase_perclone_percbarplot.pdf",width=10,height=10)


#box plot of wu scores
emt_and_d_scores<-dat@meta.data %>% 
  filter(assigned_celltype=="cancer") %>%
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>%
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,merged_assay_clones) 

plt<-ggplot(emt_and_d_scores, aes(y=Wu_EMTScores, x=merged_assay_clones)) + 
  geom_boxplot(width=1) +
  geom_jitter(width = 1) +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="emtscore_perclone_boxplot.pdf",width=10,height=10)


plt<-ggplot(emt_and_d_scores, aes(y=Wu_DScores, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")

ggsave(plt,file="dscore_perclone_boxplot.pdf",width=10,height=10)


#check esr1 expression on clones
gene_expression<-FetchData(object = dat, vars = c("ESR1","KRT5"), assay="SCT")
chromatin_expression<-FetchData(object = dat, vars = c("MA0112.3"), assay="chromvar")
dat_tmp<-AddMetaData(dat,gene_expression)
dat_tmp<-AddMetaData(dat_tmp,chromatin_expression,col.name="ESR1_TF")

esr1_clones<-dat_tmp@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
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

plt3<-ggplot(esr1_clones, aes(y=ESR1_TF, x=merged_assay_clones,fill=assigned_celltype,color=assigned_celltype)) + 
  geom_boxplot(width=1,outlier.shape = NA,fill=NA,color="black") +
  geom_jitter(width = 0.5) +theme(axis.text.x = element_text(angle = 90))+

  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")
ggsave(plt1/plt2/plt3,file="ESR1_perclone_boxplot.pdf",width=10,height=10)



#per clone, 50 cells per clone minimum
scsubtype_freq<-dat@meta.data %>% 
  filter(!isNA(merged_assay_clones) & !isNA(scsubtype)) %>% 
  filter(merged_assay_clones %in% clone_filter) %>%
  group_by(sample,Diagnosis,Mol_Diagnosis,scsubtype,merged_assay_clones) %>% 
  count(scsubtype,.drop=TRUE)


plt<-ggplot(scsubtype_freq, aes(fill=factor(scsubtype,levels=names(scsubtype_col)), y=n, x=merged_assay_clones)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=scsubtype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")+ theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file="scsubtype_perclone_percbarplot.pdf",width=10,height=10)




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

