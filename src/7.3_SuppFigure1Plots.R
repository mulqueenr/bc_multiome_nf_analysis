#ssh mulqueen@arc-infra-1
#screen
#srun --partition=guest --cpus-per-task=30 --time=12:00:00 --mem=300G --nodes=1 --pty /bin/bash
#srun --partition=interactive --cpus-per-task=30 --time=12:00:00 --mem=300G --nodes=1 --pty /bin/bash
#srun --partition=cedar --cpus-per-task=30 --time=12:00:00 --mem=300G --nodes=1 --pty /bin/bash

#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome --bind /home/groups/MohammedLab $sif

library(Signac)
library(Seurat)
library(circlize)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
library(plyr)
library(optparse)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(harmony) #local
library(dendextend)
set.seed(1234)
setwd("/home/groups/MohammedLab/bc_multiome/seurat_object")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="6_merged.celltyping.SeuratObject.rds", 
              help="Input seurat object", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat <- readRDS(file=opt$object_input)
dir.create("/home/groups/MohammedLab/bc_multiome/suppfig1")


####################################################
#           Supp Fig 1 Feature Plot                  #
###################################################

#snRNA markers
#from Kumar et al.
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT")
hbca_snmarkers[["basal"]]=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2") 
hbca_snmarkers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
#hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
hbca_snmarkers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
#hbca_snmarkers[["pdc"]]=c("IGKC","PTGDS","IRF8","DNASE1L3","LGALS2","C1orf54","CLIC3")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
features<-llply(hbca_snmarkers, unlist)

Idents(dat)<-factor(dat$assigned_celltype,levels=c("cancer","luminal_hs","luminal_asp","basal_myoepithelial",
"adipocyte","endothelial_vascular","endothelial_lymphatic","pericyte","fibroblast",
"myeloid","bcell","plasma","tcell"))
plt<-DotPlot(subset(dat,cells=names(Idents(dat))),features=features,cluster.idents=FALSE,dot.scale=8)+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))+
  theme(axis.text.x = element_text(angle=90))

ggsave(plt,file=paste0("/home/groups/MohammedLab/bc_multiome/suppfig1/suppfig1_assigned_celltypes.features.pdf"),height=10,width=40,limitsize=F)


####################################################
#           Supp Fig 1 Heatmap                         #
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
plot_top_tf_markers_nocnv<-function(x=out_subset,colfun,group_by,prefix,n_markers=20,order_by_idents=FALSE,plot_by,outdir){

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
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_rna),row.names(tf_ga)))
    average_matrix=(tf_rna+tf_motif+tf_ga)/3. #matrix averages for clustering
    #average_matrix=tf_motif #just cluster only on cnvs

    #set up heatmap seriation and order by average z score
    #first cluster columns, then sort row orders by average

    o_cols =t(average_matrix) %>% dist()  %>% 
                          hclust() %>%
                          as.dendrogram() %>%
                          ladderize() %>% labels()

    o_rows = average_matrix %>% dist()  %>% 
                          hclust() %>%
                          as.dendrogram() %>%
                          ladderize() %>% labels()

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
    colfun_ga<-colfun
    colfun_rna<-colfun
    colfun_motif<-colfun

    tf_rna<-tf_rna[o_rows,o_cols]
    tf_ga<-tf_ga[o_rows,o_cols]
    tf_motif<-tf_motif[o_rows,o_cols]

    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), 
                                            labels =row.names(tf_rna),
                                            labels_gp=gpar(fontsize=6)),
                            motifs = anno_image(motif_plots))


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
    print(draw(ga_plot+rna_plot+motif_plot,row_title=prefix))
    dev.off()
    print(paste("Plotted... ",paste0(outdir,"/",paste0(prefix,".tf.heatmap.pdf"))))
}



colfun <- colorRamp2(
  #breaks = c(-3,-1,-0.5,0,0.5,1,3),
  breaks = c(-3,-2,-1,0,1,2,3),
  colors = c("#053061","#487590","#d1e5f0","#f7f7f7","#fddbc7","#a76146","#67001f"))

#all cells by cell types
plot_top_tf_markers_nocnv(x=dat,
        group_by="assigned_celltype",
        plot_by="assigned_celltype",
        prefix="celltypes",
        colfun=colfun,
        n_markers=10,
        order_by_idents=TRUE,
        outdir="/home/groups/MohammedLab/bc_multiome/suppfig1")

markers <- presto:::wilcoxauc.Seurat(X = dat, group_by = "assigned_celltype", 
  groups_use=unname(unlist(unique(dat@meta.data$assigned_celltype))),
  y=unname(unlist(unique(dat@meta.data$assigned_celltype))), 
  assay = 'data', seurat_assay = "ATAC")

library(dplyr)
markers<-markers %>% filter(pval<0.05)
saveRDS(markers,file="/home/groups/MohammedLab/bc_multiome/suppfig1/celltypes.markers.peaks.df.rds")
write.table(markers,file="/home/groups/MohammedLab/bc_multiome/suppfig1/celltypes.markers.peaks.df.tsv",col.names=T,sep="\t",quote=F,row.names=F)


celltype_freq<-dat@meta.data %>% 
  group_by(sample,Diagnosis,Mol_Diagnosis,assigned_celltype) %>% 
  count(assigned_celltype,.drop=FALSE)

plt<-ggplot(celltype_freq, aes(fill=factor(assigned_celltype,levels=names(celltype_col)), y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity",width = 1) +
  scale_fill_manual(values=celltype_col)+ 
  facet_grid(~paste(Diagnosis,Mol_Diagnosis),scale="free_x",space="free")+ theme(axis.text.x = element_text(angle = 90))

ggsave(plt,file=paste0(outdir,"/","celltype_percbarplot.pdf"),width=10,height=10)
print(paste0(outdir,"/","celltype_percbarplot.pdf"))



#plot atac markers
outdir="/home/groups/MohammedLab/bc_multiome/suppfig1"
DefaultAssay(dat)<-"ATAC"
atac_markers<-c("FOXA1","ELF5","KRT17","ACACB","VWF","LAMA2","CD163","IL7R")
plt_list<-lapply(atac_markers,function(gene) {
  CoveragePlot(
      object = dat, 
      region = gene, 
      upstream=2000,downstream=2000,
      features = gene, combine=F,
      expression.assay = "RNA",
      idents = Idents(dat))+
    scale_fill_manual(celltype_col)
    })



ggsave(plt,file=paste0(outdir,"/","celltype_atac_markers.pdf"),width=length(atac_markers)*10,height=10)
print(paste0(outdir,"/","celltype_atac_markers.pdf"))
