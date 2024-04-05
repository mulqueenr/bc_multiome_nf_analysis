#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/harmony.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(reshape2)
library(optparse)
library(plyr)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
#opt$object_input="merged.clustered.SeuratObject.rds"


#split by immune, epithelial, stromal
#label transfer with more cell state refinement
#https://www.nature.com/articles/s41588-024-01688-9/figures/2

dat=readRDS(opt$object_input)

immune<-subset(dat,cells=row.names(dat@meta.data)[dat$lineage %in% c("immune")])
stromal<-subset(dat,cells=row.names(dat@meta.data)[dat$lineage %in% c("stromal")])
epithelial<-subset(dat,cells=row.names(dat@meta.data)[dat$lineage %in% c("epithelial")])


label_transfer<-function(dat,ref_obj,ref_prefix){
  ref_obj<-SCTransform(ref_obj, verbose = FALSE)
  transfer.anchors <- FindTransferAnchors(
    reference = ref_obj,
    reference.assay="SCT",
    query = dat,
    query.assay="SCT",
    features=VariableFeatures(ref_obj),
    verbose=T,
    dims=1:30
  )

  predictions<- TransferData(
    anchorset = transfer.anchors,
    refdata = ref_obj$cell_type,
  )
  colnames(predictions)<-paste0(ref_prefix,"_",colnames(predictions))

  dat<-AddMetaData(dat,metadata=predictions)
  
  return(dat)
  }

cell_lineage_clustering<-function(x,outprefix,res,markers,ref){
  x<-label_transfer(dat=x,ref_obj=ref,ref_prefix="reed")
  DefaultAssay(x)<-"peaks"
  x <- RunTFIDF(x)
  x <- FindTopFeatures(x, min.cutoff = 50)
  x <- RunSVD(x, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
  x <- RunUMAP(x, reduction = 'lsi', dims = 2:30)
  p1 <- DimPlot(x, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Unintegrated ATAC")

  x <- RunHarmony(
    object = x,
    group.by.vars = 'sample',
    reduction.use = 'lsi',
    assay.use = 'peaks',
    reduction.save= 'harmony_atac',
    project.dim = FALSE, ncores=10
  )

  # re-compute the UMAP using corrected LSI embeddings
  x <- RunUMAP(x, dims = 2:30, reduction = 'harmony_atac')
  p2 <- DimPlot(x, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Harmony ATAC")

  DefaultAssay(x)<-"SCT"
  x <- RunPCA(x)
  x <- RunUMAP(x, reduction = 'pca', dims = 1:30)
  p3 <- DimPlot(x, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Unintegrated RNA")

  x <- RunHarmony(
    object = x,
    group.by.vars = 'sample',
    reduction.use = 'pca',
    assay.use = 'SCT',
    reduction.save= 'harmony_rna',
    project.dim = FALSE, ncores=10
  )

  # re-compute the UMAP using corrected LSI embeddings
  x <- RunUMAP(x, dims = 1:30, reduction = 'harmony_rna')
  p4 <- DimPlot(x, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Harmony RNA")


  plt<-(p1 + p2)/(p3+p4)
  ggsave(plt,file=paste0(outprefix,"_harmony_integrations.pdf"),height=20,width=20)

  x <- FindMultiModalNeighbors(
  object = x,
  reduction.list = list("harmony_rna", "harmony_atac"),
  dims.list = list(1:50, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
  )


  x <- RunUMAP(
  object = x,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  assay = "RNA",
  verbose = TRUE
  )

  x$diag_moldiag<-paste(x$Diagnosis,x$Mol_Diagnosis,sep="_")

  x<- FindNeighbors(x)
  x<- FindClusters(x,graph.name='wsnn',resolution=res)

  plt2<-DimPlot(x,reduction = "wnn.umap", group.by = c('sample','Diagnosis','HBCA_predicted.id','diag_moldiag','seurat_clusters','reed_predicted.id'),ncol=2)

  ggsave(plt2,file=paste0(outprefix,"_harmony_integration.coembedded.pdf"),width=20,height=30)
  Idents(x)<-x$seurat_clusters
  plt1<-DotPlot(dat,assay = "SCT",markers, dot.scale = 10,cluster.idents = TRUE)+ RotatedAxis()+ scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))

  ggsave(plt1,file=paste0(outprefix,"_markers.pdf"),height=10,width=30,limitsize = FALSE)

  return(x)
}


#markers for cell states from https://www.nature.com/articles/s41588-024-01688-9/figures/2
immune_markers=list()
immune_markers[["immune"]]=c("PTPRC","CXCR4")
immune_markers[["tcell_lineage"]]=c("CD3D","CD3G")
immune_markers[["cd4_all"]]=c("CD4","CCR4")
immune_markers[["cd4_naive"]]=c("SELL","CD40LG")
immune_markers[["cd4_th"]]=c("CCL2O","PTPN13")
immune_markers[["cd8_t_all"]]=c("CD8A","CD8B")
immune_markers[["cd8_tem"]]=c("GZMK","TNFRSF9")
immune_markers[["cd8_trm"]]=c("KLRC1","ITGA1")
immune_markers[["cd8_tc1"]]=c("IFNG","TNF")
immune_markers[["nk_nkt"]]=c("GZMH","GZMB")
immune_markers[["nk"]]=c("KLRC3","XCL2")
immune_markers[["ilc"]]=c("IL1R1","KIT")
immune_markers[["bcell"]]=c("MS4A1","CD79A")
immune_markers[["b_naive"]]=c("TCL1A","YBX3")
immune_markers[["b_mem_switched"]]=c("LTB","TNFRSF13C")
immune_markers[["b_mem_unswitched"]]=c("MYC","EGR3")
immune_markers[["plasma"]]=c("IGKC","JCHAIN")
immune_markers[["macrophage"]]=c("FCER1G","C1QB")
immune_markers[["macrophage_lipid_assoc"]]=c("APOE","GPNMB")
immune_markers[["dendritic"]]=c("CPVL","ID2")
immune_markers<-llply(immune_markers, unlist)
reed_immune<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/ref/reed/immune.rds")
immune<-cell_lineage_clustering(immune,outprefix="immune",res=0.05,markers=immune_markers,ref=reed_immune)

stromal_markers=list()
stromal_markers[["fibro_all"]]=c("DCN","LUM")
stromal_markers[["fibro_1"]]=c("COL1A2","VEGFD")
stromal_markers[["fibro_2"]]=c("ARSG","MMP3")
stromal_markers[["fibro_3"]]=c("GREM1","CLU")
stromal_markers[["fibro_4"]]=c("MTSS1","SLC2A1")
stromal_markers[["perivasc"]]=c("PDGFRB","NOTCH3","MYH11")
stromal_markers[["perivasc_1"]]=c("ADRA2C","CTSC")
stromal_markers[["perivasc_2"]]=c("GPAT2","PGAP1")
stromal_markers[["perivasc_3"]]=c("RERGL","ACTA2")
stromal_markers[["perivasc_4"]]=c("ADGRL3","GGT5")
stromal_markers[["perivasc_5"]]=c("RGS5","THBS4")
stromal_markers[["vasc_endo"]]=c("PLVAP","PECAM1","CLDN5")
stromal_markers[["vasc_endo_vein"]]=c("ACKR1","SELP")
stromal_markers[["vasc_endo_capillary"]]=c("CA4","CD36")
stromal_markers[["vasc_endo_artery"]]=c("SOX17","EFNB2")
stromal_markers[["vasc_endo_angioterm"]]=c("PXDN","ANGPT2")
stromal_markers[["lymph_endo"]]=c("MMRN1","PDPN","CCL21")
stromal_markers[["lymph_endo_1"]]=c("SEMA3D","B3GNT7")
stromal_markers[["lymph_endo_2"]]=c("LYVE1","TXB1")
stromal_markers<-llply(stromal_markers, unlist)
reed_stromal<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/ref/reed/stroma.rds")
cell_lineage_clustering(stromal,outprefix="stromal",res=0.05,markers=stromal_markers,ref=reed_immune)


epithelial_markers=list()
epithelial_markers[["epi_all"]]=c("EPCAM","CDH1","KRT17")
epithelial_markers[["lactocyte"]]=c("CSN2","CSN3","LALBA")
epithelial_markers[["LASP"]]=c("KIT","ALDH1A3","SLPI") #lumSec progenitors
epithelial_markers[["LASP_1"]]=c("SERPINB3","MUC16") #lumSec
epithelial_markers[["LASP_2"]]=c("KRT16","KRT6B") #lumSec
epithelial_markers[["LASP_3"]]=c("RASSF1","COBL") #lumSec
epithelial_markers[["LASP_4"]]=c("ROPN1","ELF5") #lumSec
epithelial_markers[["LASP_5"]]=c("MKI67","TOP2A","BRCA1","BRCA2") #lumSec
epithelial_markers[["LHS"]]=c("AREG","ANKRD30A","FOXA1") #lumHR
epithelial_markers[["LHS_1"]]=c("TMEM45B","CACNG4") #lumHR
epithelial_markers[["LHS_2"]]=c("ESR1","PGR") #lumHR
epithelial_markers[["LHS_3"]]=c("SERPINA1","PIP") #lumHR
epithelial_markers[["BMYO"]]=c("ACTG2","ACTA2","TAGLN") #basal
epithelial_markers[["BMYO_1"]]=c("KRT14","KRT5") #basal
epithelial_markers[["BMYO_2"]]=c("HSPA1B","OXTR") #basal
reed_stromal<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/ref/reed/epithelialrds")
cell_lineage_clustering(epithelial,outprefix="epithelial",res=0.05,markers=epithelial_markers,ref=reed_epithelial)


#########NONEPI###########
