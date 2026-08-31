#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

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
set.seed(1234)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="5_merged.geneactivity.SeuratObject.rds", 
              help="Input seurat object", metavar="character")
); 
 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat<-readRDS(file=opt$object_input)

if (!dir.exists("7_celltyping")) {
  dir.create("7_celltyping")
}

#clustering function for RNA/ATAC/RNA+ATAC
multimodal_cluster<-function(dat=dat,res=0.5,prefix="allcells",rna_pcs=1:50,atac_pcs=50,harmony_integrate=FALSE){
  # Perform standard analysis of each modality independently 
  # if harmony_integrate is set to TRUE run harmony on each modality then merge via multimodal neighbors
    #https://github.com/satijalab/seurat/issues/6094

  #RNA analysis
  DefaultAssay(dat) <- 'RNA'

  dat<-NormalizeData(dat) %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs=max(rna_pcs))
  
  if(harmony_integrate){
    #integrate PCA over layers (samples)
    dat <- IntegrateLayers(object = dat, assay="RNA", method = HarmonyIntegration, 
                            orig.reduction = "pca",
                            new.reduction = 'pca', verbose = TRUE)
    }

  dat <- RunUMAP(dat, 
    reduction="pca", 
    dims = rna_pcs, 
    reduction.name = paste(prefix,"umap","rna",sep="."),
    reduction.key = "rnaUMAP_")

  #ATAC analysis
  DefaultAssay(dat)<- 'ATAC'
  dat <- RunTFIDF(dat) %>%  FindTopFeatures() %>% RunSVD(n=atac_pcs)

  if(harmony_integrate){
      dat<-RunHarmony(
            object = dat,
            group.by.vars = 'sample',
            reduction = 'lsi',
            assay.use = 'ATAC',
            project.dim = FALSE,
            reduction.save = "lsi")
    }

  dat <- RunUMAP(dat, 
    reduction = "lsi", 
    dims = c(2:atac_pcs), 
    reduction.name=paste(prefix,"umap","atac",sep="."), 
    reduction.key = "atacUMAP_")

  # build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(object = dat,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:max(rna_pcs), 2:max(atac_pcs)),
    modality.weight.name = "RNA.weight",
    weighted.nn.name=paste(prefix,"weighted.nn",sep="."),
    snn.graph.name=paste(prefix,"wsnn",sep="."),
    verbose = TRUE)

  dat <- RunUMAP(dat, 
    nn.name = paste(prefix,"weighted.nn",sep="."), 
    reduction.name = paste(prefix,"wnn.umap",sep="."), 
    n.neighbors=200, ###UP THIS FOR LESS STRANDYNESS IN UMAP
    n.epochs=200,
    min.dist=0.1,
    reduction.key = "wnnUMAP_")
  
  dat <- FindClusters(dat, 
    graph.name = paste(prefix,"wsnn",sep="."), 
    algorithm = 3, 
    resolution = res, 
    verbose = FALSE)


 return(dat)
}

umap_plotting<-function(dat,metadat_column,prefix="allcells",dotsize=1){
  p1<-DimPlot(dat, pt.size=dotsize,group.by=metadat_column,label = TRUE, repel = TRUE, reduction = paste(prefix,"umap.rna",sep="."),raster=T) + ggtitle(paste(metadat_column,prefix,"rna"))
  p2<-DimPlot(dat, pt.size=dotsize,group.by=metadat_column,label = TRUE, repel = TRUE, reduction = paste(prefix,"umap.atac",sep="."),raster=T) + ggtitle(paste(metadat_column,prefix,"atac"))
  p3<-DimPlot(dat, pt.size=dotsize,group.by=metadat_column,label = TRUE, repel = TRUE, reduction = paste(prefix,"wnn.umap",sep="."),raster=T) + ggtitle(paste(metadat_column,prefix,"rna+atac"))
return(p1|p2|p3)
}

dat<-multimodal_cluster(dat,rna_pcs=1:30,atac_pcs=30,res=0.8)

predicted_id_list<-colnames(dat@meta.data)[endsWith(colnames(dat@meta.data),suffix="predicted.id")]
predicted_id_list<-c("seurat_clusters","Diagnosis","Manuscript_Name","Mol_Diagnosis",predicted_id_list)
umap_plot<-lapply(predicted_id_list,function(x) umap_plotting(dat,metadat_column=x))

plt_out<-patchwork::wrap_plots(umap_plot, ncol = 1)
plt_qc<-FeaturePlot(dat,features=c("scrublet_Scores","nCount_SCT","nCount_ATAC"),reduction = "allcells.wnn.umap",ncol=3)

ggsave(plt_out+plt_qc,file=paste0(getwd(),"/7_celltyping/","allcells.umap.pdf"),width=40,height=length(predicted_id_list)*10,limitsize=F)


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

DefaultAssay(dat)<-"RNA"
Idents(dat)<-dat$seurat_clusters
plt<-DotPlot(subset(dat,cells=names(Idents(dat))),features=features,cluster.idents=TRUE,dot.scale=8)+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))+
  theme(axis.text.x = element_text(angle=90))

ggsave(plt,file=paste0("./7_celltyping/","seuratclusters_celltypes.features.pdf"),height=10,width=40,limitsize=F)
plt_cluster<-DimPlot(dat,group.by="seurat_clusters",reduction = "allcells.wnn.umap",label=TRUE)
ggsave(plt_cluster,file=paste0(getwd(),"/7_celltyping/","seuratclusters_celltypes.dimplot.pdf"),height=10,width=10,limitsize=F)
paste0(getwd(),"/7_celltyping/","seuratclusters_celltypes.dimplot.pdf")
#just top level of HBCA cell types
#https://navinlabcode.github.io/HumanBreastCellAtlas.github.io/assets/svg/celltype_tree.svg

#hard coded because of seed setting
dat$assigned_celltype<-"cancer"
dat@meta.data[dat$seurat_clusters %in% c("39","12","18","21"),]$assigned_celltype<-"luminal_hs"
dat@meta.data[dat$seurat_clusters %in% c("22","20"),]$assigned_celltype<-"luminal_asp"
dat@meta.data[dat$seurat_clusters %in% c("29","16","28"),]$assigned_celltype<-"basal_myoepithelial"

dat@meta.data[dat$seurat_clusters %in% c("40"),]$assigned_celltype<-"adipocyte"
dat@meta.data[dat$seurat_clusters %in% c("11"),]$assigned_celltype<-"endothelial_vascular"
dat@meta.data[dat$seurat_clusters %in% c("37"),]$assigned_celltype<-"endothelial_lymphatic"
dat@meta.data[dat$seurat_clusters %in% c("33"),]$assigned_celltype<-"pericyte"
dat@meta.data[dat$seurat_clusters %in% c("13","8"),]$assigned_celltype<-"fibroblast"

dat@meta.data[dat$seurat_clusters %in% c("7"),]$assigned_celltype<-"myeloid"
dat@meta.data[dat$seurat_clusters %in% c("34"),]$assigned_celltype<-"bcell"
dat@meta.data[dat$seurat_clusters %in% c("24"),]$assigned_celltype<-"plasma"
dat@meta.data[dat$seurat_clusters %in% c("15"),]$assigned_celltype<-"tcell"
dat$seurat_clusters_cellassignment<-dat$seurat_clusters
dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)

saveRDS(dat,file="6_merged.celltyping.SeuratObject.rds")
