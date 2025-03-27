#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
library(plyr)
library(optparse)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
set.seed(1234)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="5_merged.geneactivity.SeuratObject.rds", 
              help="Input seurat object", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat<-readRDS(file=opt$object_input)

#clustering function for RNA/ATAC/RNA+ATAC
multimodal_cluster<-function(dat=dat,res=0.5,prefix="allcells"){
  # Perform standard analysis of each modality independently 
  #RNA analysis
  DefaultAssay(dat) <- 'RNA'
  dat<-NormalizeData(dat) %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  dat <- RunUMAP(dat, 
    reduction="pca", 
    dims = 1:30, 
    reduction.name = paste(prefix,"umap","rna",sep="."),
    reduction.key = "rnaUMAP_")

  #ATAC analysis
  DefaultAssay(dat)<- 'ATAC'
  dat<-RunTFIDF(dat) %>%  FindTopFeatures() %>% RunSVD()
  dat <- RunUMAP(dat, 
    reduction = "lsi", 
    dims = 2:30, 
    reduction.name=paste(prefix,"umap","atac",sep="."), 
    reduction.key = "atacUMAP_")

  # build a joint neighbor graph using both assays
  dat <- FindMultiModalNeighbors(object = dat,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:40),
    modality.weight.name = "RNA.weight",
    weighted.nn.name=paste(prefix,"weighted.nn",sep="."),
    snn.graph.name=paste(prefix,"wsnn",sep="."),
    verbose = TRUE)

  dat <- RunUMAP(dat, 
    nn.name = paste(prefix,"weighted.nn",sep="."), 
    reduction.name = paste(prefix,"wnn.umap",sep="."), 
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

dat<-multimodal_cluster(dat)

predicted_id_list<-colnames(dat@meta.data)[endsWith(colnames(dat@meta.data),suffix="predicted.id")]
predicted_id_list<-c("seurat_clusters","Diagnosis","Manuscript_Name","Mol_Diagnosis",predicted_id_list)
umap_plot<-lapply(predicted_id_list,function(x) umap_plotting(dat,metadat_column=x))

plt_out<-patchwork::wrap_plots(umap_plot, ncol = 1)
plt_qc<-FeaturePlot(dat,features=c("scrublet_Scores","nCount_SCT","nCount_ATAC"),reduction = "allcells.wnn.umap",ncol=3)

ggsave(plt_out+plt_qc,file="allcells.umap.pdf",width=40,height=length(predicted_id_list)*10,limitsize=F)

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

ggsave(plt,file="seuratclusters_celltypes.features.pdf",height=10,width=40,limitsize=F)

#just top level of HBCA cell types
#https://navinlabcode.github.io/HumanBreastCellAtlas.github.io/assets/svg/celltype_tree.svg

dat$assigned_celltype<-"cancer"
dat@meta.data[dat$seurat_clusters %in% c("39","12","18","21"),]$assigned_celltype<-"luminal_hs"
dat@meta.data[dat$seurat_clusters %in% c("22","20"),]$assigned_celltype<-"luminal_asp"
dat@meta.data[dat$seurat_clusters %in% c("29","16","28"),]$assigned_celltype<-"basal_myoepithelial"

dat@meta.data[dat$seurat_clusters %in% c("40"),]$assigned_celltype<-"adipocyte"
dat@meta.data[dat$seurat_clusters %in% c("11"),]$assigned_celltype<-"endothelial_vascular"
dat@meta.data[dat$seurat_clusters %in% c("37"),]$assigned_celltype<-"endothelial_lymphatic"
dat@meta.data[dat$seurat_clusters %in% c("33"),]$assigned_celltype<-"pericyte"
dat@meta.data[dat$seurat_clusters %in% c("28","13","8"),]$assigned_celltype<-"fibroblast"

dat@meta.data[dat$seurat_clusters %in% c("7"),]$assigned_celltype<-"myeloid"
dat@meta.data[dat$seurat_clusters %in% c("34"),]$assigned_celltype<-"bcell"
dat@meta.data[dat$seurat_clusters %in% c("24"),]$assigned_celltype<-"plasma"
dat@meta.data[dat$seurat_clusters %in% c("15"),]$assigned_celltype<-"tcell"
dat$seurat_clusters_cellassignment<-dat$seurat_clusters
dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)

####################################################
#           Fig 1 Feature Plot                  #
###################################################

Idents(dat)<-factor(dat$assigned_celltype,levels=c("cancer","luminal_hs","luminal_asp","basal_myoepithelial",
"adipocyte","endothelial_vascular","endothelial_lymphatic","pericyte","fibroblast",
"myeloid","bcell","plasma","tcell"))
plt<-DotPlot(subset(dat,cells=names(Idents(dat))),features=features,cluster.idents=FALSE,dot.scale=8)+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))+
  theme(axis.text.x = element_text(angle=90))

ggsave(plt,file="FIG1_assigned_celltypes.features.pdf",height=10,width=40,limitsize=F)


####################################################
#           Fig 1 Colors                            #
###################################################

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
"bcell"="B cells",
"plasma"="#742b8c",
"tcell"="#003147")


hist_col=c("NAT"="#99CCFF",
"DCIS"="#CCCCCC",
"IDC"="#FF9966",
"ILC"="#006633")

clin_col=c("IDC ER+/PR−/HER2+"="#f37872", 
"DCIS"="#cccccb", 
"IDC ER+/PR−/HER2−"="#7fd0df", 
"IDC ER+/PR+/HER2−"="#8d86c0", 
"ILC ER+/PR−/HER2−"="#b9db98", 
"ILC ER+/PR+/HER2−"="#0f6bea1", 
"NAT NA"="#c1e7f6")

sampled_col=c("Primary"="#8A4C80",
"Metastasis of Lymph Node"="#4c9173","NAT"="#99CCFF")


assay_col=c("0"="white","1"="black")

####################################################
#           Fig 1 UMAP                            #
###################################################

p1<-DimPlot(dat,group.by="seurat_clusters",reduction = "allcells.wnn.umap")
p2<-DimPlot(dat,cols=celltype_col,group.by="assigned_celltype",reduction = "allcells.wnn.umap",col=celltype_col)
p3<-DimPlot(dat,cols=hist_col,group.by="Diagnosis",reduction = "allcells.wnn.umap")
p4<-DimPlot(dat,cols=clin_col,group.by="Diag_MolDiag",reduction = "allcells.wnn.umap")
p5<-DimPlot(dat,group.by="sample",reduction = "allcells.wnn.umap")

ggsave(p1/p2/p3/p4/p5,file="FIG1_umap_assigned_celltype.pdf",width=10,height=50)

####################################################
#           Fig 1 Sample Heatmap                  #
###################################################
met<-dat@meta.data
met<-met[!duplicated(met$sample),]
age=c('DCIS_01'='31', 'DCIS_02'='49', 'DCIS_03'='61', 'IDC_01'='75', 'IDC_10'='68', 'IDC_11'='NA', 'IDC_12'='NA', 'IDC_02'='51', 'IDC_03'='74', 'IDC_04'='67', 'IDC_05'='34', 'IDC_06'='76', 'IDC_07'='44', 'IDC_08'='63', 'IDC_09'='63', 'ILC_01'='57', 'ILC_02'='64','NAT_11'='37', 'NAT_14'='50', 'NAT_04'='67', 'IDC_13'='68', 'IDC_14'='40', 'IDC_15'='43', 'IDC_16'='75', 'ILC_03'='71', 'ILC_04'='65', 'ILC_05'='34')
plot_order=c('DCIS_01'='1', 'DCIS_02'='2', 'DCIS_03'='3', 'IDC_01'='4', 'IDC_02'='5', 'IDC_03'='6', 'IDC_04'='7', 'IDC_16'='9', 'ILC_04'='10', 'IDC_05'='11', 'IDC_06'='12', 'IDC_07'='13', 'IDC_08'='14', 'IDC_09'='15', 'IDC_10'='16', 'IDC_11'='17', 'IDC_12'='18', 'IDC_13'='19', 'IDC_15'='20', 'ILC_01'='21', 'ILC_02'='22', 'ILC_03'='23', 'ILC_05'='24', 'NAT_04'='25', 'NAT_11'='26', 'NAT_14'='27') 
spatial_atac=c('DCIS_01'='0', 'DCIS_02'='0', 'DCIS_03'='0', 'IDC_01'='1', 'IDC_02'='1', 'IDC_03'='0', 'IDC_04'='0', 'IDC_16'='0', 'ILC_04'='1', 'IDC_05'='0', 'IDC_06'='1', 'IDC_07'='1', 'IDC_08'='1', 'IDC_09'='1', 'IDC_10'='0', 'IDC_11'='0', 'IDC_12'='0', 'IDC_13'='0', 'IDC_15'='0', 'ILC_01'='0', 'ILC_02'='0', 'ILC_03'='0', 'ILC_05'='0', 'NAT_04'='0', 'NAT_11'='0', 'NAT_14'='0') 
wgs=c('DCIS_01'='1', 'DCIS_02'='1', 'DCIS_03'='1', 'IDC_01'='1', 'IDC_02'='1', 'IDC_03'='1', 'IDC_04'='1', 'IDC_16'='1', 'ILC_04'='1', 'IDC_05'='0', 'IDC_06'='1', 'IDC_07'='1', 'IDC_08'='1', 'IDC_09'='1', 'IDC_10'='1', 'IDC_11'='1', 'IDC_12'='0', 'IDC_13'='1', 'IDC_15'='1', 'ILC_01'='0', 'ILC_02'='1', 'ILC_03'='1', 'ILC_05'='1', 'NAT_04'='0', 'NAT_11'='1', 'NAT_14'='1') 

met$age<-as.numeric(age[met$sample])
met$wgs<-as.numeric(wgs[met$sample])
met$spatial_atac<-as.numeric(spatial_atac[met$sample])
met$plot_order<-as.numeric(plot_order[met$sample])
met$multiome<-1

sample_heatmap<-met[c("Manuscript_Name","age","sample_weight","Diagnosis","Mol_Diagnosis","sampled_site","plot_order","multiome","wgs","spatial_atac")]
row.names(sample_heatmap)<-sample_heatmap$Manuscript_Name
sample_heatmap$Diag_MolDiag<-paste(sample_heatmap$Diagnosis,sample_heatmap$Mol_Diagnosis)
sample_heatmap<-sample_heatmap[order(sample_heatmap$plot_order),]
sample_heatmap<-sample_heatmap[c("age","sample_weight","Diagnosis","Diag_MolDiag","sampled_site","multiome","wgs","spatial_atac")]
sample_heatmap[startsWith(sample_heatmap$sampled_site,"NAT"),]$sampled_site<-"NAT"
age_col=colorRamp2(breaks=c(min(sample_heatmap$age,na.rm=T),max(sample_heatmap$age,na.rm=T)),c("#d789d7","#2a3d66"))
mass_col=colorRamp2(breaks=c(min(sample_heatmap$sample_weight,na.rm=T),max(sample_heatmap$sample_weight,na.rm=T)),c("#f2fc9f","#b05977"))

#plot metadata
ha = rowAnnotation(age=sample_heatmap$age,
                      mass=sample_heatmap$sample_weight,
                      histological_type=sample_heatmap$Diagnosis,
                      molecular_type=sample_heatmap$Diag_MolDiag,
                      sampled_site=sample_heatmap$sampled_site,
                      col = list(age=age_col,
                                  mass=mass_col,
                                    histological_type =hist_col,
                                    clinical_subtype=clin_col,
                                    sampled_site=sampled_col))

pdf("FIG1_sample_metadata.heatmap.pdf")
plt<-Heatmap(sample_heatmap[c("multiome","wgs","spatial_atac")],
 cluster_columns=F,cluster_rows=F,
 left_annotation=ha,
 col=assay_col)
print(plt)
dev.off()

####################################################
#           Fig 1 Stacked Celltype ID             #
###################################################
#Make stacked barplot on identities per cluster

DF<-as.data.frame(dat@meta.data %>% group_by(assigned_celltype,sample) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=assigned_celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_color_discrete(values=celltype_col)
ggsave(plt1,file="FIG1_allcells.assigned_celltype_barplots.pdf",width=50,limitsize=F)

#plot cellcount (no epi distinction)
DF<-as.data.frame(dat@meta.data %>% group_by(sample) %>% tally())
DF$sample<-factor(DF$sample,levels=names(plot_order))
DF$log_count<-log10(DF$n)
plt1<-ggplot(DF)+geom_bar(aes(x=sample,y=log_count),stat="identity",position="dodge")+theme_minimal()
ggsave(plt1,file="FIG1_allcells.cellcount_barplots.pdf",width=50,limitsize=F)

saveRDS(dat,file="6_merged.celltyping.SeuratObject.rds")
