#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggalluvial)
library(reshape2)
library(optparse)


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.SeuratObject.rds"

###########PANEL B UPDATE ################
dat<-readRDS(opt$object_input)

# Perform standard analysis of each modality independently 
#RNA analysis
dat[["RNA"]] <- as(dat[["RNA"]], Class = "Assay5")
DefaultAssay(dat) <- 'RNA'
dat<-NormalizeData(dat) %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA()
dat <- RunUMAP(dat, reduction="pca", dims = 1:30, reduction.name = "umap.rna",reduction.key = "rnaUMAP_")

#ATAC analysis
DefaultAssay(dat)<- 'peaks'
dat<-RunTFIDF(dat) %>%  FindTopFeatures() %>% RunSVD()
dat <- RunUMAP(dat, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
dat <- FindMultiModalNeighbors(
  object = dat,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

dat <- RunUMAP(dat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)


p1<-DimPlot(dat, group.by='HBCA_predicted.id',label = TRUE, repel = TRUE, reduction = "umap.rna") + NoLegend()
p2<-DimPlot(dat, group.by='HBCA_predicted.id',label = TRUE, repel = TRUE, reduction = "umap.atac") + NoLegend()
p3<-DimPlot(dat, group.by='HBCA_predicted.id',label = TRUE, repel = TRUE, reduction = "wnn.umap") + NoLegend()

p4<-DimPlot(dat, group.by = 'seurat_clusters',label = TRUE, repel = TRUE, reduction = "umap.rna") + NoLegend()
p5<-DimPlot(dat, group.by = 'seurat_clusters',label = TRUE, repel = TRUE, reduction = "umap.atac") + NoLegend()
p6<-DimPlot(dat, group.by = 'seurat_clusters',label = TRUE, repel = TRUE, reduction = "wnn.umap") + NoLegend()


plt<-(p1 + p2 + p3)/(p4 + p5 + p6)
ggsave(plt,file="umap.pdf",height=20,width=30)


############PANEL B UMAP #########################
dat=readRDS(opt$object_input)
DefaultAssay(dat)<-"peaks"
dat <- RunTFIDF(dat)
dat <- FindTopFeatures(dat, min.cutoff = 50)
dat <- RunSVD(dat, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
dat <- RunUMAP(dat, reduction = 'lsi', dims = 2:30)
p1 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Unintegrated ATAC")

dat <- RunHarmony(
  object = dat,
  group.by.vars = 'sample',
  reduction.use = 'lsi',
  assay.use = 'peaks',
  reduction.save= 'harmony_atac',
  project.dim = FALSE, ncores=10
)

# re-compute the UMAP using corrected LSI embeddings
dat <- RunUMAP(dat, dims = 2:30, reduction = 'harmony_atac')
p2 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Harmony ATAC")

DefaultAssay(dat)<-"SCT"
dat <- RunPCA(dat)
dat <- RunUMAP(dat, reduction = 'pca', dims = 1:50)
p3 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Unintegrated RNA")

dat <- RunHarmony(
  object = dat,
  group.by.vars = 'sample',
  reduction.use = 'pca',
  assay.use = 'RNA',
  reduction.save= 'harmony_rna',
  project.dim = FALSE, ncores=10
)

# re-compute the UMAP using corrected LSI embeddings
dat <- RunUMAP(dat, dims = 1:30, reduction = 'harmony_rna')
p4 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Harmony RNA")


plt<-(p1 + p2)/(p3+p4)
ggsave(plt,file="harmony_integrations.pdf",height=20,width=20)


dat <- FindMultiModalNeighbors(
object = dat,
reduction.list = list("harmony_rna", "harmony_atac"),
dims.list = list(1:50, 2:30),
modality.weight.name = "SCT.weight",
verbose = TRUE
)

dat <- RunUMAP(
object = dat,
nn.name = "weighted.nn",
reduction.name = "wnn.umap",
assay = "RNA",
verbose = TRUE
)

#generate nonharmony integrated multiomics clustering
dat <- FindMultiModalNeighbors(
object = dat,
reduction.list = list("pca", "lsi"),
dims.list = list(1:50, 2:30),
modality.weight.name = "RNA.weight",
weighted.nn.name = "weighted.nn.unintegrated",
verbose = TRUE
)

dat <- RunUMAP(
object = dat,
nn.name = "weighted.nn.unintegrated",
reduction.name = "wnn.umap.unintegrated",
assay = "RNA",
verbose = TRUE
)

dat$diag_moldiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis,sep="_")
plt2<-DimPlot(dat,reduction = "wnn.umap", group.by = c('sample','Diagnosis','HBCA_predicted.id','diag_moldiag'),raster=FALSE)
ggsave(plt2,file="harmony_integration.coembedded.pdf",width=20,height=20)

plt3<-DimPlot(dat,reduction = "wnn.umap.unintegrated", group.by = c('sample','Diagnosis','HBCA_predicted.id','diag_moldiag'),raster=FALSE)
ggsave(plt3,file="multiome_unintegrated.coembedded.pdf",width=20,height=20)
#ended up using multiome_unintegrated.coembedded.pdf for final figure.

dat <- FindNeighbors(dat,reduction='wnn.umap',dims=1:2,graph.name="umap.snn")
dat <- FindClusters(dat,resolution=0.1,graph="umap.snn")

dat$diag_moldiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis,sep="_")
dat$low_res_all_cells_cluster<-dat$seurat_clusters
dat$epithelial<-"FALSE"
dat@meta.data[dat$low_res_all_cells_cluster %in% c("0","1","2","3","4","5","7","8","9","10","13","14","17","19"),]$epithelial<-"TRUE"


#refine cell typing with marker genes
#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT")
hbca_snmarkers[["basal"]]=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2")
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
features<-llply(hbca_snmarkers, unlist)

#Idents(dat)<-factor(dat$EMBO_predicted.id,levels=rev(c("epithelial","cycling.epithelial","CAFs","Endothelial","Pericytes","Myeloid","TAMs","TAMs_2","Plasma.cells","B.cells","T.cells","NA")))
#Idents(dat)<-factor(dat$low_res_all_cells_cluster,levels=c("luminal epithelial cell of mammary gland","basal cell","fibroblast","endothelial cell of lymphatic vessel", "endothelial cell of vascular tree","pericyte","myeloid cell","T cell","mast cell","adipocyte of breast" ))
Idents(dat)<-dat$low_res_all_cells_cluster
plt1<-DotPlot(dat,assay = "SCT",features, dot.scale = 10,cluster.idents = TRUE)+ RotatedAxis()+ scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))

ggsave(plt1,file="hbca_snRNA_markers.low_res.pdf",height=10,width=30,limitsize = FALSE)

met<-as.data.frame(dat@meta.data)
met_celltypes<-melt(as.data.frame(table(
  met$HBCA_predicted.id,
  met$low_res_all_cells_cluster)))
met_celltypes$Var1<-gsub(met_celltypes$Var1,pattern="\\.",replace="-")

hbca_order<-c("luminal epithelial cell of mammary gland","basal cell",
  "fibroblast","pericyte","adipocyte of breast",
  "endothelial cell of vascular tree","endothelial cell of lymphatic vessel",
  "T cell","myeloid cell","mast cell")

cluster_order<-c("0","1","2","3","4","5","7","8","9","10","13","14","17","19","16","6","15","12","18","11")

met_celltypes$Var1<-factor(met_celltypes$Var1,levels=hbca_order)
met_celltypes$Var2<-factor(met_celltypes$Var2,levels=cluster_order)

plt1<-ggplot(as.data.frame(met_celltypes),
       aes(y = value, axis1 = Var1, axis2 = Var2)) +
  geom_alluvium(aes(fill = Var1), width = 1/12) +
  geom_stratum(width = 1/12, aes(fill = Var1)) +
  ggrepel::geom_text_repel(
    aes(label =Var1),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  ) +
  ggrepel::geom_text_repel(
    aes(label = Var2),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5
  ) +
  scale_x_discrete(limits = c("Var1", "Var2"), expand = c(.05, .05)) 
ggsave(plt1,file="cell_type_assignment.alluvial.low_res.pdf",width=25,limitsize=F)
#cluster with/without epithelial cells

dat$lineage<-"epithelial"
dat@meta.data[dat$low_res_all_cells_cluster %in% c("6","15"),]$lineage<-"stromal"
dat@meta.data[dat$low_res_all_cells_cluster %in% c("12","18","11"),]$lineage<-"immune"


plt2<-DimPlot(dat,reduction = "wnn.umap", group.by = c('sample','seurat_clusters','HBCA_predicted.id','diag_moldiag','epithelial','lineage'),raster=FALSE,ncol=2)

ggsave(plt2,file="harmony_integration.coembedded.pdf",width=20,height=30)

saveRDS(dat,file="merged.clustered.SeuratObject.rds")

#########NONEPI###########

nonepi<-subset(dat,cell=row.names(dat@meta.data)[dat$epithelial=="FALSE"])
DefaultAssay(nonepi)<-"peaks"
nonepi <- RunTFIDF(nonepi)
nonepi <- FindTopFeatures(nonepi, min.cutoff = 50)
nonepi <- RunSVD(nonepi, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
nonepi <- RunUMAP(nonepi, reduction = 'lsi', dims = 2:30)
p1 <- DimPlot(nonepi, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Unintegrated ATAC")

nonepi <- RunHarmony(
  object = nonepi,
  group.by.vars = 'sample',
  reduction.use = 'lsi',
  assay.use = 'peaks',
  reduction.save= 'harmony_atac',
  project.dim = FALSE, ncores=10
)

# re-compute the UMAP using corrected LSI embeddings
nonepi <- RunUMAP(nonepi, dims = 2:30, reduction = 'harmony_atac')
p2 <- DimPlot(nonepi, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Harmony ATAC")

DefaultAssay(nonepi)<-"SCT"
nonepi <- RunPCA(nonepi)
nonepi <- RunUMAP(nonepi, reduction = 'pca', dims = 1:50)
p3 <- DimPlot(nonepi, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Unintegrated RNA")

nonepi <- RunHarmony(
  object = nonepi,
  group.by.vars = 'sample',
  reduction.use = 'pca',
  assay.use = 'RNA',
  reduction.save= 'harmony_rna',
  project.dim = FALSE, ncores=10
)

# re-compute the UMAP using corrected LSI embeddings
nonepi <- RunUMAP(nonepi, dims = 1:30, reduction = 'harmony_rna')
p4 <- DimPlot(nonepi, group.by = 'HBCA_predicted.id', pt.size = 0.5) + ggplot2::ggtitle("Harmony RNA")


plt<-(p1 + p2)/(p3+p4)
ggsave(plt,file="harmony_integrations.nonepi.pdf",height=20,width=20)


nonepi <- FindMultiModalNeighbors(
object = nonepi,
reduction.list = list("harmony_rna", "harmony_atac"),
dims.list = list(1:50, 2:30),
modality.weight.name = "SCT.weight",
verbose = TRUE
)


nonepi <- RunUMAP(
object = nonepi,
nn.name = "weighted.nn",
reduction.name = "wnn.umap",
assay = "RNA",
verbose = TRUE
)


#generate nonharmony integrated multiomics clustering
nonepi <- FindMultiModalNeighbors(
object = nonepi,
reduction.list = list("pca", "lsi"),
dims.list = list(1:50, 2:30),
modality.weight.name = "RNA.weight",
weighted.nn.name = "weighted.nn.unintegrated",
verbose = TRUE
)

nonepi <- RunUMAP(
object = nonepi,
nn.name = "weighted.nn.unintegrated",
reduction.name = "wnn.umap.unintegrated",
assay = "RNA",
verbose = TRUE
)


nonepi$diag_moldiag<-paste(nonepi$Diagnosis,nonepi$Mol_Diagnosis,sep="_")
plt2<-DimPlot(nonepi,reduction = "wnn.umap", group.by = c('sample','Diagnosis','HBCA_predicted.id','diag_moldiag'))
ggsave(plt2,file="harmony_integration.coembedded.nonepi.pdf",width=20,height=20)


plt3<-DimPlot(nonepi,reduction = "wnn.umap.unintegrated", group.by = c('sample','Diagnosis','HBCA_predicted.id','diag_moldiag'),raster=FALSE)
ggsave(plt3,file="multiome_unintegrated.coembedded.nonepi.pdf",width=20,height=20)
#ended up using multiome_unintegrated.coembedded.pdf for final figure.

############PANEL C ALLUVIAL PLOTS#########################

embo_cell_cols<-c("epithelial"="#DC3977","T.cells"="#003147","TAMs"="#E9E29C","Plasma.cells"="#B7E6A5","CAFs"="#E31A1C","B.cells"="#089099","NA"="grey","Endothelial"="#EEB479", "Pericytes"= "#F2ACCA", "TAMs_2"="#e9e29c","cycling.epithelial"="#591a32", "Myeloid"="#dbc712")    

met<-as.data.frame(dat@meta.data)
met_celltypes<-melt(as.data.frame(table(met$EMBO_predicted.id,
  met$HBCA_predicted.id,
  met$swarbrick_predicted.id,
  met$seurat_clusters)))
met_celltypes$Var1<-gsub(met_celltypes$Var1,pattern="\\.",replace="-")

met_celltypes$Var1<-paste("pal_",met_celltypes$Var1)
met_celltypes$Var2<-paste("hbca_",met_celltypes$Var2)
met_celltypes$Var3<-paste("wu_",met_celltypes$Var3)

pal_order<-paste("pal_",c("cycling epithelial","epithelial",
  "NA","CAFs","Pericytes","Endothelial",
  "T cells","B cells","Plasma cells","Myeloid","TAMs_2","TAMs"))
hbca_order<-paste("hbca_",c("luminal epithelial cell of mammary gland","basal cell",
  "adipocyte of breast","fibroblast","pericyte",
  "endothelial cell of vascular tree","endothelial cell of lymphatic vessel",
  "T cell","myeloid cell","mast cell"))
wu_order<-paste("wu_",c("Cancer Epithelial","Normal Epithelial",
  "CAFs","PVL","Endothelial",
  "T-cells","B-cells","Plasmablasts","Myeloid"))
met_celltypes$Var1<-factor(met_celltypes$Var1,levels=pal_order)
met_celltypes$Var2<-factor(met_celltypes$Var2,levels=hbca_order)
met_celltypes$Var3<-factor(met_celltypes$Var3,levels=wu_order)
plt1<-ggplot(as.data.frame(met_celltypes),
       aes(y = value, axis1 = Var4, axis2 = Var2)) +
  geom_alluvium(aes(fill = Var2), width = 1/12) +
  geom_stratum(width = 1/12, aes(fill = Var2)) +
  ggrepel::geom_text_repel(
    aes(label =Var4),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  ) +
  ggrepel::geom_text_repel(
    aes(label = Var2),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5
  ) +
  scale_x_discrete(limits = c("Var4", "Var2"), expand = c(.05, .05)) 
ggsave(plt1,file="cell_type_assignment.alluvial.pdf",width=25,limitsize=F)
#cluster with/without epithelial cells


############PANEL D STACKED BAR PLOTS#########################
library(dplyr) 


###########Color Schema#################
type_cols<-c(
#epithelial
"Cancer Epithelial" = "#7C1D6F", "Normal Epithelial" = "#DC3977", #immune
"B-cells" ="#089099", "T-cells" ="#003147","Myeloid" ="#E9E29C", "Plasmablasts"="#B7E6A5", #other
"CAFs" ="#E31A1C", "Endothelial"="#EEB479",  "PVL" ="#F2ACCA")

embo_cell_cols<-c("epithelial"="#DC3977","T.cells"="#003147","TAMs"="#E9E29C","Plasma.cells"="#B7E6A5","CAFs"="#E31A1C","B.cells"="#089099","NA"="grey","Endothelial"="#EEB479", "Pericytes"= "#F2ACCA", "TAMs_2"="#e9e29c","cycling.epithelial"="#591a32", "Myeloid"="#dbc712")    
       
diag_cols<-c("IDC"="red", "DCIS"="grey")

molecular_type_cols<-c("DCIS"="grey", "er+_pr+_her2-"="#EBC258", "er+_pr-_her2-"="#F7B7BB")
########################################

#Set up metadata and set up facet labels as factors for ordering
metadat<-as.data.frame(dat@meta.data)
hbca_order<-c("luminal epithelial cell of mammary gland","basal cell",
  "adipocyte of breast","fibroblast","pericyte",
  "endothelial cell of vascular tree","endothelial cell of lymphatic vessel",
  "T cell","myeloid cell","mast cell")
metadat$HBCA_predicted.id<-factor(metadat$HBCA_predicted.id,levels=hbca_order)

metadat$diagnosis = factor(metadat$Diagnosis, levels=c("NAT","DCIS","IDC","ILC"), labels=c("NAT","DCIS","IDC","ILC")) 
metadat$molecular_type = factor(metadat$Mol_Diagnosis, levels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2+","ER+/PR-/HER2-"), labels=c("NA","DCIS","ER+/PR+/HER2-","ER+/PR-/HER2+","ER+/PR-/HER2-")) 

#Cells PF
metadat$epi<-"Nonepi"
metadat[metadat$HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"),]$epi<-"Epi"
DF<-as.data.frame(metadat %>% group_by(Diagnosis, Mol_Diagnosis,sample,epi) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=epi,y=n))+geom_bar(stat="identity")+theme_minimal()+facet_grid(.~Diagnosis+Mol_Diagnosis,scales="free_x",space="free") #+ scale_y_continuous(trans='log10')
ggsave(plt1,file="barplot_qc_cellcount.pdf")

#Cell types (stacked bar)
DF<-as.data.frame(metadat %>% group_by(Diagnosis, Mol_Diagnosis,sample,HBCA_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=HBCA_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~Diagnosis ~ Mol_Diagnosis,scales="free_x",space="free")
ggsave(plt1,file="hbca_barplot_qc_celltype.pdf")

#Cell types (stacked bar)
DF<-as.data.frame(metadat[metadat$epi=="Nonepi",] %>% group_by(Diagnosis, Mol_Diagnosis,sample,HBCA_predicted.id) %>% tally())
plt1<-ggplot(DF,aes(x=sample,fill=HBCA_predicted.id,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+facet_grid(.~Diagnosis ~ Mol_Diagnosis,scales="free_x",space="free")
ggsave(plt1,file="hbca_barplot_qc_celltype_nonepi.pdf")


############PANEL E MARKER GENES#########################
library(plyr)
library(dplyr) 

#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4")
hbca_snmarkers[["lumsec"]]=c("AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT")
hbca_snmarkers[["basal"]]=c("AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2")
hbca_snmarkers[["fibro"]]=c("LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
hbca_snmarkers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
hbca_snmarkers[["myeloid"]]=c("F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247")
hbca_snmarkers[["mast"]]=c("NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC")
hbca_snmarkers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
features<-llply(hbca_snmarkers, unlist)

#Idents(dat)<-factor(dat$EMBO_predicted.id,levels=rev(c("epithelial","cycling.epithelial","CAFs","Endothelial","Pericytes","Myeloid","TAMs","TAMs_2","Plasma.cells","B.cells","T.cells","NA")))
Idents(dat)<-factor(dat$HBCA_predicted.id,levels=c("luminal epithelial cell of mammary gland",
  "basal cell","fibroblast","endothelial cell of lymphatic vessel",
  "endothelial cell of vascular tree","pericyte","myeloid cell","T cell","mast cell","adipocyte of breast" ))

plt1<-DotPlot(dat,assay = "SCT",features, dot.scale = 10,cluster.idents = FALSE)+ RotatedAxis()+ scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3))

ggsave(plt1,file="hbca_snRNA_markers.pdf",height=10,width=30,limitsize = FALSE)

#apriori chomvar markers
Idents(dat)<-dat$HBCA_predicted.id
tf_motifs<-FindAllMarkers(dat,
assay="chromvar",    
only.pos = TRUE,
test.use = 'LR',
latent.vars = 'nCount_peaks')

#get top 5 per cluster
DefaultAssay(dat)<-"peaks"
top_5<-as.data.frame(tf_motifs %>% arrange(-desc(p_val_adj)) %>% group_by(cluster) %>% slice_head(n = 10) %>% ungroup)
top_5$tf_gene<-unlist(llply(top_5$gene,function(x)ConvertMotifID(object=dat,id=x)))
top_5<-top_5[which(!duplicated(top_5$gene)),]
top_5_gene<-split(top_5$gene,f=top_5$cluster)
top_5_genename<-split(top_5$tf_gene,f=top_5$cluster)
top_5$cluster<-factor(top_5$cluster,levels=c("luminal epithelial cell of mammary gland",
  "basal cell","fibroblast","endothelial cell of lymphatic vessel",
  "endothelial cell of vascular tree","pericyte","myeloid cell","T cell","mast cell","adipocyte of breast" ))
top_5<-top_5[order(top_5$cluster),]
top_5<-top_5[top_5$tf_gene %in% row.names(GetAssayData(dat,"SCT")),]

pltrna<-DotPlot(dat,assay = "SCT",top_5$tf_gene,dot.scale = 20,cluster.idents = FALSE)+ RotatedAxis() + scale_color_gradient2(low="#7f3b08",mid="white",high="#542788")+scale_x_discrete(labels=top_5$tf_gene)
plt1<-DotPlot(dat,assay = "chromvar",top_5$gene,dot.scale = 20,cluster.idents = FALSE)+ RotatedAxis() + scale_color_gradient2(low="#7f3b08",mid="white",high="#542788")+scale_x_discrete(labels=top_5$tf_gene)+scale_y_discrete()
plt1$data$pct.exp<-pltrna$data$pct.exp
ggsave(plt1,file="apriori_chromvar_tfs.pdf",height=20,width=30,limitsize = FALSE)

DefaultAssay(dat)<-"peaks"
ga_markers<-FindAllMarkers(dat,
assay="GeneActivity",    
only.pos = TRUE,
test.use = 'LR',
latent.vars = 'nCount_peaks')

for(i in unlist(ga_markers$gene)){
plt<-CoveragePlot(
  object = dat,
  region=i,
  extend.upstream=5000,
  expression.assay="SCT",
  extend.downstream=5000,
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = TRUE)
ggsave(plt,file="feature.pdf")
}


