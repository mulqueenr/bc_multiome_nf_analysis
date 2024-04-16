#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/harmony.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(optparse)


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
#opt$object_input="merged.geneactivity.SeuratObject.rds"
dat=readRDS(opt$object_input)
DefaultAssay(dat)<-"peaks"
dat <- RunTFIDF(dat)
dat <- FindTopFeatures(dat, min.cutoff = 50)
dat <- RunSVD(dat, n = 30, reduction.name = 'lsi', reduction.key = 'LSI_')
dat <- RunUMAP(dat, reduction = 'lsi', dims = 2:30,reduction.name="unintegrated_atac_umap")
p1 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5,reduction="unintegrated_atac_umap",raster = FALSE) + ggplot2::ggtitle("Unintegrated ATAC")

dat <- RunHarmony(
  object = dat,
  group.by.vars = 'sample',
  reduction.use = 'lsi',
  assay.use = 'peaks',
  reduction.save= 'harmony_atac',
  project.dim = FALSE, ncores=10
)

# re-compute the UMAP using corrected LSI embeddings
dat <- RunUMAP(dat, dims = 2:30, reduction = 'harmony_atac',reduction.name="integrated_atac_umap")
p2 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5,reduction="integrated_atac_umap",raster = FALSE) + ggplot2::ggtitle("Harmony ATAC")

DefaultAssay(dat)<-"SCT"
dat <- RunPCA(dat)
dat <- RunUMAP(dat, reduction = 'pca', dims = 1:50,reduction.name="unintegrated_rna_umap")
p3 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5,reduction="unintegrated_rna_umap",raster = FALSE) + ggplot2::ggtitle("Unintegrated RNA")

dat <- RunHarmony(
  object = dat,
  group.by.vars = 'sample',
  reduction.use = 'pca',
  assay.use = 'RNA',
  reduction.save= 'harmony_rna',
  project.dim = FALSE, ncores=10
)

# re-compute the UMAP using corrected LSI embeddings
dat <- RunUMAP(dat, dims = 1:30, reduction = 'harmony_rna',reduction.name="integrated_rna_umap")
p4 <- DimPlot(dat, group.by = 'HBCA_predicted.id', pt.size = 0.5,reduction = "integrated_rna_umap",raster = FALSE) + ggplot2::ggtitle("Harmony RNA")


plt<-(p1 + p2)/(p3+p4)
ggsave(plt,file="harmony_integrations.pdf",height=10,width=20)

#joint with both
dat <- FindMultiModalNeighbors(
object = dat,
reduction.list = list("harmony_rna", "harmony_atac"),
dims.list = list(1:50, 2:30),
modality.weight.name = "RNA.weight",
verbose = TRUE
)

dat <- RunUMAP(
object = dat ,
nn.name = "weighted.nn",
reduction.name = "wnn.umap",
assay = "RNA",
verbose = TRUE
)
plt2<-DimPlot(dat ,reduction = "wnn.umap", group.by = c('HBCA_predicted.id','sample','Diagnosis','Mol_Diagnosis'),pt.size=0.5,raster=FALSE)
ggsave(plt2,file="harmony_integrations_multimodal.pdf",height=20,width=20)


#cluster with/without epithelial cells