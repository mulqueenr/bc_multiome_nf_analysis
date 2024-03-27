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


plt<-(p1 + p2)/(p3+p4)+ plot_layout(guides = "collect")
ggsave(plt,file="harmony_integrations.pdf",height=10,width=20)


#cluster with/without epithelial cells