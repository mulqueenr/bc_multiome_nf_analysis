#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(optparse)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="4_merged.chromvar.SeuratObject.rds", 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

gene_activity<-GeneActivity(dat,process_n=10000)

dat[["GeneActivity"]]<-CreateAssayObject(counts=gene_activity)
dat<- NormalizeData(
  object = dat,
  assay = "GeneActivity",
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_GeneActivity)
)
saveRDS(dat,file="5_merged.geneactivity.SeuratObject.rds")


