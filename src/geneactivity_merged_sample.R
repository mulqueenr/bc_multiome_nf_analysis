library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
outdir<-opt$plot_output_directory

gene_activity<-GeneActivity(dat,process_n=10000)
saveRDS(gene_activity,file=paste0(outname,".GeneActivity.rds"))

dat[["GeneActivity"]]<-CreateAssayObject(counts=gene_activity)
dat<- NormalizeData(
  object = dat,
  assay = "GeneActivity",
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_GeneActivity)
)
saveRDS(dat,file=paste0(outname,".geneactivity.SeuratObject.rds"))


