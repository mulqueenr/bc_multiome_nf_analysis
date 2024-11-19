#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
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


dat=readRDS(opt$object_input)

pseudobulk_dat <- AggregateExpression(dat, assays = "RNA", return.seurat = T, group.by = c("sample"))
saveRDS(pseudobulk_dat,file="sample.pseudobulk_pam50.rds")
