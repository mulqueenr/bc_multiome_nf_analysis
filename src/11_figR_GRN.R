#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
# Load Packages
library(Signac)
library(Seurat)
library(FigR)
library(GenomicRanges)
library(harmony)
library(ggplot2)
library(patchwork)
library(optparse)
library(doParallel)
registerDoParallel(cores = 4)

setwd("/home/groups/MohammedLab/bc_multiome/seurat_object")

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
dat<-subset(dat,cells=row.names(dat@meta.data)[is.na(dat@meta.data$merged_assay_clones) | dat@meta.data$merged_assay_clones != "contamination"])
table(dat$sample,dat$assigned_celltype)

cell_metadata <- dat@meta.data

#prepare atac
atac_matrix <- GetAssayData(dat, assay = "ATAC", layer = "counts")
atac_peaks <- granges(dat[["ATAC"]])

ATAC.se <- SummarizedExperiment(
  assays = list(counts = atac_matrix),
  rowRanges = atac_peaks,
  colData = cell_metadata)

#prepare RNA
dat <- JoinLayers(dat, assay = "RNA")
rna_matrix <- GetAssayData(dat, assay = "RNA", layer = "data")
rna_matrix <- rna_matrix[Matrix::rowSums(rna_matrix) > 0, ]

#use clustering for knn graph #CHECK THIS PART
knn_graph <- dat@graphs$allcells.wsnn
knn_matrix <- as(knn_graph, "dgCMatrix")

cis_correlations <- runGenePeakcorr(
  ATAC.se = ATAC.se,
  RNAmat = rna_matrix,
  genome = "hg38",
  nCores = 4,
  p.cut = NULL,
  normalizeATACmat = TRUE)

#filter
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

#plot genes by peak associations
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                         cutoff = 10, # No. sig peaks needed to be called a DORC
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=2)

dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 4)

dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 4)
