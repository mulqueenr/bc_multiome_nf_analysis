#https://www.sciencedirect.com/science/article/pii/S2666979X22001082
#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
# Load Packages
library(FigR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects")
dat_in<-readRDS("merged.public_transfer.SeuratObject.rds")



ATAC.SE<-SummarizedExperiment(assays=dat_in[["peaks"]]@counts,
                          rowRanges=dat_in[["peaks"]]@ranges,
                          metadata=dat_in@meta.data,
                          checkDimnames=TRUE)
assayNames(ATAC.SE)<-"counts"
rnaMat<-as(object = JoinLayers(dat_in[["RNA"]])$counts, Class = "Matrix")



# Run using multiple cores if parallel support
cisCor <- runGenePeakcorr(ATAC.se = ATAC.SE,
                           RNAmat = rnaMat,
                           genome = "hg38", # Also supports mm10 and hg38
                           nCores = 10, 
                           p.cut=NULL)

# Filter peak-gene correlations by p-value                    
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)

# Determine DORC genes
dorcGenes <- cisCor.filt %>% dorcJPlot(cutoff=7, returnGeneList = TRUE)

# Get DORC scores
dorcMat <- getDORCScores(ATAC.SE,dorcTab=cisCor.filt,geneList=dorcGenes,nCores=10)

# Smooth DORC scores (using cell KNNs)
dorcMat.smooth <- smoothScoresNN(NNmat=cellKNN.mat,mat=dorcMat,nCores=10)

rnaMat.smooth <- smoothScoresNN(NNmat=cellKNN.mat,mat=rnaMat,nCores=10)

# Run FigR
fig.d <- runFigRGRN(ATAC.se=ATAC.SE,
                    rnaMat=rnaMat.smooth, # Smoothed RNA matrix using paired cell kNNs
                    dorcMat=dorcMat.smooth,
                    dorcTab=cisCor.filt,
                    genome="hg38",
                    dorcGenes=dorcGenes,
                    nCores=4)