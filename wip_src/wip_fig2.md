```bash
#module load singularity
cd /home/groups/CEDAR/mulqueen/bc_multiome
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome multiome_bc.sif
R
```

# GRN inferences
## Pando
https://github.com/quadbio/Pando
```R
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
registerDoParallel(20)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE) #updating signac to make this work with seurat5
library(Signac)

# Get motif data
data('motifs')
data('motif2tf')


dat<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")

# Select variable features
dat <- Seurat::FindVariableFeatures(dat, assay='RNA')

dat[["peaks"]]<-as(object = dat[["peaks"]], Class = "Assay")
# Initiate GRN object and select candidate regions
dat <- initiate_grn(dat,
	rna_assay = "RNA",
    peak_assay = "ATAC")

# Scan candidate regions for TF binding motifs
dat <- find_motifs(
    dat,
    pfm = motifs,
    motif_tfs=motif2tf,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# Infer gene regulatory network
dat <- infer_grn(dat,peak_to_gene_method='GREAT',parallel=T)

# Print inferred coefficients
coef(dat)

# Find gene and regulatory modules 
test_srt <- find_modules(dat)

# Print modules
NetworkModules(test_srt)
```


## FigR
https://buenrostrolab.github.io/FigR/
```R
library(FigR)
library(Seurat)
library(Signac)

dat<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")


# Create a summarized experiment
ATAC.SE <- SummarizedExperiment(
rowRanges = dat@assays$peaks@ranges,
colData = dat@meta.data,
assays = list(counts = dat@assays$peaks@counts)
)

# Run using multiple cores if parallel support
cisCor <- runGenePeakcorr(ATAC.se = ATAC.SE,
                        RNAmat = dat[["RNA"]]@counts,
                        genome = "hg38",
                        nCores = 5, 
                        p.cut=NULL)

        
# Filter peak-gene correlations by p-value                    
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)

# Determine DORC genes
dorcGenes <- cisCor.filt %>% dorcJPlot(cutoff=7, # Default
                                       returnGeneList = TRUE)

# Get DORC scores
dorcMat <- getDORCScores(ATAC.SE,dorcTab=cisCor.filt,geneList=dorcGenes,nCores=5)

# Smooth DORC scores (using cell KNNs)
dorcMat.smooth <- smoothScoresNN(NNmat=cellKNN.mat,mat=dorcMat,nCores=5)

# Run FigR
fig.d <- runFigRGRN(ATAC.se=,ATAC.SE,
                    rnaMat=rnaMat.smooth, # Smoothed RNA matrix using paired cell kNNs
                    dorcMat=dorcMat.smooth,
                    dorcTab=cisCor.filt,
                    genome="hg38",
                    dorcGenes=dorcGenes,
                    nCores=5)

```

# Cell type bias inferences

## scDC

https://sydneybiox.github.io/scDC/articles/vignette.html
```R
library(scDC)
library(Signac)
dat<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")

sample<-dat$sample
cellTypes<-dat$HBCA_predicted.id
cond<-unlist(lapply(strsplit(unique(sample),"_"),"[[",1))

res_scDC_noClust <- scDC_noClustering(cellTypes, sample, calCI = TRUE, 
                                     calCI_method = c("percentile", "BCa", "multinom"),
                                     nboot = 100)


pdf("scDC_test.pdf",width=10,height=10)
barplotCI(res_scDC_noClust, cond)
dev.off()


#res_GLM <- fitGLM(res_scDC_noClust, cond, pairwise = FALSE) to run this need to fix matrix https://community.rstudio.com/t/error-in-initializeptr-function-cholmod-factor-ldeta-not-provided-by-package-matrix/178694
```


# Compare topics

## DecoupleR
