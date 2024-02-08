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
Subsetting to just epithelial cells for now. 
https://buenrostrolab.github.io/FigR/
```R
library(FigR)
library(Seurat)
library(Signac)

dat<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")

dat<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))

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
                        nCores = 1, 
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
# Plot genomic ranges with Signac

```R
#module load singularity
#cd /home/groups/CEDAR/mulqueen/bc_multiome
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome multiome_bc.sif

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)

dat<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")

dat_downsamp<- dat[, sample(colnames(dat), size =10000, replace=F)]
plt1<-DimPlot(dat_downsamp,reduction="umap",group.by="Diagnosis")
DefaultAssay(dat)<-"RNA"
plt3<-DimPlot(dat_downsamp,reduction="umap",group.by="HBCA_predicted.id")

ggsave(plt1/plt3,file="umap.pdf",width=10,height=20)


library(plyr)
#snRNA markers
hbca_snmarkers=list()
hbca_snmarkers[["lumhr"]]=c("ANKRD30A","AFF3","ERBB4")
hbca_snmarkers[["lumsec"]]=c("COBL","ELF5","KIT")
hbca_snmarkers[["basal"]]=c("CARMN","ACTA2","KLHL29")
hbca_snmarkers[["fibro"]]=c("LAMA2","ANK2","SLIT2")
hbca_snmarkers[["lymphatic"]]=c("KLHL4","RHOJ","MMRN1")
hbca_snmarkers[["vascular"]]=c("MECOM","VWF","LDB2")
hbca_snmarkers[["perivasc"]]=c("RGS6","COL25A1","ADGRL3")
hbca_snmarkers[["myeloid"]]=c("F13A1","CD163","RAB31")
hbca_snmarkers[["tcells"]]=c("PTPRC","THEMIS","CD247")
hbca_snmarkers[["mast"]]=c("SLC24A3","HPGD","HDC")
hbca_snmarkers[["adipo"]]=c("PCDH9","CLSTN2","TRHDE")
features<-llply(hbca_snmarkers, unlist)

#Idents(dat)<-factor(dat$EMBO_predicted.id,levels=rev(c("epithelial","cycling.epithelial","CAFs","Endothelial","Pericytes","Myeloid","TAMs","TAMs_2","Plasma.cells","B.cells","T.cells","NA")))
Idents(dat)<-factor(dat$HBCA_predicted.id,levels=c("luminal epithelial cell of mammary gland",
  "basal cell","fibroblast","endothelial cell of lymphatic vessel",
  "endothelial cell of vascular tree","pericyte","myeloid cell","T cell","mast cell","adipocyte of breast" ))

plt1<-DotPlot(dat, features = features,dot.scale = 20,cluster.idents = FALSE) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +RotatedAxis()+
  scale_color_gradient2(low="#313695",mid="#ffffbf",high="#a50026",limits=c(-1,3)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

ggsave(plt1,file="hbca_markers.pdf",height=10,width=30,limitsize = FALSE)

outdir="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/plots"
dat_backup<-dat



dat<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
Idents(dat)<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
DefaultAssay(dat)<-"peaks"
cov_plot <- CoveragePlot(
  object = dat,
  region = "ESR1",extend.upstream = 2000, extend.downstream = 2000,
  annotation = TRUE,
  peaks = TRUE
)
ggsave(cov_plot,file="coverage.pdf")
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


```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome
module load singularity
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif


library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)

dat<-readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.public_transfer.SeuratObject.rds")
```

# Compare topics

## DecoupleR
