#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects")
dat_in<-readRDS("merged.public_transfer.SeuratObject.rds")

# Get motif data
data(motifs)
#convert assays to seurat v3 assays
options(Seurat.object.assay.version = "v3")
dat <- CreateSeuratObject(counts = dat_in)

dat[["RNA"]] <- as(object = dat[["RNA"]], Class = "Assay")
dat[["peaks"]] <- as(object = dat[["peaks"]], Class = "ChromatinAssay")

# Select variable features
dat <- Seurat::FindVariableFeatures(dat, assay='RNA')

# Initiate GRN object and select candidate regions
dat <- initiate_grn(dat) #####this is broken, have to convert to v3

# Scan candidate regions for TF binding motifs
seurat_object <- find_motifs(
    seurat_object,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# Infer gene regulatory network
seurat_object <- infer_grn(seurat_object)

# Print inferred coefficients
coef(seurat_object)

# Find gene and regulatory modules 
test_srt <- find_modules(test_srt)

# Print modules
NetworkModules(test_srt)