library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)
library(BiocParallel)
library(universalmotif)
library(GenomicRanges)
library(patchwork)
library(optparse)
register(SerialParam()) #using single core mode


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
outdir<-opt$plot_output_directory

DefaultAssay(dat)<-"ATAC"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species =9606, all_versions = FALSE))

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(dat[["ATAC"]]))) %in% main.chroms)
dat[["ATAC"]] <- subset(dat[["ATAC"]], features = rownames(dat[["ATAC"]][keep.peaks]))

# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
peaks<-granges(dat[["ATAC"]])

motif.matrix.hg38 <- CreateMotifMatrix(features = peaks, 
  pwm = pfm, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  use.counts = FALSE)

motif.hg38 <- CreateMotifObject(data = motif.matrix.hg38, 
  pwm = pfm)

dat <- SetAssayData(object = dat, 
  assay = 'ATAC', 
  slot = 'motifs', 
  new.data = motif.hg38)

dat <- RegionStats(object = dat, 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="ATAC")

dat <- RunChromVAR( object = dat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay="ATAC")

saveRDS(dat,file="merged.chromvar.SeuratObject.rds")

#Many transcription factors share the same motifs. To account for this, we are also going to perform chromvar across TF families.
#download cluster root motifs
system("wget --no-check-certificate https://jaspar2020.genereg.net/static/clustering/2020/vertebrates/CORE/interactive_trees/JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.tf") #use JASPAR2020 motif clusters
system("wget --no-check-certificate https://jaspar2020.genereg.net/static/clustering/2020/vertebrates/CORE/interactive_trees/JASPAR_2020_matrix_clustering_vertebrates_central_motifs_IDs.tab")
tf<-read_transfac("JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.tf") #read in transfac format
tf_cluster_names<-read.table("JASPAR_2020_matrix_clustering_vertebrates_central_motifs_IDs.tab",sep="\t",header=F)
#set up PWMatrix-List
pfm<-lapply(tf,function(x) convert_motifs(x,class="TFBSTools-PWMatrix"))
names(pfm)<-lapply(pfm,function(x) x@name)
pfm<-do.call(PWMatrixList,pfm)
names(pfm)<-tf_cluster_names[match(names(pfm),tf_cluster_names$V1),]$V3 #use readable names from jaspar

#Run regular chromvar
# Scan the DNA sequence of each peak for the presence of each motif, using orgo_atac for all objects (shared peaks)
DefaultAssay(dat)<-"ATAC"
motif.matrix <- CreateMotifMatrix(features = granges(dat), pwm = pfm, genome = 'hg38', use.counts = FALSE)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)

# Add the Motif object to the assays and run ChromVar
dat_chrom<- SetAssayData(object = dat, assay = 'ATAC', slot = 'motifs', new.data = motif) #write to dat_chrom so full motif list is not overwritten
dat_chrom<- RegionStats(object = dat_chrom, genome = BSgenome.Hsapiens.UCSC.hg38)
dat_chrom<- RunChromVAR(object = dat_chrom,genome = BSgenome.Hsapiens.UCSC.hg38,new.assay.name="jaspar_tffamily")

dat[["jaspar_tffamily"]]<-dat_chrom@assays$jaspar_tffamily
saveRDS(dat,file="merged.chromvar.SeuratObject.rds")