library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
set.seed(1234)
library(stringr)
library(CaSpER) 
args = commandArgs(trailingOnly=TRUE)

dat=readRDS(arg[1])
outname<-strsplit(arg[2],"[.]")[1]
outdir=arg[3]
ref_dir=arg[4]
cellranger_in=arg[5]
outdir_sample=paste0(outdir,"/",outname,"/casper/")
system(paste0("mkdir -p ",outdir_sample))

dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname

bam_location<-paste0(dir_in,"/gex_possorted_bam.bam")
BAFExtract_location<-paste0(ref_dir,"/BAFExtract/bin/BAFExtract")
hg38_list_location<-paste0(ref_dir,"/BAFExtract/hg38.list") #downloaded from https://github.com/akdess/BAFExtract
hg38_folder_location<-paste0(ref_dir,"/BAFExtract/hg38/")
baf_sample_directory<-paste0(dir_in,"/casper")

DefaultAssay(dat)<-"RNA"
dat$cnv_ref<-"FALSE"
dat@meta.data[!(dat$EMBO_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
control<-names(dat$cnv_ref == "TRUE") 
log.ge <- as.matrix(dat@assays$RNA@data)
genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = F)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
log.ge <- log2(log.ge +1)

system(paste0("samtools view ",bam_location," | ",BAFExtract_location," -generate_compressed_pileup_per_SAM stdin ",hg38_list_location," ",baf_sample_directory," 30 0 && wait;")) #generate BAF calls
#example of actual call: samtools view /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs/gex_possorted_bam.bam| /home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper 30 0 &

system(paste0(BAFExtract_location," -get_SNVs_per_pileup ",hg38_list_location," ",baf_sample_directory," ",hg38_folder_location," 1 1 0.1 ",baf_sample_directory,"/test.snp")) #generage snv files from BAF
#example of actual call: /home/groups/CEDAR/mulqueen/src/BAFExtract/bin/BAFExtract -get_SNVs_per_pileup /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38.list /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper /home/groups/CEDAR/mulqueen/src/BAFExtract/hg38/ 1 1 0.1 /home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_1/outs_casper/test.snp
  
loh <- readBAFExtractOutput ( path=baf_sample_directory, sequencing.type="bulk") 
names(loh) <- gsub(".snp", "", names(loh))
load(paste0(hg38_folder_location,"/maf.rda")) ## from https://github.com/akdess/CaSpER/blob/master/data/maf.rda
loh<- list()
loh[[1]] <- maf
names(loh) <- sample_name
loh.name.mapping <- data.frame (loh.name= sample_name , sample.name=colnames(log.ge))

#analysis demonstration: https://rpubs.com/akdes/673120
object <- CreateCasperObject(raw.data=log.ge,
loh.name.mapping=loh.name.mapping, 
sequencing.type="single-cell", 
cnv.scale=3, 
loh.scale=3, 
expr.cutoff=0.1, 
filter="median", 
matrix.type="normalized",
annotation=annotation, 
method="iterative", 
loh=loh, 
control.sample.ids=control, 
cytoband=cytoband)

saveRDS(object,paste0(outdir_sample,"/casper.initialobj.rds"))

## runCaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")
saveRDS(final.objects,paste0(outdir_sample,"/casper.finalobj.rds"))

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents(final.objects, thr=0.75)
final.obj <- final.objects[[9]]

saveRDS(final.obj,paste0(outdir_sample,"/casper.finalobj.rds"))
saveRDS(finalChrMat,paste0(outdir_sample,"/casper.finalchrmat.rds"))

#Segmentations
gamma <- 6
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary(final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]

#summrize segmentation across genes
all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann) #just need to fix genes.ann
saveRDS(rna.matrix, paste0(outdir_sample,"/casper.finalgenemat.rds"))
