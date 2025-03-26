#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(parallel)
library(EnsDb.Hsapiens.v86)
library(optparse)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="4_merged.chromvar.SeuratObject.rds", 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

#filter unneeded fragment files
frag_file_list<-unlist(lapply(1:length(Fragments(dat)), function(x) if(length(Fragments(dat)[[x]]@cells)>0) {return(x)}))
filtered_frag_files<-lapply(frag_file_list, function(x) Fragments(dat)[[x]])
Fragments(dat)<-NULL
Fragments(dat)<-filtered_frag_files

#filter and create final fragments file
tmpf <- tempfile(fileext = ".gz")

#run 5 cores to correct fragment files for cell renaming
frag_update<-mclapply(1:length(Fragments(dat)), function(x) {
    frag_in<-Fragments(dat)[[x]]
    cells_filt<-names(frag_in@cells)[names(frag_in@cells) %in% colnames(dat)]
    if(length(cells_filt)>0){
    outpath=paste(dirname(frag_in@path),"passqc.atac_fragments.tsv.gz",sep="/")
    FilterCells(frag_in@path,
    cells=unlist(lapply(strsplit(cells_filt,"_"),"[",3)),
    outfile = outpath, buffer_length = 256L, verbose = TRUE)
    frag_update <- CreateFragmentObject(path = outpath,
    cells = setNames(unlist(lapply(strsplit(cells_filt,"_"),"[",3)),nm=cells_filt),
    validate.fragments = TRUE)
    return(frag_update)
    }
},mc.cores=5)

Fragments(dat) <- NULL
Fragments(dat) <- frag_update

gene_activity<-GeneActivity(dat,process_n=10000)

dat[["GeneActivity"]]<-CreateAssayObject(counts=gene_activity)
dat<- NormalizeData(
  object = dat,
  assay = "GeneActivity",
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_GeneActivity)
)
saveRDS(dat,file="5_merged.geneactivity.SeuratObject.rds")


