```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(Seurat)
library(Signac,lib="/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3")
library(ggplot2)
library(patchwork)
library(dplyr)
library(optparse)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")


option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="merged.clone_annot.passqc.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

DefaultAssay(dat)<-"peaks"
#filter and create final fragments file
#tmpf <- tempfile(fileext = ".gz")

#run 5 cores to correct fragment files for cell renaming
#frag_update<-mclapply(2:length(Fragments(dat)), function(x) {
#    frag_in<-Fragments(dat)[[x]]
#    cells_filt<-names(frag_in@cells)[names(frag_in@cells) %in% colnames(dat)]
#    if(length(cells_filt)>0){
#    outpath=paste(dirname(frag_in@path),"passqc.atac_fragments.tsv.gz",sep="/")
#    FilterCells(frag_in@path,
#    cells=unlist(lapply(strsplit(cells_filt,"_"),"[",3)),
#    outfile = outpath, buffer_length = 256L, verbose = TRUE)
#    frag_update <- CreateFragmentObject(path = outpath,
#    cells = setNames(unlist(lapply(strsplit(cells_filt,"_"),"[",3)),nm=cells_filt),
#    validate.fragments = TRUE)
#    return(frag_update)
#    }
#},mc.cores=5)

#Fragments(dat) <- NULL
#Fragments(dat) <- frag_update
#saveRDS(dat,file=opt$object_input)
#cell_idx<-unlist(lapply(1:length(frag_update), function(x) names(frag_update[[x]]@cells)))
#dat<-subset(dat,cells=cell_idx)


dat$Diag_MolDiag<-paste(dat$Diagnosis,dat$Mol_Diagnosis)
dat<-ScaleData(dat,assay="RNA",features=row.names(dat@assays$RNA@data),vars.to.regress="nCount_RNA")

dat_epi<-subset(dat,cell=row.names(dat@meta.data)[dat$assigned_celltype %in% c("cancer_luminal_epithelial")])
dat_epi<-ScaleData(dat_epi,features=c("ERBB2","CDH1"))
Idents(dat_epi)<-dat_epi$Diagnosis
feat="ERBB2"
cov_plot <- CoveragePlot(object = dat_epi,
                        region=feat,
                        features=feat,
                        expression.slot="scale.data",
                        assay.scale="common",
                        extend.upstream=50000,
                        extend.downstream=50000,
                        split.by="Diag_MolDiag",
                        ncol=1)
ggsave(cov_plot,file="coverage.ERBB2.pdf")


feat="CDH1"
cov_plot <- CoveragePlot(object = dat_epi,
                        region=feat,
                        features=feat,
                        expression.slot="scale.data",
                        assay.scale="common",
                        extend.upstream=50000,
                        extend.downstream=50000,
                        split.by="Diag_MolDiag",
                        ncol=1)
ggsave(cov_plot,file="coverage.CDH1.pdf")
