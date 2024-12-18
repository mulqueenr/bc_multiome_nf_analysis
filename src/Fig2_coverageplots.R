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
tmpf <- tempfile(fileext = ".gz")

#run 5 cores
frag_update<-mclapply(2:length(Fragments(dat)), function(x) {
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

#rename cells to sample id
#TO HERE
Fragments(dat) <- NULL
Fragments(dat) <- frag_update

saveRDS(dat,file=opt$object_input)


dat$ploidy<-"aneuploid"

cell_idx<-unlist(lapply(1:length(frag_update), function(x) names(frag_update[[x]]@cells)))
DefaultAssay(dat)<-"peaks"
dat<-subset(dat,cells=cell_idx)
cov_plot <- CoveragePlot(object = dat,
                        region="TP63",
                        features="TP63",
                        assay.scale="common",
                        extend.upstream=5000,
                        extend.downstream=5000,
                        ncol=1)
ggsave(cov_plot,file="coverage.pdf",height=10)


ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = TRUE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  all.groups <- as.character(x = unique(x = groups))
  if (any(is.na(x = groups))) {
    all.groups <- c(all.groups, NA)
  }
  ngroup <- length(x = all.groups)
  npos <- ncol(x = mat)

  group <- unlist(
    x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
  )
  position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
  count <- vector(mode = "numeric", length = npos * ngroup)

  for (i in seq_along(along.with = all.groups)) {
    grp <- all.groups[[i]]
    if (is.na(x = grp)) {
      pos.cells <- names(x = groups)[is.na(x = groups)]
    } else {
      pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    }
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    count[((i - 1) * npos + 1):((i * npos))] <- totals
  }

  # construct dataframe
  coverages <- data.frame(
    "group" = group, "position" = position, "count" = count,
    stringsAsFactors = FALSE
  )

  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}
