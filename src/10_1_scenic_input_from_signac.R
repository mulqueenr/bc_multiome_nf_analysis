#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif
#cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

library(Seurat)
library(Signac)
library(ggplot2)
library(optparse)
library(dplyr)
library(Matrix)

set.seed(123)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)


### RUN FOR ALL CELLS ###
#write out fragment paths per cell
fragment_paths<-lapply(1:length(dat@assays$ATAC@fragments),function(x) cbind(dat@assays$ATAC@fragments[[x]]@path,unlist(names(dat@assays$ATAC@fragments[[x]]@cells)),unlist(unname(dat@assays$ATAC@fragments[[x]]@cells))))
fragment_paths<-as.data.frame(do.call("rbind",fragment_paths))
row.names(fragment_paths)<-fragment_paths$V2
write.table(fragment_paths,file="frag_paths.csv",sep=",",col.names=F,row.names=T)

#write out metadata
write.table(dat@meta.data,file="metadata.csv",col.names=T,row.names=T,sep=",")

#write out ATAC counts matrix
write.table(colnames(dat@assays$ATAC@counts),file="atac_counts.cells.csv",sep=",")
write.table(row.names(dat@assays$ATAC@counts),file="atac_counts.peaks.csv",sep=",")
writeMM(dat@assays$ATAC@counts,file="atac_counts.mtx")

#write out RNA counts
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])
dat[["RNA"]]<-as(object = dat[["RNA"]], Class = "Assay")
write.table(colnames(dat@assays$RNA@counts),file="rna_counts.cells.csv",sep=",")
write.table(row.names(dat@assays$RNA@counts),file="rna_counts.genes.csv",sep=",")
writeMM(dat@assays$RNA@counts,file="rna_counts.mtx")

### RUN FOR EPI AND CANCER ###
dat=readRDS(opt$object_input)

dat<-subset(dat,assigned_celltype %in% c("cancer","luminal_asp","luminal_hs","basal_myoepithelial"))
#write out fragment paths per cell
fragment_paths<-lapply(1:length(dat@assays$ATAC@fragments),function(x) cbind(dat@assays$ATAC@fragments[[x]]@path,unlist(names(dat@assays$ATAC@fragments[[x]]@cells)),unlist(unname(dat@assays$ATAC@fragments[[x]]@cells))))
fragment_paths<-as.data.frame(do.call("rbind",fragment_paths))
row.names(fragment_paths)<-fragment_paths$V2
write.table(fragment_paths,file="epi_frag_paths.csv",sep=",",col.names=F,row.names=T)

#write out metadata
write.table(dat@meta.data,file="epi_metadata.csv",col.names=T,row.names=T,sep=",")

#write out ATAC counts matrix
write.table(colnames(dat@assays$ATAC@counts),file="epi_atac_counts.cells.csv",sep=",")
write.table(row.names(dat@assays$ATAC@counts),file="epi_atac_counts.peaks.csv",sep=",")
writeMM(dat@assays$ATAC@counts,file="epi_atac_counts.mtx")

#write out RNA counts
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])
dat[["RNA"]]<-as(object = dat[["RNA"]], Class = "Assay")
write.table(colnames(dat@assays$RNA@counts),file="epi_rna_counts.cells.csv",sep=",")
write.table(row.names(dat@assays$RNA@counts),file="epi_rna_counts.genes.csv",sep=",")
writeMM(dat@assays$RNA@counts,file="epi_rna_counts.mtx")
