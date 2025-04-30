#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif
#cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

library(Seurat)
library(Signac)
library(ggplot2)
library(optparse)
library(dplyr)
library(patchwork)
set.seed(123)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="7_merged.scsubtype.SeuratObject.rds", 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

clin_col=c(
"DCIS DCIS"="#cccccb", 
"ILC ER+/PR+/HER2-"="#f6bea1", 
"ILC ER+/PR-/HER2-"="#b9db98", 
"IDC ER+/PR-/HER2+"="#f37872", 
"IDC ER+/PR+/HER2-"="#8d86c0", 
"IDC ER+/PR-/HER2-"="#7fd0df", 
"NAT NA"="#c2d9ea")


#write out DA sites for tornado plots
Idents(dat)<-dat$Diag_MolDiag
dat_cancer<-subset(dat,cells=colnames(dat)[dat$assigned_celltype %in% c("cancer")])

Idents(dat_cancer)<-factor(dat_cancer$Diag_MolDiag,
    levels=c("DCIS DCIS","ILC ER+/PR+/HER2-", "ILC ER+/PR-/HER2-", "IDC ER+/PR-/HER2+", "IDC ER+/PR+/HER2-", "IDC ER+/PR-/HER2-", "NAT NA"))

da_peaks<-FindAllMarkers(dat_cancer,assay="ATAC")
write.table(da_peaks,file="diagnosis_da_peaks.csv",col.names=T,row.names=T,sep=",")

#setup outnames to be more nice
outnames=gsub(x=names(clin_col),pattern="[+]",replacement="plus")
outnames=gsub(x=outnames,pattern="[-]",replacement="neg")
outnames=gsub(x=outnames,pattern="[/]",replacement="_")
outnames=gsub(x=outnames,pattern="[ ]",replacement="_")
outnames=setNames(nm=names(clin_col),outnames)

#making it smaller by removing other assays
DefaultAssay(dat_cancer)<-"ATAC"

#subsetting object to just chromatin assay
dat_cancer<-subset(dat,cells=colnames(dat)[dat$assigned_celltype %in% c("cancer")])
atac_assay<-dat_cancer[["ATAC"]]
meta<-dat_cancer@meta.data
dat_cancer <- CreateSeuratObject(
  counts = atac_assay,
  assay = 'ATAC',
  project = 'ATAC',
  meta.data = meta
)
Idents(dat_cancer)<-dat_cancer$Diag_MolDiag
dat_cancer<-subset(x = dat_cancer, downsample = 500)
table(Idents(dat_cancer))
#all DA peaks
tornado_plot<-function(obj=dat_cancer,da_peak_set=da_peaks,i=da_peaks$cluster,peak_count=100){
    print(i)
    top_peaks <- da_peak_set %>% filter(cluster==i) %>% arrange(p_val_adj) %>% slice_head(n=peak_count)
    print(head(top_peaks))

    obj<-RegionMatrix(obj,key="DA_mat",
    regions=StringToGRanges(top_peaks$gene),
    upstream=5000,downstream=5000,
    assay="ATAC")

    plt<-RegionHeatmap(obj,key="DA_mat",
    upstream=5000,downstream=5000,
    order=TRUE, window=(10000)/100,
    assay="ATAC", idents=levels(Idents(obj)),
    nrow=length(unique(da_peak_set$cluster)))+ 
    scale_fill_gradient2(low="black",mid=unname(clin_col[i]),midpoint=0.08,high="white",na.value="black",breaks=seq(0.02,0.1,0.01))+
    ggtitle(i)
    #cols=setNames(nm=c("ATAC"),c(unname(clin_col[i]))),

    print("Returning plot...")
    return(plt)
}


plt<-tornado_plot(obj=dat_cancer,da_peak_set=da_peaks,i=x)
ggsave(plt,file="tornado.all_da_peaks.diag.pdf",width=20,height=20)


#all peaks
plt_list<-lapply(unique(da_peaks$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1,guides='collect')
ggsave(plt,file="tornado.all_da_peaks.diag.pdf",width=20,height=20)

#plt<-tornado_plot(obj=dat_cancer,da_peak_set=da_peaks,i=x)
#ggsave(plt,file="tornado.all_da_peaks.diag.pdf",width=20,height=20)

#da peaks that overlap with ESR1 motif only
motif_name="ESR1"
motif<-names(dat@assays$ATAC@motifs@motif.names[which(dat@assays$ATAC@motifs@motif.names==motif_name)])
da_peaks_motif_filt<-da_peaks[da_peaks$gene %in% names(which(dat@assays$ATAC@motifs@data[,motif])),]

plt_list<-lapply(unique(da_peaks_motif_filt$cluster),function(x) {
    tornado_plot(obj=dat_cancer,da_peak_set=da_peaks_motif_filt,i=x)
    })
plt<-wrap_plots(plt_list,nrow=1)
ggsave(plt,file="tornado.ESR1_overlap_da_peaks.diag.pdf",width=20,height=20)
