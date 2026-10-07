Download public data for comparisons
```bash
#also plot metabric and tcga samples
cd /home/groups/CEDAR/mulqueen/bc_multiome/ref
mkdir -p metabric_breast; cd metabric_breast
wget https://datahub.assets.cbioportal.org/brca_metabric.tar.gz 
tar -xvf brca_metabric.tar.gz
ls /home/groups/CEDAR/mulqueen/bc_multiome/ref/metabric_breast/brca_metabric


cd /home/groups/CEDAR/mulqueen/bc_multiome/ref
mkdir -p tcga_breast; cd tcga_breast
wget https://datahub.assets.cbioportal.org/brca_tcga_pan_can_atlas_2018.tar.gz
tar -xvf brca_tcga_pan_can_atlas_2018.tar.gz
ls /home/groups/CEDAR/mulqueen/bc_multiome/ref/tcga_breast/brca_tcga_pan_can_atlas_2018
```

```bash

sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell \
--bind /home/groups/CEDAR/mulqueen/bc_multiome \
$sif


cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

```


```R

library(Seurat)
library(Signac)
library(optparse)
library(rtracklayer)
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
library(dplyr)
library(parallel)
library(ggplot2)
library(ggrepel)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])
dat_sub<-subset(dat,cells=row.names(dat@meta.data[!is.na(dat@meta.data$merged_assay_clones),]))


#get CN bins that are variable
cn <- GetAssayData(dat_sub,assay="cnv",layer="data")
bin_var<-apply(cn,1,var)
bin_filter<-names(bin_var>0)
cn<-cn[bin_filter,]
bin_filter_granges<-StringToGRanges(bin_filter, sep = c("-", "-"))
bin_filter_granges$cn_win<-bin_filter

cellid_list<-colnames(cn)


#CN vs Motif and RNA
#motif
chromvar<-GetAssayData(dat_sub,assay="ATAC",layer="motifs")
chromvar_names<-chromvar@motif.names
chromvar_dat<-GetAssayData(dat_sub,assay="chromvar",layer="data")
row.names(chromvar_dat)<-unlist(chromvar_names[row.names(chromvar_dat)])
chromvar_dat<-chromvar_dat[,cellid_list]

#rna
rna_dat<-GetAssayData(dat_sub,layer="data",assay="RNA")
rna_dat<-rna_dat[,cellid_list]
rna_dat<-rna_dat[row.names(rna_dat) %in% row.names(chromvar_dat),]

chromvar_dat<-chromvar_dat[row.names(chromvar_dat) %in% row.names(rna_dat),]

#assign genes to cn bins
rna_ref <- import("/home/groups/CEDAR/mulqueen/bc_multiome/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
rna_ref <- rna_ref %>% as.data.frame() %>% filter(!duplicated(gene_name)) %>% filter(gene_name %in% row.names(chromvar_dat)) %>% GRanges()
rna_hits <- findOverlaps(rna_ref, bin_filter_granges, select="first")
rna_ref$cn_bin<-row.names(cn)[rna_hits]
rna_ref<-rna_ref[!is.na(rna_ref$cn_bin),]

cn_cor<-lapply(rna_ref$gene_name,function(gene){
  cn_bin<-rna_ref %>% as.data.frame() %>% filter(gene_name==gene) %>% select(cn_bin) %>% unlist()
  cn_tmp<-unlist(cn[cn_bin,cellid_list])
  cn_var<-var(cn_tmp,na.rm=T)

  rna_tmp<-unlist(rna_dat[gene,cellid_list])
  rna_var<-var(rna_tmp,na.rm=T)

  motif_tmp<-unlist(chromvar_dat[gene,cellid_list])
  motif_var<-var(motif_tmp,na.rm=T)

  cn_rna_cor<-cor.test(cn_tmp,rna_tmp,use="pairwise.complete",method="spearman")
  cn_motif_cor<-cor.test(cn_tmp,motif_tmp,use="pairwise.complete",method="spearman")
  motif_rna_cor<-cor.test(motif_tmp,rna_tmp,use="pairwise.complete",method="spearman")

  return(setNames(
      c(gene,cn_bin,cn_var,rna_var,motif_var,cn_rna_cor$p.value,cn_motif_cor$p.value,motif_rna_cor$p.value,cn_rna_cor$estimate,cn_motif_cor$estimate,motif_rna_cor$estimate),
      nm=c("gene","cn_bin","cn_var","rna_var","motif_var","cn_rna_cor_pval","cn_motif_cor_pval","motif_rna_cor_pval","cn_rna_cor_rho","cn_motif_cor_rho","motif_rna_cor_rho")))
})

cn_cor<-as.data.frame(do.call("rbind",cn_cor))
cn_cor<-cn_cor[complete.cases(cn_cor),] %>% mutate(across(c(cn_rna_cor_pval, cn_motif_cor_pval, motif_rna_cor_pval,
                                                            cn_rna_cor_rho, cn_motif_cor_rho,motif_rna_cor_rho,
                                                            cn_var,rna_var,motif_var), as.numeric))

cn_cor$cn_rna_cor_qval<-p.adjust(cn_cor$cn_rna_cor_pval,method="bonferroni")
cn_cor$cn_motif_cor_qval<-p.adjust(cn_cor$cn_motif_cor_pval,method="bonferroni")
cn_cor$motif_rna_cor_qval<-p.adjust(cn_cor$motif_rna_cor_pval,method="bonferroni")

cn_cor$label<-NA
cn_cor$min_qval<-apply(cn_cor %>% select(cn_motif_cor_qval,cn_rna_cor_qval),1,min)
label_genes<-unname(unlist(cn_cor %>% slice_min(min_qval,n=50) %>% select(gene)))
cn_cor$label<-ifelse(cn_cor$gene %in% label_genes, cn_cor$gene, NA) 

cn_cor$col<-"black"
rna_col<-unname(unlist(cn_cor %>% filter(cn_rna_cor_qval<0.05) %>% select(gene)))
motif_col<-unname(unlist(cn_cor %>% filter(cn_motif_cor_qval<0.05) %>% select(gene)))
both_col<-unname(unlist(intersect(rna_sig,motif_sig)))
cn_cor$col<-unlist(lapply(1:nrow(cn_cor), function(x){ifelse(cn_cor$gene[x] %in% rna_col,"red",cn_cor$col[x])}))
cn_cor$col<-unlist(lapply(1:nrow(cn_cor), function(x){ifelse(cn_cor$gene[x] %in% motif_col,"blue",cn_cor$col[x])}))
cn_cor$col<-unlist(lapply(1:nrow(cn_cor), function(x){ifelse(cn_cor$gene[x] %in% both_col,"purple",cn_cor$col[x])}))

plt1<-ggplot(cn_cor,aes(x=cn_rna_cor_rho,y=rna_var,size=abs(cn_var),label=label,color=col))+geom_point()+theme_minimal()+geom_text_repel(max.overlaps=Inf)+scale_color_identity()
plt2<-ggplot(cn_cor,aes(x=cn_rna_cor_rho,y=cn_motif_cor_rho,size=abs(motif_rna_cor_rho),label=label,color=col))+geom_point()+theme_minimal()+geom_text_repel(max.overlaps=Inf)+scale_color_identity()
plt3<-ggplot(cn_cor,aes(x=motif_var,y=cn_motif_cor_rho,size=abs(cn_var),label=label,color=col))+geom_point()+theme_minimal()+geom_text_repel(max.overlaps=Inf)+scale_color_identity()
ggsave(plt1+plot_spacer()+plt2+plt3+plot_layout(ncol=2,guide="collect"),file="test.pdf",width=30,height=30)
#run per bin, correlate bin to chromvar and rna data
```

```R
#and run metabric and tcga as validation
#metabric
  metabric_cnv<-read.table(sep="\t",file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/metabric_breast/brca_metabric/data_cna.txt",header=T)
  metabric_rna<-read.table(sep="\t",file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/metabric_breast/brca_metabric/data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt",header=T)
  metabric_meta<-read.table(sep="\t",file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/metabric_breast/brca_metabric/data_clinical_sample.txt",header=T)

#filter to match cancer types (ER+ IDC/ILC)
  metabric_meta<-metabric_meta %>% 
    filter(CANCER_TYPE_DETAILED %in% c("Breast Invasive Ductal Carcinoma","Breast Invasive Lobular Carcinoma","Breast Mixed Ductal and Lobular Carcinoma")) %>% 
    filter(ER_STATUS %in% c("Positive"))
  metabric_meta$sample<-gsub(metabric_meta$SAMPLE_ID,pattern="-",replacement=".")
  metabric_cnv<-metabric_cnv[,colnames(metabric_cnv) %in% c("Hugo_Symbol",unique(metabric_meta$sample))]
  metabric_rna<-metabric_rna[,colnames(metabric_rna) %in% c("Hugo_Symbol",unique(metabric_meta$sample))]
  metabric_rna<-metabric_rna[!duplicated(metabric_rna$Hugo_Symbol),]
  metabric_cnv<-metabric_cnv[!duplicated(metabric_cnv$Hugo_Symbol),]

  metabric_columns_to_keep<-intersect(colnames(metabric_cnv),colnames(metabric_rna))
  metabric_rows_to_keep<-intersect(metabric_cnv$Hugo_Symbol,metabric_rna$Hugo_Symbol)
  metabric_cnv<-metabric_cnv[metabric_cnv$Hugo_Symbol %in% metabric_rows_to_keep,colnames(metabric_cnv) %in% metabric_columns_to_keep]
  metabric_rna<-metabric_rna[metabric_rna$Hugo_Symbol %in% metabric_rows_to_keep,colnames(metabric_rna) %in% metabric_columns_to_keep]

dim(metabric_cnv)
dim(metabric_rna)


#tcga
  tcga_cnv<-read.table(sep="\t",file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/tcga_breast/brca_tcga_pan_can_atlas_2018/data_cna.txt",header=T)
  tcga_rna<-read.table(sep="\t",file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/tcga_breast/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt",header=T)
  tcga_meta<-read.table(sep="\t",file="/home/groups/CEDAR/mulqueen/bc_multiome/ref/tcga_breast/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt",header=T)

#filter to match cancer types (IDC/ILC)
#tcga doesnt report ER STATUS
  tcga_meta<-tcga_meta %>% 
    filter(CANCER_TYPE_DETAILED %in% c("Breast Invasive Ductal Carcinoma","Breast Invasive Lobular Carcinoma","Breast Mixed Ductal and Lobular Carcinoma")) 
  tcga_meta$sample<-gsub(tcga_meta$SAMPLE_ID,pattern="-",replacement=".")
  tcga_cnv<-tcga_cnv[,colnames(tcga_cnv) %in% c("Hugo_Symbol",unique(tcga_meta$sample))]
  tcga_rna<-tcga_rna[,colnames(tcga_rna) %in% c("Hugo_Symbol",unique(tcga_meta$sample))]
  tcga_rna<-tcga_rna[!duplicated(tcga_rna$Hugo_Symbol),]
  tcga_cnv<-tcga_cnv[!duplicated(tcga_cnv$Hugo_Symbol),]

  tcga_columns_to_keep<-intersect(colnames(tcga_cnv),colnames(tcga_rna))
  tcga_rows_to_keep<-intersect(tcga_cnv$Hugo_Symbol,tcga_rna$Hugo_Symbol)
  tcga_cnv<-tcga_cnv[tcga_cnv$Hugo_Symbol %in% tcga_rows_to_keep,colnames(tcga_cnv) %in% tcga_columns_to_keep]
  tcga_rna<-tcga_rna[tcga_rna$Hugo_Symbol %in% tcga_rows_to_keep,colnames(tcga_rna) %in% tcga_columns_to_keep]

#filter to shared genes across tcga and metabric
rows_to_keep<-row.names(chromvar_dat)
rows_to_keep<-rows_to_keep[rows_to_keep %in% tcga_rna$Hugo_Symbol]
rows_to_keep<-rows_to_keep[rows_to_keep %in% metabric_rna$Hugo_Symbol]



  tcga_cnv<-tcga_cnv[tcga_cnv$Hugo_Symbol %in% rows_to_keep,]
  tcga_rna<-tcga_rna[tcga_rna$Hugo_Symbol %in% rows_to_keep,]
  row.names(tcga_cnv)<-tcga_cnv$Hugo_Symbol; tcga_cnv<-tcga_cnv[,2:ncol(tcga_cnv)]
  row.names(tcga_rna)<-tcga_rna$Hugo_Symbol; tcga_rna<-tcga_rna[,2:ncol(tcga_rna)]

  metabric_cnv<-metabric_cnv[metabric_cnv$Hugo_Symbol %in% rows_to_keep,]
  metabric_rna<-metabric_rna[metabric_rna$Hugo_Symbol %in% rows_to_keep,]
  row.names(metabric_cnv)<-metabric_cnv$Hugo_Symbol; metabric_cnv<-metabric_cnv[,2:ncol(metabric_cnv)]
  row.names(metabric_rna)<-metabric_rna$Hugo_Symbol; metabric_rna<-metabric_rna[,2:ncol(metabric_rna)]
  

#note that for bulk public data CN is reported by gene level already
#also note CN is reported with 0 being diploid, which is a different center than our single cell data
#also also note RNA data is much different in bulk

cn_cor_tcga<-lapply(rows_to_keep,function(gene){
  if(gene %in% rna_ref$gene_name){
  cn_bin<-rna_ref %>% as.data.frame() %>% filter(gene_name==gene) %>% select(cn_bin) %>% unlist()
  cn_tmp<-unlist(tcga_cnv[gene,])
  cn_var<-var(cn_tmp,na.rm=T)

  rna_tmp<-unlist(tcga_rna[gene,colnames(tcga_rna) %in% names(cn_tmp)])
  rna_var<-var(rna_tmp,na.rm=T)

  cn_rna_cor<-cor.test(cn_tmp,rna_tmp,use="pairwise.complete",method="spearman")

  return(setNames(c(gene,cn_bin,cn_var,rna_var,cn_rna_cor$p.value,cn_rna_cor$estimate),
      nm=c("gene","cn_bin","cn_var","rna_var","cn_rna_cor_pval","cn_rna_cor_rho")))
}})
cn_cor_tcga<-as.data.frame(do.call("rbind",cn_cor_tcga))
colnames(cn_cor_tcga)<-paste0("tcga_",colnames(cn_cor_tcga))


cn_cor_metabric<-lapply(rows_to_keep,function(gene){
  if(gene %in% rna_ref$gene_name){
  cn_bin<-rna_ref %>% as.data.frame() %>% filter(gene_name==gene) %>% select(cn_bin) %>% unlist()
  cn_tmp<-unlist(metabric_cnv[gene,])
  cn_var<-var(cn_tmp,na.rm=T)

  rna_tmp<-log10(unlist(metabric_rna[gene,]))
  rna_var<-var(rna_tmp,na.rm=T)

  cn_rna_cor<-cor.test(cn_tmp,rna_tmp,use="pairwise.complete",method="spearman")

  return(setNames(
      c(gene, cn_bin, cn_var, rna_var, cn_rna_cor$p.value, cn_rna_cor$estimate),
      nm=c("gene","cn_bin","cn_var","rna_var","cn_rna_cor_pval","cn_rna_cor_rho")))
}})

cn_cor_metabric<-as.data.frame(do.call("rbind",cn_cor_metabric))
colnames(cn_cor_metabric)<-paste0("metabric_",colnames(cn_cor_metabric))

bulk_dat<-merge(cn_cor_tcga,cn_cor_metabric,by.x="tcga_gene",by.y="metabric_gene")
bulk_dat<-bulk_dat[complete.cases(bulk_dat),] %>% mutate(across(c(metabric_cn_rna_cor_pval, tcga_cn_rna_cor_pval,
                                                            metabric_cn_rna_cor_rho, tcga_cn_rna_cor_rho,
                                                            metabric_cn_var,metabric_rna_var,
                                                            tcga_cn_var,tcga_rna_var), as.numeric))

#labelling same genes as single cell data
bulk_dat$label<-NA
bulk_dat$label<-ifelse(bulk_dat$tcga_gene %in% label_genes, bulk_dat$tcga_gene, NA) 

#save into empty slot of previous plot
plt4<-ggplot(dat=bulk_dat,aes(x=metabric_cn_rna_cor_rho,y=tcga_cn_rna_cor_rho,label=label))+geom_point()+theme_minimal()+geom_text_repel()
ggsave(plt1+plt4+plt2+plt3+plot_layout(ncol=2,guide="collect"),file="test.pdf",width=30,height=30)

```

Additional plots of specific genes

```R
gene_list=c("GATA3", "GRHL1", "FOXA1", "SNAI2", "ZEB1","FOXO3","MEF2A","LMX1B","ESR1","SOX10","HNF4A","SREBF2")

meta_dat<-dat_sub@meta.data
meta_dat<-setNames(nm=row.names(dat_sub@meta.data),paste(dat_sub@meta.data$Diagnosis,dat_sub@meta.data$Mol_Diagnosis))

plt_list<-lapply(gene_list, function(gene) {
gene_rna<-rna_dat[gene,]
gene_motif<-chromvar_dat[gene,]
gene_cn<-cn[rna_ref[rna_ref$gene_name==gene,]$cn_bin,]
gene_dat<-data.frame(rna=gene_rna,motif=gene_motif[names(gene_rna)],cn=gene_cn[names(gene_rna)],meta=meta_dat[names(gene_rna)])

plt1<-ggplot(gene_dat,aes(x=paste(factor(cn),meta),y=rna,color=meta))+geom_jitter(alpha=0.2,size=0.2)+geom_violin(fill=NA)+geom_boxplot(fill=NA,outlier.shape = NA)+ggtitle(paste(gene,"RNA"))+facet_grid(.~cn,scales="free_x",space = "free_x")
plt2<-ggplot(gene_dat,aes(x=paste(factor(cn),meta),y=motif,color=meta))+geom_jitter(alpha=0.2,size=0.2)+geom_violin(fill=NA)+geom_boxplot(fill=NA,outlier.shape = NA)+ggtitle(paste(gene,"motif"))+facet_grid(.~cn,scales="free_x",space = "free_x")
return(plt1+plt2)
})

ggsave(wrap_plots(plt_list,ncol=1,guide="collect")*theme_minimal(),file="test2.pdf",height=length(gene_list)*3,width=30)


```



#run cistopic on cells with cnv clone assignment
cistopic_outdir=paste0(getwd(),"/cistopic_clones")
system(paste0("mkdir -p ", cistopic_outdir))

### RUN FOR ALL CELLS ASSIGNED CLONES###
#write out fragment paths per cell
fragment_paths<-lapply(1:length(dat@assays$ATAC@fragments),function(x) cbind(dat@assays$ATAC@fragments[[x]]@path,unlist(names(dat@assays$ATAC@fragments[[x]]@cells)),unlist(unname(dat@assays$ATAC@fragments[[x]]@cells))))
fragment_paths<-as.data.frame(do.call("rbind",fragment_paths))
row.names(fragment_paths)<-fragment_paths$V2
write.table(fragment_paths,file=paste0(cistopic_outdir,"/","frag_paths.csv"),sep=",",col.names=F,row.names=T)

#write out metadata
write.table(dat@meta.data,file=paste0(cistopic_outdir,"/","metadata.atac.csv"),col.names=T,row.names=T,sep=",")
write.table(dat@meta.data,file=paste0(cistopic_outdir,"/","metadata.rna.csv"),col.names=T,row.names=T,sep=",")

#write out ATAC counts matrix
write.table(colnames(dat@assays$ATAC@counts),file=paste0(cistopic_outdir,"/","atac_counts.cells.csv"),sep=",")
write.table(row.names(dat@assays$ATAC@counts),file=paste0(cistopic_outdir,"/","atac_counts.peaks.csv"),sep=",")
Matrix::writeMM(dat@assays$ATAC@counts,file=paste0(cistopic_outdir,"/","atac_counts.mtx"))

#write out RNA counts
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])
dat[["RNA"]]<-as(object = dat[["RNA"]], Class = "Assay")
write.table(colnames(dat@assays$RNA@counts),file=paste0(cistopic_outdir,"/","rna_counts.cells.csv"),sep=",")
write.table(row.names(dat@assays$RNA@counts),file=paste0(cistopic_outdir,"/","rna_counts.genes.csv"),sep=",")
Matrix::writeMM(dat@assays$RNA@counts,file=paste0(cistopic_outdir,"/","rna_counts.mtx"))
```

Load scenicplus singularity file for processing
```bash
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/scenicplus.sif
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/cistopic_clones
```


```python
import os
from pycisTopic.cistopic_class import *
from scipy.io import mmread
import pickle
import scanpy as sc
from pycisTopic.lda_models import run_cgs_models_mallet
import argparse


parser = argparse.ArgumentParser(
    description="Function to run single cistopic model.")

parser.add_argument("-f", "--frag_path", default = 'frag_paths.csv', help = "List of fragment locations, csv format")
parser.add_argument("-n", "--atac_counts", default = 'atac_counts.mtx', help = "Raw counts, mtx")
parser.add_argument("-c", "--atac_cells", default = 'atac_counts.cells.csv', help = "List of cells, csv")
parser.add_argument("-p", "--atac_peaks", default = 'atac_counts.peaks.csv', help = "List of peaks, csv")
parser.add_argument("-m", "--meta", default = 'metadata.atac.csv', help = "Metadata table, csv")
parser.add_argument("-o", "--outDir", default ="./", help = "Output Directory")

args = parser.parse_args()

# Project directory and files
frag_path=args.frag_path
meta=args.meta
atac_counts=args.atac_counts
atac_cells=args.atac_cells
atac_peaks=args.atac_peaks
outDir=args.outDir

# Create cisTopic object
atac_counts = mmread(atac_counts)
atac_peaks =  pd.read_csv(atac_peaks)
atac_peaks = [peak.replace("-",":",1) for peak in atac_peaks['x']]  #reformat name
atac_cells =  pd.read_csv(atac_cells)
cell_data =  pd.read_csv(meta,dtype="string")

cistopic_obj = create_cistopic_object(fragment_matrix=atac_counts.tocsr(),
    cell_names=atac_cells['x'].tolist(),
    region_names=atac_peaks,
    tag_cells=False)

# Adding cell information
cistopic_obj.add_cell_data(cell_data)
pickle.dump(
    cistopic_obj,
    open(os.path.join(outDir, "scenicplus_"+args.atac_counts.split('_')[0]+"_cistopic_obj.pkl"), "wb")
)

parser.add_argument("-m", "--memory", required = False, default = '400G', help = "Memory for MALLET_MEMORY")
parser.add_argument("-t", "--taskCpus", required = False, default = 10, help = "CPUS to run")
parser.add_argument("-M", "--mallet", required = False, default ="/container_mallet/bin/mallet", help = "Mallet path, built into container")

args = parser.parse_args()

#tmp dir
tmpDir = './'

# Run models with mallet
os.environ['MALLET_MEMORY'] = '300G'

# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=list(range(5, 40, 5)), 
    n_cpu=int(10),
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=tmpDir,
    save_path=tmpDir,
    mallet_path='/container_mallet/bin/mallet',
)


```

```R
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(dplyr)
library(optparse)
library(org.Hs.eg.db)
library(dendextend)
library(msigdbr,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(fgsea,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(grImport,lib.loc = "/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3/") #local
library(cisTopic)
library(SeuratWrappers)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AUCell)
set.seed(1234)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="8_merged.cnv_clones.SeuratObject.rds", 
              help="Sample input seurat object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)

dat<-subset(dat,merged_assay_clones != "contamination")
dat<-subset(dat,cells=names(which(!is.na(dat$merged_assay_clones))))
dat[["RNA"]]<-JoinLayers(dat[["RNA"]])

#run cistopic on cells with cnv clone assignment
cistopic_outdir=paste0(getwd(),"/cistopic_clones")
system(paste0("mkdir -p ", cistopic_outdir))
out_seurat_object<-paste0(cistopic_outdir,"/","9_merged.cnv_clones_cistopic.SeuratObject.rds")
model_selection_out<-paste0(cistopic_outdir,"/","cistopic.model_selection.pdf")
out_cistopic_obj<-paste0(cistopic_outdir,"/","cistopic.cistopicObject.rds")
umap_out<-paste0(cistopic_outdir,"/","cistopic.umap.pdf")
cistopic_counts_frmt<-dat@assays$ATAC@counts
row.names(cistopic_counts_frmt)<-sub("-", ":", row.names(cistopic_counts_frmt))
sub_cistopic<-cisTopic::createcisTopicObject(cistopic_counts_frmt)
print("Made cistopic object")
sub_cistopic_models<-cisTopic::runModels(sub_cistopic,
  topic=seq(from=10, to=30, by=5),
  nCores=1,
  addModels=FALSE) #using v2 of cistopic (we are only using single core anyway)
sub_cistopic_models<-cisTopic::addCellMetadata(sub_cistopic_models, cell.data =dat@meta.data)
sub_cistopic_models<- cisTopic::selectModel(sub_cistopic_models, type='derivative')
print("Finshed running cistopic")
#Add cell embeddings into seurat
cell_embeddings<-as.data.frame(sub_cistopic_models@selected.model$document_expects)
colnames(cell_embeddings)<-sub_cistopic_models@cell.names
n_topics<-nrow(cell_embeddings)
row.names(cell_embeddings)<-paste0("topic_",1:n_topics)
cell_embeddings<-as.data.frame(t(cell_embeddings))
#Add feature loadings into seurat
feature_loadings<-as.data.frame(sub_cistopic_models@selected.model$topics)
row.names(feature_loadings)<-paste0("topic_",1:n_topics)
feature_loadings<-as.data.frame(t(feature_loadings))
#combined cistopic results (cistopic loadings and umap with seurat object)
cistopic_obj<-CreateDimReducObject(embeddings=as.matrix(cell_embeddings),loadings=as.matrix(feature_loadings),assay="peaks",key="topic_")
print("Cistopic Loading into Seurat")
dat@reductions$cistopic<-cistopic_obj
n_topics<-ncol(Embeddings(dat,reduction="cistopic")) #add scaling for ncount peaks somewhere in here
print("Running UMAP")
dat<-RunUMAP(dat,reduction="cistopic",dims=1:n_topics)
dat <- FindNeighbors(object = dat, reduction = 'cistopic', dims = 1:n_topics ) 
dat <- FindClusters(object = dat, verbose = TRUE, graph.name="peaks_snn", resolution=0.2 ) 
print("Plotting UMAPs")
plt1<-DimPlot(dat,reduction="umap",group.by=c("HBCA_predicted.id"))
pdf(umap_out,width=10)
print(plt1)
dev.off()
saveRDS(sub_cistopic_models,file=out_cistopic_obj)
saveRDS(dat,file=out_seurat_object)


#run titan on cells with cnv clone assignment
model_maker <- function(topics,outDir,cellList,Genes,iterations,burnin,alpha,beta) {
    selected.Model <- lda::lda.collapsed.gibbs.sampler(
      cellList, topics, Genes, 
      num.iterations = iterations, 
      alpha = alpha, eta = beta, 
      compute.log.likelihood = TRUE, 
      burnin = burnin)[-1]
      saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
    }

RPC_calculation<-function(model_file,outDir,data.use,topic_numbers) {
    topic_num <- as.numeric(gsub("[^0-9]+([0-9]+).*", "\\1", model_file))
    topic_numbers <- c(topic_numbers, topic_num)
    model <- readRDS(paste0(outDir, "/", model_file))
    docterMat <- t(as.matrix(data.use))
    docterMat <- as(docterMat, "sparseMatrix")
    topworddist <- normalize(model$topics, byrow = T)
    doctopdist <- normalize(t(model$document_sums), byrow = T)
    perp <- text2vec::perplexity(docterMat, topworddist, doctopdist)
    return(c(topic_num,perp))
}
   
single_sample_titan_generation<- function(Object,
  assay="RNA",
  seed.number=123,
  iterations=500,
  burnin=250,
  alpha_val=50,
  beta_val=0.1,
  nfeat=10000,
  outDir="./TITAN_LDA",
  topic_counts=seq(from=10, to=30, by=5),
  epithelial_only=FALSE,
  sample_in){
  
  print(paste0("Running TITAN on ",sample_in," ..."))
  set.seed(seed.number)
  print(paste0("Setting ",assay," as assay..."))
  if(epithelial_only){
    Object<-subset(Object,sample==sample_in)
    Object<-subset(Object,assigned_celltype %in% c("cancer_luminal_epithelial","luminal_epithelial","basal_epithelial"))
    out_seurat_object<-paste0(outDir,"/",sample_in,".titan_epithelial.SeuratObject.rds")
    out_titan_obj<-paste0(outDir,"/",sample_in,".titan_epithelial.titanObject.rds")
    elbow_out<-paste0(outDir,"/",sample_in,".titan_epithelial.elbow.pdf")
    umap_out<-paste0(outDir,"/",sample_in,".titan_epithelial.umap.pdf")
    model_outdir<-paste0(outDir,"/",sample_in,"_epithelial")
  }
  else {
    Object<-subset(Object,sample==sample_in)
    out_seurat_object<-paste0(outDir,"/",sample_in,".titan.SeuratObject.rds")
    out_titan_obj<-paste0(outDir,"/",sample_in,".titan.titanObject.rds")
    umap_out<-paste0(outDir,"/",sample_in,".titan.umap.pdf")
    elbow_out<-paste0(outDir,"/",sample_in,".titan.elbow.pdf")
    model_outdir<-paste0(outDir,"/",sample_in)
  }

  #skip titan if cell count too low
  if(sum(Object$nCount_RNA>500)<200){
    print("Cell count for RNA seq is too low...")
    saveRDS(Object,file=out_seurat_object)


  }else{
    system(paste0("mkdir -p ",model_outdir))
    Object<-subset(Object,nCount_RNA>500)
    Object[[assay]]<- as(object = Object[[assay]], Class = "Assay") #enforcing v3 assay style

    print(paste0("Finding ",as.character(nfeat)," variable features..."))
    Object <- FindVariableFeatures(
      Object, 
      selection.method = "vst",
      nfeatures = nfeat, 
      assay=assay)

    print(paste0("Setting up data for TITAN LDA..."))
    Object.sparse <- GetAssayData(Object, slot = "data", assay = assay)
    Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assay), ]
    data.use <- Matrix::Matrix(Object.sparse, sparse = T)
    data.use <- data.use * 10 #not sure if needed
    data.use <- round(data.use) #not sure if needed
    data.use <- Matrix::Matrix(data.use, sparse = T)
    sumMat <- Matrix::summary(data.use)
    cellList <- split(as.integer(data.use@i), sumMat$j)
    ValueList <- split(as.integer(sumMat$x), sumMat$j)
    cellList <- mapply(rbind, cellList, ValueList, SIMPLIFY = F)
    Genes <- rownames(data.use)
    cellList <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1, ] + 1]; x})

    print(paste0("Running ",as.character(length(topic_counts)), " topic models..."))
    lda_out<-mclapply(topic_counts, 
      function(x) 
      model_maker(topics=x,outDir=model_outdir,
        cellList=cellList,Genes=Genes,
        iterations=iterations,burnin=burnin,alpha=alpha_val,beta=beta_val), 
      mc.cores = length(topic_counts))

    files <- list.files(path = model_outdir, pattern = "Model_")
    perp_list <- NULL
    topic_numbers <- NULL
    RPC <- NULL
    files <- files[order(nchar(files), files)]

    print(paste0("Generating perplexity estimate for model selection..."))
    perp_out<-as.data.frame(do.call("rbind",
      lapply(files,function(x) RPC_calculation(model_file=x,outDir=model_outdir,data.use=data.use,topic_numbers=topic_numbers))))
    colnames(perp_out) <- c("Topics", "RPC")
    perp_out$Topics<-as.numeric(as.character(perp_out$Topics))
    perp_out$RPC<-as.numeric(as.character(perp_out$RPC))

    rpc_dif<-diff(perp_out$RPC)
    topic_dif<-diff(perp_out$Topics)
    perp_out$perp<-c(NA,abs(rpc_dif)/topic_dif)

    #select topics from model based on elbow
    elbow_topic<-perp_out$Topics[which(min(diff(perp_out$perp),na.rm=T)==diff(perp_out$perp))+1]
    print(paste0("Found ",as.character(elbow_topic), " topics as best topic model based on elbow plot..."))
    plt1 <- ggplot(data = perp_out, aes(x = Topics, y = RPC, group = 1)) + geom_line() + geom_point() +geom_vline(xintercept=elbow_topic,color="red")
    plt2 <- ggplot(data = perp_out, aes(x = Topics, y = perp, group = 1)) + geom_line() + geom_point()+geom_vline(xintercept=elbow_topic,color="red")
    ggsave(plt1/plt2,file=elbow_out)

    top_topics<-elbow_topic #set this up as autoselect based on elbow of plt2 in future
    top_model<-readRDS(paste0(model_outdir, "/", "Model_",as.character(top_topics),"topics.rds"))

    print(paste0("Adding topic model to Seurat Object as lda reduction..."))
    Object <- addTopicsToSeuratObject(model = top_model, Object = Object)
    #GeneDistrubition <- GeneScores(top_model)

    print(paste0("Running UMAP and clustering..."))
    Object<-RunUMAP(Object,reduction="lda",dims=1:top_topics,reduction.name="lda_umap")
    Object <- FindNeighbors(object = Object, reduction = 'lda', dims = 1:top_topics ,graph.name="lda_snn")
    Object <- FindClusters(object = Object, verbose = TRUE, graph.name="lda_snn", resolution=0.2 ) 
    print("Plotting UMAPs...")
    plt1<-DimPlot(Object,reduction="lda_umap",group.by=c("HBCA_predicted.id"))
    ggsave(plt1,file=umap_out,width=10)
    print("Done!")
    saveRDS(Object,out_seurat_object)
    saveRDS(top_model,out_titan_obj)
  }
}

#run pca on cnvs

#correlate
