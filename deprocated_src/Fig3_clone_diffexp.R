```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```

```R

library(Signac)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(grid)
library(dplyr) 
library(ggplot2)
library(ggrepel)
library(patchwork)
library(presto)
library(seriation)
library(org.Hs.eg.db)
library(ggtern)
library(dendextend)
library(optparse)
library(msigdbr) #local
library(fgsea) #local
library(GeneNMF,lib="/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3") #local

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.clone_annot.passqc.SeuratObject.rds"
dat=readRDS(opt$object_input)
dat[["RNA"]]<-as(object = dat[["RNA"]], Class = "Assay")


Idents(dat)<-dat$RNA_cluster

#H is hallmark cancer #C1 is cytoband #C3 is TF targets
run_gsea_enrichment <- function(obj,sample="IDC_12",assay="RNA",category="H") {
    obj<-subset(dat,cells=row.names(dat@meta.data)[dat$sample==sample])
    de_genes<-FindAllMarkers(obj,assay=assay,group_by="RNA_cluster")

    gsea<-lapply(unique(de_genes$cluster), 
    function(x){
        program_genes=unlist(de_genes[de_genes$cluster==x,]$gene)
        out<-runGSEA(program_genes, universe=row.names(obj@assays[[assay]]), category = category) 
        out<-as.data.frame(out)
        out$cluster<-x
        return(out)})
    
    gsea_out<-do.call("rbind",gsea)
    pltdat<-gsea_out %>% group_by(cluster) %>% slice_max(order_by = -padj, n = 5)
    plt<-ggplot(pltdat,aes(x=cluster,y=pathway))+
    geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
    theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(plt,file=paste(sample,"clone_cnv",category,"RNA.pdf",sep="_"),width=20)
}

run_gsea_enrichment(obj=dat,sample="IDC_12",assay="RNA",category="H")
run_gsea_enrichment(obj=dat,sample="IDC_12",assay="RNA",category="C1")
run_gsea_enrichment(obj=dat,sample="IDC_12",assay="RNA",category="C3")

run_gsea_enrichment(obj=dat,sample="ILC_01",assay="RNA",category="H")
run_gsea_enrichment(obj=dat,sample="ILC_01",assay="RNA",category="C1")
run_gsea_enrichment(obj=dat,sample="ILC_01",assay="RNA",category="C3")