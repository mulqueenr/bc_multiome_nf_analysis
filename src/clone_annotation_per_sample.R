
```bash
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome \
--bind /home/groups/CEDAR/scATACcnv/Hisham_data \
--bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(Seurat)
library(Signac)
library(optparse)

option_list = list(
  make_option(c("-c", "--cistopic"), type="character", default=NULL, 
              help="List of sample cisTopic RDS files", metavar="character"),
  make_option(c("-t", "--titan"), type="character", default=NULL, 
              help="List of sample TITAN RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3")
opt$object_input="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects/merged.geneactivity.SeuratObject.rds"
dat=readRDS(opt$object_input)
#correct dat cell names
dat<-RenameCells(dat,new.names=paste0(dat$sample,"_",unlist(lapply(strsplit(dat$cellID,"_"),"[",3))))

cluster_idents_samples<-lapply(
    list.files(path="/home/groups/CEDAR/scATACcnv/Hisham_data/new_seq/1MB/consensus",pattern="*_consensus6_clust.csv",full.names=T), 
    function(x){
    in_clus=read.csv(x)
    in_clus$sample=substr(basename(x),1,nchar(basename(x))-nchar("_consensus6_clust.csv"))
    in_clus$cellID=paste(in_clus$sample,in_clus$X,sep="_")
    in_clus$cellID<-gsub('\\.', '-', in_clus$cellID)
    row.names(in_clus)<-in_clus$cellID
    return(in_clus)}
)

cluster_in<-do.call("rbind",cluster_idents_samples)
dat<-AddMetaData(dat,cluster_in)
saveRDS(dat,opt$object_input)

diploid_clusters<-read.csv("NormalClusters.csv",col.names=c("sample","diploid_cluster")) #copy and pasted from a TM
met<-dat@meta.data
met$cnv_ploidy<-NA
met[!is.na(met$cluster),]$cnv_ploidy<-"aneuploid"

for(x in 1:nrow(diploid_clusters)){
    tmp<-diploid_clusters[x,]
    if(!is.na(tmp$diploid_cluster)){
        if(is.integer(tmp$diploid_cluster)){
            met[(met$sample %in% tmp$sample) && (met$cluster %in% tmp$diploid_cluster),]$cnv_ploidy<-"diploid"
        } else if(tmp$diploid_cluster=="All"){
            met[(met$sample==tmp$sample),]$cnv_ploidy<-"diploid"
        }
    }
}

cnv_ploidy<-setNames(nm=row.names(met),met$cnv_ploidy)
dat<-AddMetaData(dat,cnv_ploidy,col.name="cnv_defined_ploidy")
saveRDS(dat,opt$object_input)
