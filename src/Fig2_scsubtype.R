
```bash
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_nmf.sif"
singularity shell \
--bind /home/groups/CEDAR/scATACcnv/Hisham_data \
--bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
```


```R

library(Seurat)
library(Signac)
library(ggplot2)
library(optparse)
set.seed(123)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.clone_annot.passqc.SeuratObject.rds"
dat=readRDS(opt$object_input)
dat_epi<-subset(dat,cell=row.names(dat@meta.data)[dat$reclust %in% c("luminal_epithelial","basal_epithelial")])
dat_epi<-JoinLayers(dat_epi,assay="RNA")

#downloading PAM50 gene list from SCSubtype git repo
if(!file.exists("NatGen_Supplementary_table_S4.csv")){
  system("wget https://raw.githubusercontent.com/Swarbricklab-code/BrCa_cell_atlas/main/scSubtype/NatGen_Supplementary_table_S4.csv")
}
# read in scsubtype gene signatures
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv",col.names=c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

module_feats<-list()
module_feats[["Basal_SC"]]=as.vector(sigdat[,"Basal_SC"])
module_feats[["Her2E_SC"]]=as.vector(sigdat[,"Her2E_SC"])
module_feats[["LumA_SC"]]=as.vector(sigdat[,"LumA_SC"])
module_feats[["LumB_SC"]]=as.vector(sigdat[,"LumB_SC"])
module_feats<-lapply(module_feats,function(x) {x[x!=""]})

dat_epi<-AddModuleScore(dat_epi,module_feats,
  assay = "RNA",
  name = paste0("SC_Subtype_",c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC")),
  search = TRUE)

scsubtype<-dat_epi@meta.data[grep("SC_Subtype",colnames(dat_epi@meta.data))]
dat<-AddMetaData(dat,scsubtype)
saveRDS(dat,file=opt$object_input)

````