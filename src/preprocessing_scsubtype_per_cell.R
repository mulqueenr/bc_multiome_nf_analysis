```bash
module load singularity
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome \
$sif
```

```R
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(ggtern)
library(patchwork)
library(optparse)
set.seed(1234)

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
dat_epi<-subset(dat,cnv_defined_ploidy=="aneuploid" & Diagnosis=="IDC")
pseudo_clones <- AggregateExpression(dat_epi, assays = "RNA", return.seurat = T, group.by = c("sample","cluster"))

system("wget https://raw.githubusercontent.com/Swarbricklab-code/BrCa_cell_atlas/main/scSubtype/NatGen_Supplementary_table_S4.csv")
# read in scsubtype gene signatures
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv",col.names=c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

#Read in the single cell RDS object as 'Mydata'
Mydata<-pseudo_clones
DefaultAssay(Mydata)<-"RNA"
Mydata[["RNA"]] <- as(object = Mydata[["RNA"]], Class = "Assay")
Mydata <- ScaleData(Mydata, features=temp_allgenes,assay="RNA") #running only on aneuploid epithelial cells
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data) #using RNA (not SoupXRNA corrected read counts)


#calculate mean scsubtype scores
outdat <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
  outdat[i,]<-as.numeric(temp)
}

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

##Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames
finalm.sweep.t$sample <- unlist(lapply(strsplit(row.names(finalm.sweep.t),"_"),"[",1))
finalm.sweep.t$sample<-unlist(lapply(finalm.sweep.t$sample, function(x) gsub("-","_",x)))

##Writing out output files (rownames remain the same for both)
write.table(finalm.sweep.t, "SCSubtype_Scores.txt", sep="\t")
write.table(Finalnames, "SCSubtype_calls.txt", sep="\t")

met<-dat@meta.data
met<-met[met$Diagnosis %in% c("IDC"),]
met$group<-paste(met$Diagnosis,met$Mol_Diagnosis)
met<-met[!duplicated(met$sample),]
sctype<-merge(met,finalm.sweep.t,by="sample")
table(sctype$group,sctype$SCSubtypeCall)

plt<-ggtern(data=sctype,aes(x=LumA_SC.y, z=LumB_SC.y,y=Basal_SC.y,color=group)) +
  geom_point(size=0.2,alpha=0.5) + theme_minimal()
plt<-ggtern(data=sctype,aes(x=LumA_SC.y, z=LumB_SC.y,y=Basal_SC.y,color=SCSubtypeCall.y)) +
  geom_point(size=0.2,alpha=0.5) + theme_minimal()
gsave(plt,file="scsubtype_ternary.pdf")



dat_epi<-subset(dat,cnv_defined_ploidy=="aneuploid")
module_scores<-AddModuleScore(dat_epi,features=module_feats,assay="SoupXRNA",search=TRUE,name=names(module_feats)) #use add module function to add cell scores
module_scores<-module_scores@meta.data[seq(ncol(module_scores@meta.data)-(length(module_feats)-1),ncol(module_scores@meta.data))]
colnames(module_scores)<-names(module_feats) #it adds a number at the end to each name by default, which I don't like

dat<-AddMetaData(dat,metadata=module_scores)


sample_names<-paste(unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(1))),
  unlist(lapply(strsplit(colnames(dat[["RNA"]]@counts),"_"),"[",c(2))),sep="_")
counts<-as.data.frame(t(dat[["RNA"]]@counts)) 
counts<-cbind(counts,sample_names)
counts<-as.data.frame(counts %>% group_by(sample_names) %>% summarize_all(funs(sum)))
row.names(counts)<-counts$sample_name
counts<-counts[,2:ncol(counts)]
counts<-counts[,colSums(counts)>0]
dat_in<-counts
dat_in<-dat_in[!(row.names(dat_in) %in% c("RM_4","sample_15","sample_19")),] #exclude NAT samples
dat_in<-as.data.frame(t(dat_in))

#set up matrix by unique entrez gene names
dat_in<-dat_in[!duplicated(mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')),]
dat_in<-dat_in[!isNA(mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')),]
row.names(dat_in)<-mapIds(org.Hs.eg.db, row.names(dat_in), 'ENTREZID', 'SYMBOL')
myresults <- applySSP(gex=as.matrix(dat_in), id=row.names(dat_in), ssp.name="ssp.pam50",id.type="EntrezGene",report=TRUE)

dat_pam50<-setNames(nm=row.names(dat@meta.data),myresults[match(dat@meta.data$sample,row.names(myresults)),1])
dat<-AddMetaData(dat,dat_pam50,col.name="pseudobulk_sspbc_PAM50")
saveRDS(dat,file="phase2.QC.filt.SeuratObject.rds")
````
