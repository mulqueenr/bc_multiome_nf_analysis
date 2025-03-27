#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Seurat)
library(Signac)
library(ggplot2)
library(optparse)
library(dplyr)
library(ComplexHeatmap)
library(dendextend)
library(ggdendro)
library(circlize)
library(ggtern)

set.seed(123)
option_list = list(
  make_option(c("-i", "--object_input"), type="character", default="5_merged.geneactivity.SeuratObject.rds", 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dat=readRDS(opt$object_input)


#downloading PAM50 gene list from SCSubtype git repo
if(!file.exists("/ref/NatGen_Supplementary_table_S4.csv")){
  system("wget https://raw.githubusercontent.com/Swarbricklab-code/BrCa_cell_atlas/main/scSubtype/NatGen_Supplementary_table_S4.csv")
  sigdat <- read.csv("NatGen_Supplementary_table_S4.csv",col.names=c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))
} else {
  sigdat <- read.csv("/ref/NatGen_Supplementary_table_S4.csv",col.names=c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"))
}

# read in scsubtype gene signatures
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

#group by scsubtype average scores per cell
scsubtype<-dat_epi@meta.data[grep("SC_Subtype",colnames(dat_epi@meta.data))]
scsubtype<-cbind(scsubtype,dat_epi@meta.data[,c("merge_cluster_50min","ploidy","Diag_MolDiag")])
scsubtype[is.na(scsubtype$merge_cluster_50min),]$merge_cluster_50min<-"NaN"
scsubtype<-scsubtype[startsWith(prefix=c("IDC"),scsubtype$merge_cluster_50min),] 
scsubtype_sum<-as.data.frame(scsubtype %>% 
            group_by(merge_cluster_50min) %>% 
            summarize(ploidy=first(ploidy),diag_moldiag=first(Diag_MolDiag),basal=mean(SC_Subtype_Basal_SC1),
            her2=mean(SC_Subtype_Her2E_SC2),
            lumA=mean(SC_Subtype_LumA_SC3),lumB=mean(SC_Subtype_LumB_SC4)))
row.names(scsubtype_sum)<-scsubtype_sum$merge_cluster_50min

#cluster by all marker genes
sum_da_dend <- as.matrix(dist(scsubtype_sum[4:ncol(scsubtype_sum)])) %>% as.dist %>% hclust(method="ward.D2") %>% as.dendrogram 
k_in<-find_k(sum_da_dend,krange=2:5)
sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
member_split<-cutree(sum_da_dend,k_in$k)
plot_order<-labels(sum_da_dend)

ha = rowAnnotation(
  Diag_MolDiag=scsubtype_sum$diag_moldiag,
  Ploidy=scsubtype_sum$ploidy)

#ha = rowAnnotation(foo = anno_density(m, type = "violin", p = gpar(fill = 1:10)))

pdf("scsubtype_plot.pdf")
Heatmap(scale(scsubtype_sum[4:ncol(scsubtype_sum)],scale=T),cluster_rows=sum_da_dend,  left_annotation=ha,)
dev.off()


#plot per gene for subtypes for clones
gene_class <- as.data.frame(rbind(cbind(as.vector(sigdat[,"Basal_SC"]),c("basal")),
                   cbind(as.vector(sigdat[,"Her2E_SC"]),c("her2")),
                   cbind(as.vector(sigdat[,"LumA_SC"]),c("luma")),
                   cbind(as.vector(sigdat[,"LumB_SC"]),c("lumb"))))
gene_class  <- unique(gene_class [!gene_class [,1] == "",])
colnames(gene_class)<-c("gene","subtype")
#group by scsubtype average scores per cell
Mydata<-AggregateExpression(obj,assay="RNA",return.seurat=T,features=gene_class$gene,group.by="merge_cluster_50min")
Mydata <- ScaleData(Mydata,assay="RNA") #running only on aneuploid epithelial cells
Mydata[["RNA"]] <- as(object = Mydata[["RNA"]], Class = "Assay")
#cluster by all marker genes
sum_da_dend <- as.matrix(dist(Mydata@assays$RNA@scale.data)) %>% as.dist %>% hclust(method="ward.D2") %>% as.dendrogram 
k_in<-find_k(sum_da_dend,krange=2:5)
sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
member_split<-cutree(sum_da_dend,k_in$k)
plot_order<-labels(sum_da_dend)
colnames(Mydata@assays$RNA@scale.data) == row.names(Mydata@meta.data)
ha = rowAnnotation(
  gene_class=gene_class[gene_class$gene %in% row.names(Mydata@assays$RNA@scale.data),]$subtype)

pdf("scsubtype_plot.genes.pdf",height=15)
Heatmap(Mydata@assays$RNA@data, 
  row_order=which(row.names(Mydata@assays$RNA@scale.data) %in% gene_class$gene),  
  left_annotation=ha,
  row_names_gp = grid::gpar(fontsize = 4)
)
dev.off()

#plot per gene for subtypes for samples (PAM50 built on bulk data)

#plot per gene for subtypes for clones
gene_class <- as.data.frame(rbind(cbind(as.vector(sigdat[,"Basal_SC"]),c("basal")),
                   cbind(as.vector(sigdat[,"Her2E_SC"]),c("her2")),
                   cbind(as.vector(sigdat[,"LumA_SC"]),c("luma")),
                   cbind(as.vector(sigdat[,"LumB_SC"]),c("lumb"))))
gene_class  <- unique(gene_class [!gene_class [,1] == "",])
colnames(gene_class)<-c("gene","subtype")
#group by scsubtype average scores per cell
Mydata<-AggregateExpression(dat,assay="RNA",return.seurat=T,features=gene_class$gene,group.by="sample")
Mydata <- ScaleData(Mydata,assay="RNA") #running only on aneuploid epithelial cells
Mydata[["RNA"]] <- as(object = Mydata[["RNA"]], Class = "Assay")

#cluster by all marker genes
sum_da_dend <- as.matrix(dist(Mydata@assays$RNA@scale.data)) %>% as.dist %>% hclust(method="ward.D2") %>% as.dendrogram 
k_in<-find_k(sum_da_dend,krange=2:5)
sum_da_dend<-sum_da_dend %>% set("branches_k_color", k = k_in$k)
member_split<-cutree(sum_da_dend,k_in$k)
plot_order<-labels(sum_da_dend)
colnames(Mydata@assays$RNA@scale.data) == row.names(Mydata@meta.data)
ha = rowAnnotation(
  gene_class=gene_class[gene_class$gene %in% row.names(Mydata@assays$RNA@scale.data),]$subtype)

pdf("scsubtype_plot.genes_by_sample.pdf",height=15)
Heatmap(Mydata@assays$RNA@scale.data, 
  row_order=which(row.names(Mydata@assays$RNA@scale.data) %in% gene_class$gene),  
  left_annotation=ha,
  row_names_gp = grid::gpar(fontsize = 4)
)
dev.off()




````


