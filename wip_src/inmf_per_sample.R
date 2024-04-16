library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
library(TITAN)
library(ComplexHeatmap)
library(reshape2)
library(optparse)
library(circlize)

option_list = list(
  make_option(c("-c", "--cistopic"), type="character", default=NULL, 
              help="List of sample cisTopic RDS files", metavar="character"),
  make_option(c("-t", "--titan"), type="character", default=NULL, 
              help="List of sample TITAN RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#cistopic=readRDS(opt$cistopic)
#titan=readRDS(opt$titan)

cistopic_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/cistopic_objects"
titan_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/titan_objects"
cistopic_objs<-list.files(path=cistopic_path,pattern="*cistopic.cistopicObject.rds$",full.names=TRUE)
titan_objs<-list.files(path=titan_path, pattern="*titan.titanObject.rds$",full.names=TRUE)

cistopic_epi_objs<-list.files(path=cistopic_path,pattern="*cistopic_epithelial.cistopicObject.rds$",full.names=TRUE)
titan_epi_objs<-list.files(path=titan_path, pattern="*titan_epithelial.titanObject.rds$",full.names=TRUE)

titan_top_genes<-function(x){
  titan<-readRDS(x)
  TitanTopicGenes <- as.list(as.data.frame(TopTopicGenes(titan, ngenes = 50)))
  return(TitanTopicGenes)
}

titan_topic_genes<-lapply(titan_objs,function(x) titan_top_genes(x)) #can parallelize this probably
names(titan_topic_genes)<-unlist(lapply(titan_objs,function(x) strsplit(basename(x),"[.]")[[1]][1]))


compare_topics<-function(x_name,y_name){
  x=titan_topic_genes[x_name]
  y=titan_topic_genes[y_name]
  lapply(names(titan_topic_genes[[x_name]]),function(xi){
    lapply(names(titan_topic_genes[[y_name]]),function(yj){
      c(x_name,y_name,xi,yj,length(intersect(unlist(titan_topic_genes[[x_name]][xi]),unlist(titan_topic_genes[[y_name]][yj]))))
      })
      }) 
}

out<-lapply(names(titan_topic_genes),function(x){lapply(names(titan_topic_genes),function(y){compare_topics(x,y)})})
out<-as.data.frame(do.call("rbind",do.call("rbind",do.call("rbind",do.call("rbind",out))))) #undo the crazy nested list I made
colnames(out)<-c("dat1","dat2","topic_dat1","topic_dat2","overlap")
out$overlap<-as.numeric(out$overlap)
out_mat<-as.data.frame(dcast(dat1+topic_dat1~dat2+topic_dat2,data=out,value.var="overlap",max))
row.names(out_mat)<-unique(paste(out$dat1,out$topic_dat1))#check that row is correct for this
out_mat<-out_mat[,3:ncol(out_mat)]

ha = rowAnnotation(diag = unlist(lapply(strsplit(row.names(out_mat),split="_"),"[",1)),
    col = list(diag = c("DCIS" = "#CCCCCC", "IDC" = "#FF6699", "ILC" = "#006633", "NAT" = "#99CCFF")))

pdf("test.pdf")
Heatmap(out_mat, name = "TITAN Topics", 
    col=colorRamp2(c(0,50),c("white", "purple")),
    row_names_gp = gpar(fontsize = 0.5),
    column_names_gp = gpar(fontsize =0.5),
    width=20,height=20, right_annotation = ha
)
dev.off()

    ,
    #column_split=unlist(lapply(strsplit(colnames(out_mat),split="\t"),"[",1)),

system("slack -F test.pdf ryan_todo")



cistopic <- getRegionsScores(cistopic, method='NormTop', scale=TRUE)
cistopic <- binarizecisTopics(cistopic, thrP=0.975, plot=TRUE)
cistopic <- GREAT(cistopic, genome='hg38', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
#https://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_10X_workflow.html
#Enrichment of epigenomic signatures in the cells

#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#system("mkdir -p output/cisTopics_asBW")
#getBigwigFiles(cistopic, path='output/cisTopics_asBW', seqlengths=seqlengths(txdb))