library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cisTopic)
library(TITAN)
library(optparse)
option_list = list(
  make_option(c("-c", "--cistopic"), type="character", default=NULL, 
              help="List of sample cisTopic RDS files", metavar="character")
  make_option(c("-t", "--titan"), type="character", default=NULL, 
              help="List of sample TITAN RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#cistopic=readRDS(opt$cistopic)
#titan=readRDS(opt$titan)

cistopic_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/cistopic_objects"
titan_path="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/titan_objects"
cistopic_objs<-list.files(path=cistopic_path,pattern="*CisTopicObject.rds$",full.names=TRUE)
titan_objs<-list.files(path=titan_path, pattern="*TITANObject.rds$",full.names=TRUE)


titan_top_genes<-function(x){
  titan<-readRDS(x)
  TitanTopicGenes <- as.list(as.data.frame(TopTopicGenes(titan, ngenes = 50)))
  return(TitanTopicGenes)
}

titan_topic_genes<-lapply(titan_objs,function(x) titan_top_genes(x)) #can parallelize this probably
names(titan_topic_genes)<-unlist(lapply(titan_objs,function(x) strsplit(basename(x),"[.]")[[1]][1]))

#this is defintely not the most efficient way to do this
out<-list()
for(a in 1:length(titan_topic_genes)){
  for(b in 1:length(titan_topic_genes[[a]])){
    for(c in 1:length(titan_topic_genes)){
      for(d in 1:length(titan_topic_genes[[c]])){
        oneset<-titan_topic_genes[[a]][b]
        twoset<-titan_topic_genes[[c]][d]
        out<-append(out,c(names(titan_topic_genes)[a],names(titan_topic_genes[[a]])[b],names(titan_topic_genes)[c],names(titan_topic_genes[[c]])[d],sum(oneset %in% twoset)))
}}}}


head(out)
lapply(1:length(titan_topic_genes[[1]]))
Reduce(intersect, titan_topic_genes[[1]])

x=titan_objs[1]
#then process this as nmf?


cistopic <- getRegionsScores(cistopic, method='NormTop', scale=TRUE)


cistopic <- binarizecisTopics(cistopic, thrP=0.975, plot=TRUE)


cistopic <- GREAT(cistopic, genome='hg38', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
#https://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_10X_workflow.html
#Enrichment of epigenomic signatures in the cells

#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#system("mkdir -p output/cisTopics_asBW")
#getBigwigFiles(cistopic, path='output/cisTopics_asBW', seqlengths=seqlengths(txdb))