library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")

#Generate tile references
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38" ,tileWidth = 1e6,outputDir = "/home/groups/CEDAR/mulqueen/ref/copyscat")

library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
setwd("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2")
args = commandArgs(trailingOnly=TRUE)


#initialize the environment
initialiseEnvironment(genomeFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_chrom_sizes.tsv",
                      cytobandFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=500,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#Set up copyscAT Loop per sample
copyscAT_per_sample<-function(x,prediction="EMBO",knn_in=FALSE,cores=1){
  if(x %in% 1:12){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220414_multiome_phase1/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else if(x %in% 13:20){
    sample_name<-paste0("sample_",x)
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs")
    outname<-paste0("sample_",x)
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220715_multiome_phase2/sample_",x,"/outs/sample_",x,".QC.SeuratObject.rds")
  }else{
    sample_name<-x
    wd<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs")
    outname<-x
    file_in<-paste0("/home/groups/CEDAR/mulqueen/projects/multiome/220111_multi/",x,"/outs/",x,".QC.SeuratObject.rds")
  }
  if (knn_in==TRUE){
  knn_list<-read.table(paste0("/home/groups/CEDAR/scATACcnv/Hisham_data/bed_files/WGS_eval/knn/",sample_name,"_knn5_neighbors.csv"),
    sep=",",header=T)
  knn_list<-as.data.frame(apply(knn_list, 2, function(y) gsub("[.]", "-", y)))
  print("Knn_In Found")
  #knn list is csv format <rowid><cell><neighbor1><neighbor2><neighbor3><neighbor4>
  }
  dat<-readRDS("phase2.QC.filt.SeuratObject.rds") #use QC controlled bulk seurat object as input
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
  obj_name=basename(file_in)
  dir_in=dirname(file_in)
  if (knn_in == TRUE){
  system(paste0("mkdir ",dir_in,"/copyscat_knn"))
  print("Knn_In Found")
  } else {
  system(paste0("mkdir ",dir_in,"/copyscat"))
  }
  #do python script preprocessing (basically just count fragments per window per cell)
  system(paste0("python /home/groups/CEDAR/mulqueen/ref/copyscat/process_fragment_file.py ",
  " -i ",dir_in,"/atac_fragments.tsv.gz",
  " -o ",dir_in,"/copyscat/copyscat.1mb.tsv",
  " -b ","1000000",
  " -f ","500",
  " -g ","/home/groups/CEDAR/mulqueen/ref/copyscat/hg38_chrom_sizes.tsv",
  " -c ",dir_in,"/metadata.tsv")) #modification takes in metadata table to filter cells by name, ignores -f flag
  if (knn_in==TRUE){
  setOutputFile(paste0(dir_in,"/copyscat_knn"),"copyscat_out_knn")
  } else {
  setOutputFile(paste0(dir_in,"/copyscat"),"copyscat_out")
  }

  #PART 1: INITIAL DATA NORMALIZATION
  scData<-readInputTable(paste0(dir_in,"/copyscat/copyscat.1mb.tsv"))
  #here is an if else, one python script also accounts for metacell merged cells, other is strictly single cell
  if(knn_in==TRUE){
    scData2<-as.data.frame(do.call("rbind",mclapply(1:nrow(knn_list), function(x){colSums(scData[row.names(scData) %in% unlist(knn_list[x,]),])},mc.cores=cores)))
    row.names(scData2)<-knn_list$cell
    scData<-scData2
  }

  #collapse into chromosome arm level
  summaryFunction<-cutAverage
  scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)
  scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 

  #PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
  #ALTERNATE METHOD FOR CNV CALLING (with normal cells as background)
  #Using same normal cell selection as used for CASPER and InferCNV
  dat$cnv_ref<-"FALSE"
  if(prediction=="EMBO"){
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("epithelial")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
    }else{
  dat@meta.data[dat$predicted.id %in% c("Endothelial","B-cells","Myeloid","Plasmablasts","PVL","T-cells","CAFs"),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  dat<-subset(dat,predicted.id %in% c("Cancer Epithelial","Normal Epithelial","Endothelial","T-cells","B-cells","Myeloid","Plasmablasts","PVL","CAFs"))
  } 
  control<-names(dat$cnv_ref == "TRUE") #pulling this from the inferCNV function

  #compute central tendencies based on normal cells only
  colnames(scData_collapse)<-gsub(outname,"",colnames(scData_collapse))
  control<-gsub(paste0(outname,"_"),"",control)
  control <- control[control %in% colnames(scData_collapse)] #filter control list to control cells that survived filter
  median_iqr <- computeCenters(scData_collapse %>% select(chrom,control),summaryFunction=summaryFunction)
  #setting medianQuantileCutoff to -1 and feeding non-neoplastic barcodes in as normalCells can improve accuracy of CNV calls
  candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,
    useDummyCells = FALSE,
    propDummy=0.25,
    minMix=0.01,
    deltaMean = 0.03,
    deltaBIC2 = 0.25,
    bicMinimum = 0.1,
    subsetSize=50,
    fakeCellSD = 0.09,
    uncertaintyCutoff = 0.65,
    summaryFunction=summaryFunction,
    maxClust = 4,
    mergeCutoff = 3,
    IQRCutoff = 0.25,
    medianQuantileCutoff = -1,
    normalCells=control) 
  candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)

  if(knn_in==TRUE){
  saveRDS(candidate_cnvs_clean,file=paste0(dir_in,"/copyscat_knn/",sample_name,"copyscat_cnvs_matrix_knn.rds"))
  }else{saveRDS(candidate_cnvs_clean,file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs_matrix.rds"))}

  #to save this data you can use annotateCNV4 as per usual, using normal barcodes
  final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, expectedNormals=control, saveOutput=TRUE,
    outputSuffix = "clean_cnv_b2",sdCNV = 0.6,filterResults=FALSE,filterRange=0.4,minAlteredCellProp = 0.5)

  if(knn_in==TRUE){
  saveRDS(final_cnv_list,file=paste0(dir_in,"/copyscat_knn/",sample_name,"copyscat_cnvs_knn.rds"))
  }else{saveRDS(final_cnv_list,file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs.rds"))}

  print(paste("Finished sample",sample_name))
}

copyscAT_per_sample(x=as.character(args[1]),knn=FALSE)
copyscAT_per_sample(x=as.character(args[1]),knn=TRUE)

#lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2", "RM_3","RM_4",4,10,12),copyscAT_per_sample)
#lapply(c(1,3,5,6,7,8,9,11,15,16,19,20,"RM_1","RM_2", "RM_3","RM_4",4,10,12),function(x) copyscAT_per_sample(x,knn_in=TRUE,cores=5))
lapply(c(4,11,10,12),function(x) copyscAT_per_sample(x,knn_in=TRUE,cores=5))
#copyscat_dat<-readRDS(file=paste0(dir_in,"/copyscat/",sample_name,"copyscat_cnvs_matrix.rds"))
