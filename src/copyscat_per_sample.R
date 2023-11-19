library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
dat=readRDS(args[1])
ref_dir=paste0(args[2],"/copyscat/")
cellranger_in=args[3]
outdir=args[4]

knn_supplied=FALSE

if(length(args)>4){
  knn_in=args[5]
  knn_supplied=TRUE
}

#Generate tile references
#generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38" ,tileWidth = 1e6,outputDir = ref_dir)

#initialize the environment
initialiseEnvironment(genomeFile=paste0(ref_dir,"hg38_chrom_sizes.tsv"),
                      cytobandFile=paste0(ref_dir,"hg38_1e+06_cytoband_densities_granges.tsv"),
                      cpgFile=paste0(ref_dir,"/hg38_1e+06_cpg_densities.tsv"),
                      binSize=1e6,
                      minFrags=500,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#Set up copyscAT Loop per sample
copyscAT_per_sample<-function(outname,
  knn_in=NULL,
  knn_supplied=FALSE){

  if (knn_supplied == TRUE){
  knn_list<-read.table(knn_in,sep=",",header=T)
  knn_list<-as.data.frame(apply(knn_list, 2, function(y) gsub("[.]", "-", y)))
  #knn list is csv format <rowid><cell><neighbor1><neighbor2><neighbor3><neighbor4>
  }

  if (knn_supplied == TRUE){
  outdir_sample=paste0(cellranger_in,"/copyscat_knn/")
  out_prefix="copyscat_knn"
  } else {
  outdir_sample=paste0(cellranger_in,"/copyscat/")
  out_prefix="copyscat"
  }
  system(paste0("mkdir -p",outdir_sample))

  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname


  #do python script preprocessing (basically just count fragments per window per cell)
  system(paste0("python ",ref_dir,"/process_fragment_file.py ",
  " -i ",cellranger_in,"/atac_fragments.tsv.gz",
  " -o ",outdir_sample,"/copyscat.1mb.tsv",
  " -b ","1000000",
  " -f ","500",
  " -g ",ref_dir,"hg38_chrom_sizes.tsv",
  " -c ",cellranger_in,"/metadata.tsv")) #modification takes in metadata table to filter cells by name, ignores -f flag
  
  setOutputFile(outdir_sample,paste0(out_prefix,"_out"))

  #PART 1: INITIAL DATA NORMALIZATION
  scData<-readInputTable(paste0(outdir_sample,"/copyscat.1mb.tsv"))

  #here is an if else, one python script also accounts for metacell merged cells, other is strictly single cell
  if(knn_in==TRUE){
    scData2<-as.data.frame(do.call("rbind",mclapply(1:nrow(knn_list), function(x){colSums(scData[row.names(scData) %in% unlist(knn_list[x,]),])},mc.cores=cores)))
    row.names(scData2)<-knn_list$cell
    scData<-scData2
  }

  #collapse into chromosome arm level
  summaryFunction<-cutAverage
  scData_k_norm <- normalizeMatrixN(scData,
    logNorm = FALSE,
    maxZero=2000,
    imputeZeros = FALSE,
    blacklistProp = 0.8,
    blacklistCutoff=125,
    dividingFactor=1,
    upperFilterQuantile = 0.95)

  scData_collapse<-collapseChrom3N(
    scData_k_norm,
    summaryFunction=summaryFunction,
    binExpand = 1,
    minimumChromValue = 100,
    logTrans = FALSE,
    tssEnrich = 1,
    logBase=2,
    minCPG=300,
    powVal=0.73) 

  #PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
  #ALTERNATE METHOD FOR CNV CALLING (with normal cells as background)
  #Using same normal cell selection as used for CASPER and InferCNV
  dat$cnv_ref<-"FALSE"
  dat@meta.data[!(dat$EMBO_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
  control<-names(dat$cnv_ref == "TRUE") 

  #compute central tendencies based on normal cells only
  colnames(scData_collapse)<-gsub(outname,"",colnames(scData_collapse))
  control<-gsub(paste0(outname,"_"),"",control)
  control <- control[control %in% colnames(scData_collapse)] #filter control list to control cells that survived filter
  median_iqr <- computeCenters(scData_collapse %>% select(chrom,control),summaryFunction=summaryFunction)
  
  #setting medianQuantileCutoff to -1 and feeding non-neoplastic barcodes in as normalCells can improve accuracy of CNV calls
  candidate_cnvs<-identifyCNVClusters(scData_collapse,
    median_iqr,
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
  
  candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,
  medianIQR = candidate_cnvs[[3]],
  minDiff=1.0) #= 1.5)

  saveRDS(candidate_cnvs_clean,
    file=paste0(outdir_sample,out_prefix,"_"outname,"_cnvs_matrix_knn.rds"))

  #to save this data you can use annotateCNV4 as per usual, using normal barcodes
  final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, 
    expectedNormals=control, 
    saveOutput=TRUE,
    outputSuffix = "clean_cnv_b2",
    sdCNV = 0.6,
    filterResults=FALSE,
    filterRange=0.4,
    minAlteredCellProp = 0.5)

  saveRDS(final_cnv_list,
    file=paste0(outdir_sample,out_prefix,"_"outname,"_copyscat_cnvs.rds"))
  print(paste("Finished sample",sample_name))
}

mclapply(unique(dat$sample),function(x) copykat_per_sample(dat=dat,outname=x),mc.cores=10)
