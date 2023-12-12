library(Seurat)
library(Signac)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)

args = commandArgs(trailingOnly=TRUE)
dat=readRDS(args[1]) #dat=readRDS("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds")
sample_arr=as.numeric(as.character(args[2])) #sample_arr=as.numeric(as.character(8)) 
proj_dir=args[3] #proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
ref_dir=paste0(proj_dir,"/ref/copyscat/")

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
copyscAT_per_sample<-function(dat, outname){
  dat<-subset(dat,sample==outname) #subset data to sample specified by x and outname
  dir_in=paste0(proj_dir,"/cellranger_data/second_round/",outname,"/outs")
  write.table(as.data.frame(dat@meta.data),col.names=T,row.names=T,sep="\t",file=paste0(dir_in,"/metadata.tsv"))
  bam_location<-paste0(dir_in,"/gex_possorted_bam.bam")

  #do python script preprocessing (basically just count fragments per window per cell)
  system(paste0(
    "python ",ref_dir,"/process_fragment_file.py ",
    " -i ",dir_in,"/atac_fragments.tsv.gz",
    " -o ",dir_in,"/copyscat.1mb.tsv",
    " -b ","1000000",
    " -f ","500",
    " -g ",ref_dir,"hg38_chrom_sizes.tsv",
    " -c ",dir_in,"/metadata.tsv")) #modification takes in metadata table to filter cells by name, ignores -f flag
    
  #PART 1: INITIAL DATA NORMALIZATION
  scData<-readInputTable(paste0(dir_in,"/copyscat.1mb.tsv"))

  #collapse into chromosome arm level
  summaryFunction<-cutAverage
  scData_k_norm <- normalizeMatrixN(scData,
    logNorm = FALSE,
    maxZero=3000,
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
  dat@meta.data[!(dat$HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell")),]$cnv_ref<-"TRUE" #set cnv ref by cell type
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
  
  candidate_cnvs_clean<-clusterCNV(
    initialResultList = candidate_cnvs,
    medianIQR = candidate_cnvs[[3]],
    minDiff=1.0) #= 1.5)

  saveRDS(candidate_cnvs_clean,file=paste0(outname,"_cnvs_matrix.copyscat.rds"))

  #to save this data you can use annotateCNV4 as per usual, using normal barcodes
  final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, 
    expectedNormals=control, 
    saveOutput=TRUE,
    outputSuffix = "clean_cnv_b2",
    sdCNV = 0.6,
    filterResults=FALSE,
    filterRange=0.4,
    minAlteredCellProp = 0.5)

  saveRDS(final_cnv_list,file=paste0(outname,"_copyscat_cnvs.copyscat.rds"))
}

copyscAT_per_sample(dat=dat,outname=unique(dat$sample)[sample_arr])


#proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
#seurat_obj="${proj_dir}/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"
#Rscript ${proj_dir}/src/copyscat_per_sample.R $seurat_obj 1 $proj_dir


"""
#alternative on single node (three jobs at once):
arr_in=$(seq 1 19)
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
src_dir=${proj_dir}"/src"
obj="${proj_dir}/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"
cd ${proj_dir}/nf_analysis/cnv_analysis/copyscat
parallel -j 1 Rscript ${src_dir}/copyscat_per_sample.R $obj {} $proj_dir ::: $arr_in

"""