#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif
library(optparse)
library(SoupX)


option_list = list(
  make_option(c("-o", "--out_name"), type="character", default=NULL, 
              help="Out name prefix", metavar="character"),
    make_option(c("-d", "--out_dir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

args = commandArgs(trailingOnly=TRUE)
outname=opt$out_name
wd=paste0("./",opt$out_dir,"/","outs")

sc = load10X(wd)
sc = autoEstCont(sc,tfidfMin=0.5,forceAccept=TRUE) #1 is default
out = adjustCounts(sc)
saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
print(paste("Finished:",outname))
