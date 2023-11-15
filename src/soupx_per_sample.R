library(SoupX)
args = commandArgs(trailingOnly=TRUE)
outname=args[1]
wd=paste0("./",args[2],"/","outs")

sc = load10X(wd)
sc = autoEstCont(sc,tfidfMin=0.5,forceAccept=TRUE) #1 is default
out = adjustCounts(sc)
saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
print(paste("Finished:",outname))


