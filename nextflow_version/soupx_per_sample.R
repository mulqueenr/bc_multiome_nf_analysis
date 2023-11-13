library(SoupX)
args = commandArgs(trailingOnly=TRUE)
outname=args[1]
wd=args[2]

sc = load10X(wd)
sc = autoEstCont(sc,tfidfMin=1) #1 is default
out = adjustCounts(sc)
saveRDS(out,paste0(wd,"/soupx_corrected_counts.rds"))
print(paste("Finished:",outname))


