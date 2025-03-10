# Public Datasets for Comparison 


### Additional Cell Signatures
<!-- Done -->
From https://github.com/yunshun/HumanBreast10X/tree/main/Signatures

```bash
cd /home/groups/CEDAR/mulqueen/ref/embo
#downloaded files from
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/Human-PosSigGenes.RData
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/ImmuneMarkers2.txt
#https://github.com/yunshun/HumanBreast10X/blob/main/Signatures/PAM50.txt
```

### Use EMBO and Swarbrick Paper Cell Types to Define Signatures
Using package genefu for PAM50 pseudobulk assignment.
https://www.bioconductor.org/packages/release/bioc/vignettes/genefu/inst/doc/genefu.html

Running SSpbc method as well

Using https://github.com/StaafLab/sspbc/archive/refs/heads/main.zip for multiple classifications
https://www.nature.com/articles/s41523-022-00465-3#code-availability

#wget https://github.com/StaafLab/sspbc/archive/refs/heads/main.zip
#file located in /home/groups/CEDAR/mulqueen/src/sspbc/sspbc-main/package
#R CMD INSTALL sspbc_1.0.tar.gz


