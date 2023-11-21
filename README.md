# Nextflow processing

This repository contains the nextflow pipeline and accessory scripts for processing samples from standard cellranger-arc output through final analysis for the breast cancer multiome CEDAR project. 

The file /nextflow_version/bc_multiome.nf.groovy contains the nextflow pipeline and calls R and bash scripts from src/
Processing of public data, used in cell label transfer, if contained in the processing_public_data.md.

The code below details the full processing steps for bcl files to cellranger-arc output. Which is subsequently used for nextflow input.

Example run of nextflow:
```bash
#module load nextflow #on exacloud
#cd /home/groups/CEDAR/mulqueen/bc_multiome

nextflow bc_multiome.nf.groovy \
  -with-dag bc_multiome.flowchart.pdf \
  -with-report bc_multiome.report.html \
  --merged_bed merged_500bp.bed \
  -resume
```

