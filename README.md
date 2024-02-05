# Nextflow processing

This repository contains the nextflow pipeline and accessory scripts for processing samples from standard cellranger-arc output through final analysis for the breast cancer multiome CEDAR project. 

The file /nextflow_version/bc_multiome.nf.groovy contains the nextflow pipeline and calls R and bash scripts from src/

The code in processing_10x_data.md details the full processing steps for bcl files to cellranger-arc output. Which is subsequently used for nextflow input.

Processing of public data, used in cell label transfer, if contained in the processing_public_data.md.

Example run of nextflow:
```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome #move to project directory
git clone https://github.com/mulqueenr/bc_multiome_nf_analysis.git #pull github repo
module load singularity #load singularity on exacloud
module load nextflow #load nextflow on exacloud

#run nextflow with defaults
nextflow bc_multiome/bc_multiome.nf.groovy \
-with-dag bc_multiome.flowchart.png \
-with-report bc_multiome.report.html 
```

Repository folders:
```bash
/nextflow_version contains groovy file for nextflow pipeline calling and the environment.yml for conda setup
/original_code is a copy of the markdown files hosted on my github webpage for reference.
/src contains python, bash and R scripts called by the groovy file. These are copied to /home/groups/CEDAR/mulqueen/bc_multiome/src for running on exacloud.
/wip_src scripts that are currently being worked on to be added to the pipeline.
```
