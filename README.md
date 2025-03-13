# Nextflow processing

This repository contains the nextflow pipeline and accessory scripts for processing samples from standard cellranger-arc output through final analysis for the breast cancer multiome CEDAR project. 

The file /nextflow_version/bc_multiome.nf.groovy contains the nextflow pipeline and calls R and bash scripts from src/

The code in processing_10x_data.md details the full processing steps for bcl files to cellranger-arc output. Which is subsequently used for nextflow input.

Processing of public data, used in cell label transfer, if contained in the processing_public_data.md.

The file sif_creation.md details the creation of a Singularity container using AWS (for sudo privledges) and transferring it to the appropriate location for the script to run. All nextflow processes are run within this single container.

Example run of nextflow:
```bash
#request interactive node
srun --partition=guest --time=1-12:00:00 --cpus-per-task=30 --mem=400G --nodes=1 --pty /bin/bash

#To enter SIF used for processing steps:
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

cd /home/groups/CEDAR/mulqueen/bc_multiome #move to project directory
git clone https://github.com/mulqueenr/bc_multiome_nf_analysis.git #pull github repo

#Note bed file input was originally generated from running this, 
#but it just takes like 30 hours to run the 
#merging, sorting, and peak calling, so using the already generated one.
#also copied the combined and sorted fragments file for future data publishing

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
bed_in="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/peaks/merged.nf.bed"

mkdir -p ${proj_dir}/nf_analysis_round4
cd /home/groups/CEDAR/mulqueen/bc_multiome
nextflow run bc_multiome_nf_analysis/nextflow_version/bc_multiome.nf.groovy \
--force_rewrite true \
--outdir ${proj_dir}/nf_analysis_round4 \
--sample_dir ${proj_dir}/cellranger_data/third_round \
--merged_bed ${bed_in} \
-resume
```

Repository folders:
```bash
* /nextflow_version contains groovy file for nextflow pipeline calling, the nextflow.config file for using the singularity container, and details about the creation of the container.
* /deprocated_src is a copy of the markdown files hosted on my github webpage for reference, these should be considered deprocated but they also contain archival comments and thoughts that might be helpful.
* /src contains python, bash and R scripts called by the groovy file. These are copied to /home/groups/CEDAR/mulqueen/bc_multiome/bc_multiome_nf_analysis/src for running on exacloud via the git clone command above.
* /wip_src scripts that are currently being worked on to be added to the pipeline.
```
