## Run InferCNV on RNA profiles

### Batch script for InferCNV Per Sample Processing
Calling infercnv_per_sample.R script 

infercnv_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-19
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=40 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=15gb ## request gigabyte per cpu
#SBATCH --qos=long_jobs
#SBATCH --time=120:00:00 ## ask for multiple days on the node
#SBATCH --
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
seurat_obj="${proj_dir}/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"
srun Rscript ${proj_dir}/bc_multiome_nf_analysis/src/infercnv_per_sample.R \
$seurat_obj \
$SLURM_ARRAY_TASK_ID
```

## Run CopyKat on RNA profiles

### Batch script for Copykat Per Sample Processing
Calling copykat_per_sample.R script

copykat_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-19
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=10 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --partition=exacloud
#SBATCH --

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
seurat_obj="${proj_dir}/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"
srun Rscript ${proj_dir}/src/copykat_per_sample.R \
$seurat_obj \
$SLURM_ARRAY_TASK_ID
```

## Run CaSpER on RNA profiles
Calling casper_per_sample.R script written above

casper_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-19
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=10 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=20gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 1 hour on the node
#SBATCH --partition=exacloud
#SBATCH --
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
seurat_obj="${proj_dir}/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"

srun Rscript ${proj_dir}/src/casper_per_sample.R $seurat_obj \
$SLURM_ARRAY_TASK_ID \
$proj_dir

```

## CopyscAT for ATAC CNV Calling 

Using scATAC calling algorithm copyscAT from git repo https://github.com/spcdot/CopyscAT/

### Now Running samples

Code from https://github.com/spcdot/CopyscAT/blob/master/copyscat_tutorial.R
Initialize reference genome information for CopyscAT.

```R
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)

#Generate tile references
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38" ,tileWidth = 1e6,
  outputDir = "/home/groups/CEDAR/mulqueen/bc_multiome/ref/copyscat")
```

### Batch script for copyscAT Per Sample Processing

Calling copyscat_per_sample.R script written above

copyscat_slurm.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --array=1-19
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=10 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 24 hour on the node
#SBATCH --partition=exacloud
#SBATCH --

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
seurat_obj="${proj_dir}/nf_analysis/seurat_objects/merged.geneactivity.SeuratObject.rds"

srun Rscript ${proj_dir}/src/copyscat_per_sample.R \
$seurat_obj \
$SLURM_ARRAY_TASK_ID \
$proj_dir
```

### Job submit all CNV callers via slurm on exacloud.

```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
cnv_dir="${proj_dir}/nf_analysis/cnv_analysis"
#submitting each while in the appropriat wd (i'm not sure if thats relevant)
cd ${cnv_dir}/infercnv; sbatch infer_slurm.sh
cd ${cnv_dir}/casper; sbatch casper_slurm.sh
cd ${cnv_dir}/copyscat; sbatch copyscat_slurm.sh
cd ${cnv_dir}/copykat; sbatch copykat_slurm.sh
```