#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=40 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --partition=exacloud
#SBATCH --mem-per-cpu=12gb ## request gigabyte per cpu
#SBATCH --time=36:00:00 ## ask for 3 hour on the node
#SBATCH --chdir="/home/groups/CEDAR/mulqueen/bc_multiome" ## ask for 3 hour on the node
#SBATCH --


module load singularity/3.8.0 #load singularity
module load nextflow/21.10.1 #load nextflow
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
bed="/home/groups/CEDAR/mulqueen/bc_multiome/merged.nf.bed"
#run nextflow with defaults

nextflow run \
bc_multiome_nf_analysis/nextflow_version/bc_multiome.nf.groovy \
--merged_bed ${proj_dir}/merged.nf.bed \
--outdir ${proj_dir}/nf_analysis_round3 \
--sample_dir ${proj_dir}/cellranger_data/third_round \
-resume

