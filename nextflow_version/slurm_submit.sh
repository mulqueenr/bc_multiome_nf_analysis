#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=40 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --partition=cedar
#SBATCH --mem-per-cpu=12gb ## request gigabyte per cpu
#SBATCH --time=36:00:00 ## ask for 3 hour on the node
#SBATCH --chdir="/home/groups/CEDAR/mulqueen/bc_multiome" ## ask for 3 hour on the node
#SBATCH --

module load singularity #load singularity
module load nextflow #load nextflow
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
#run nextflow with defaults

nextflow run \
bc_multiome_nf_analysis/nextflow_version/bc_multiome.nf.groovy \
--outdir ${proj_dir}/nf_analysis_round4 \
--sample_dir ${proj_dir}/cellranger_data/third_round \
-resume

