#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=40 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --partition=exacloud
#SBATCH --mem-per-cpu=10gb ## request gigabyte per cpu
#SBATCH --time=12:00:00 ## ask for 3 hour on the node
#SBATCH --chdir="/home/groups/CEDAR/mulqueen/bc_multiome" ## ask for 3 hour on the node
#SBATCH --

module load singularity #load singularity
module load nextflow #load nextflow
cd /home/groups/CEDAR/mulqueen/bc_multiome
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
bed="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/merged.nf.bed" #using established bed file

#run nextflow with defaults
srun --ntasks-per-core=1 \
--cpus-per-task=40 \
--mem-per-cpu=10gb \
--time=2160 \
--chdir="/home/groups/CEDAR/mulqueen/bc_multiome" \
--partition=exacloud \
--error=multiome_bc.err \
singularity shell /home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif \
nextflow run bc_multiome_nf_analysis/nextflow_version/bc_multiome.nf.groovy \
--merged_bed $bed \
-with-singularity $sif &

#-resume
#-with-dag bc_multiome.flowchart.png \
#-with-report bc_multiome.report.html \
