
Multiome processing for 10X multiome data on Primary Tumors (Phase 1+2 and preliminary experiment combined).
Analysis was performed on Exacloud, OHSU's HPC which uses slurm as a job scheduler. So many parallelized analyses utilize slurm batch processing.

I also set up my environment to automatically paste figures to a slack channel, so you may notice many system calls like "slack -F [file] slack-channel", these are just a convience function for myself. 


## Set up environment
Use ./nextflow_version/bc_conda_env.yml to set up proper conda environment.

```bash
#install conda https://conda.io/projects/conda/en/latest/user-guide/install/index.html
#conda install -c conda-forge mamba #use conda to install mamba (faster solver)
mamba env create -f nextflow_version/bc_conda_env.yml -n bc_multiome #use mamba to install environment
conda activate bc_multiome #activate
```

## Raw File Locations

Project Directory
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
```
Preliminary
```
${proj_dir}/sequencing_data/EXP211227HM #atac
${proj_dir}/sequencing_data/EXP211228HM #gex
```

Phase 1
```
${proj_dir}/sequencing_data/EXP220411HM #atac
${proj_dir}/sequencing_data/EXP220412HM #gex
```

Phase 2
```
${proj_dir}/sequencing_data/EXP220629HM #atac
${proj_dir}/sequencing_data/EXP220628HM #gex
```

## Sample Sheets
Samples are named by diagnosis and order of processing. NAT samples are matched in number to their respective IDC.

```bash
#Preliminary recieved demultiplexed data from core
echo """Lane,Sample,Index
*,ILC_1,SI-NA-A2
*,IDC_5,SI-NA-B2
*,IDC_12,SI-NA-C2
*,NAT_4,SI-NA-D2
""" > ${proj_dir}/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/multiome_atac.csv
echo """Lane,Sample,Index
*,ILC_1,SI-TT-A2
*,IDC_5,SI-TT-B2
*,IDC_12,SI-TT-C2
*,NAT_4,SI-TT-D2
""" > ${proj_dir}/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/multiome_rna.csv

#Phase 1 recieved raw BCL files from core
echo """Lane,Sample,Index
*,DCIS_1,SI-NA-A3
*,IDC_1,SI-NA-B3
*,IDC_2,SI-NA-C3
*,IDC_6,SI-NA-D3
*,IDC_7,SI-NA-E3
*,DCIS_2,SI-NA-F3
*,IDC_8,SI-NA-G3
*,IDC_9,SI-NA-H3
*,IDC_3,SI-NA-A4
*,IDC_4,SI-NA-B4
*,IDC_10,SI-NA-C4""" > ${proj_dir}/sequencing_data/EXP220411HM/220412_A01058_0226_AHCVKLDRX2-raw/multiome_atac.csv
echo """Lane,Sample,Index
*,DCIS_1,SI-TT-A5
*,IDC_1,SI-TT-B5
*,IDC_2,SI-TT-C5
*,IDC_6,SI-TT-D5
*,IDC_7,SI-TT-E5
*,DCIS_2,SI-TT-F5
*,IDC_8,SI-TT-G5
*,IDC_9,SI-TT-H5
*,IDC_3,SI-TT-A6
*,IDC_4,SI-TT-B6
*,IDC_10,SI-TT-C6""" > ${proj_dir}/sequencing_data/EXP220412HM/220413_A01058_0228_AHCVLCDRX2-raw/multiome_rna.csv


#Phase 2 recieved raw BCL files from core
echo """Lane,Sample,Index
*,NAT_11,SI-NA-B6
*,DCIS_3,SI-NA-B5
*,NAT_14,SI-NA-B2
*,IDC_11,SI-NA-B1""" > ${proj_dir}/sequencing_data/EXP220629HM/220713_A01058_0247_AHFJY3DRX2/multiome_atac.csv
echo """Lane,Sample,Index
1,NAT_11,SI-TT-G3
1,DCIS_3,SI-TT-H3
1,NAT_14,SI-TT-G4
1,IDC_11,SI-TT-H4""" > ${proj_dir}/sequencing_data/EXP220628HM/220713_A01058_0246_BHFMTNDRX2/multiome_rna.csv

```

# Initial Processing

Run Cellranger-arc mkfastq
${proj_dir}/sequencing_data/make_fq_phase_1_atac.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (ex:5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 24 hour on the node
#SBATCH --
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"

srun cellranger-arc mkfastq --id=phase_1_atac \
	--run=${proj_dir}/sequencing_data/EXP220411HM/220412_A01058_0226_AHCVKLDRX2-raw \
	--csv=${proj_dir}/sequencing_data/EXP220411HM/220412_A01058_0226_AHCVKLDRX2-raw/multiome_atac.csv \
	--localcores=20 \
	--localmem=80 \
	--output-dir=${proj_dir}/sequencing_data/fq
```
${proj_dir}/sequencing_data/make_fq_phase_1_rna.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (ex:5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 24 hour on the node
#SBATCH --
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"

srun cellranger-arc mkfastq --id=phase_1_rna \
	--run=${proj_dir}/sequencing_data/EXP220412HM/220413_A01058_0228_AHCVLCDRX2-raw \
	--csv=${proj_dir}/sequencing_data/EXP220412HM/220413_A01058_0228_AHCVLCDRX2-raw/multiome_rna.csv \
	--localcores=20 \
	--localmem=80 \
	--output-dir=${proj_dir}/sequencing_data/fq
```

${proj_dir}/sequencing_data/make_fq_phase_2_atac.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (ex:5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 24 hour on the node
#SBATCH --
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"

srun cellranger-arc mkfastq --id=phase_2_atac \
	--run=${proj_dir}/sequencing_data/EXP220629HM/220713_A01058_0247_AHFJY3DRX2 \
	--csv=${proj_dir}/sequencing_data/EXP220629HM/220713_A01058_0247_AHFJY3DRX2/multiome_atac.csv \
	--localcores=20 \
	--localmem=80 \
	--output-dir=${proj_dir}/sequencing_data/fq
```
${proj_dir}/sequencing_data/make_fq_phase_2_rna.sh
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 1 node
#SBATCH --tasks-per-node=1 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (ex:5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=5gb ## request gigabyte per cpu
#SBATCH --time=24:00:00 ## ask for 24 hour on the node
#SBATCH --
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"

srun cellranger-arc mkfastq --id=phase_2_rna \
	--run=${proj_dir}/sequencing_data/EXP220628HM/220713_A01058_0246_BHFMTNDRX2 \
	--csv=${proj_dir}/sequencing_data/EXP220628HM/220713_A01058_0246_BHFMTNDRX2/multiome_rna.csv \
	--localcores=20 \
	--localmem=80 \
	--use-bases-mask=Y28,I10,I10,Y90 \
	--lanes=1 \
	--output-dir=${proj_dir}/sequencing_data/fq
```

```bash
srun ${proj_dir}/sequencing_data/make_fq_phase_1_atac.sh
srun ${proj_dir}/sequencing_data/make_fq_phase_1_rna.sh
srun ${proj_dir}/sequencing_data/make_fq_phase_2_atac.sh
srun ${proj_dir}/sequencing_data/make_fq_phase_2_rna.sh
```
## Specify File Location
Generate libraries csv file specifying fastq locations for cellranger-arc. Note that preliminary batch of data is different format because the core split data rather than supplied bcl files.

```bash
cd ${proj_dir}/sequencing_data

#preliminary batch
prelim_atac_dir="${proj_dir}/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM"
prelim_rna_dir="${proj_dir}/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM"

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-A2_1,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-A2_2,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-A2_3,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-A2_4,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_1,Gene Expression""" > ILC_1.csv

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-B2_2,EXP211227HM_RM_2,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-B2_3,EXP211227HM_RM_2,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-B2_4,EXP211227HM_RM_2,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_2,Gene Expression""" > IDC_5.csv

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-C2_1,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-C2_2,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-C2_3,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-C2_4,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_3,Gene Expression""" > IDC_12.csv

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-D2_1,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-D2_2,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-D2_3,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-D2_4,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_4,Gene Expression""" > NAT_4.csv

#phase 1
for i in "DCIS_1" "IDC_1" "IDC_2" "IDC_6" "IDC_7" "DCIS_2" "IDC_8" "IDC_9" "IDC_3" "IDC_4" "IDC_10"; do
echo """fastqs,sample,library_type
${proj_dir}/sequencing_data/fq/HFJY3DRX2/${i},${i},Chromatin Accessibility
${proj_dir}/sequencing_data/fq/HCVLCDRX2/${i},${i},Gene Expression""" > ${i}.csv ; done

#phase 2
for i in "NAT_11" "DCIS_3" "NAT_14" "IDC_11"; do
echo """fastqs,sample,library_type
${proj_dir}/sequencing_data/fq/HCVKLDRX2/${i},${i},Chromatin Accessibility
${proj_dir}/sequencing_data/fq/HFMTNDRX2/${i},${i},Gene Expression""" > ${i}.csv ; done

```



Then use cellranger-arc count for processing.
```bash
sample_list=$(ls ${proj_dir}/sequencing_data/*.csv)
for i in $sample_list; do
  outname=${i::-4};
  cellranger-arc count --id=${outname} \
   --reference=${proj_dir}/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
   --libraries=${proj_dir}/cellranger_data/${i} \
   --localcores=30 \
   --localmem=90 
done &

```
