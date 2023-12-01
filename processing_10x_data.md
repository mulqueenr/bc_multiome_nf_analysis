
# Initial cellranger run of data
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

#Phase 2 recieved raw BCL files from core
echo """Lane,Sample,Index
*,NAT_11,SI-NA-B6
*,DCIS_3,SI-NA-B5
*,NAT_14,SI-NA-B2
*,IDC_11,SI-NA-B1""" > ${proj_dir}/sequencing_data/EXP220629HM/220713_A01058_0247_AHFJY3DRX2/multiome_atac.csv


```
## Generation of RNA Libraries bcl2fastq
Prelim 1 rna workflow b
https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/bcl2fastq-direct
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
run_dir="${proj_dir}/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY"
echo """
[Header]
EMFileVersion,4
 
[Reads]
28
90
 
[Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,Original_Sample_ID
,ILC_1,ILC_1,GTGGATCAAA,CAGGGTTGGC,preliminary_1_rna,SI-TT-A2
,IDC_5,IDC_5,TCTACCATTT,GACTCTCCCG,preliminary_1_rna,SI-TT-B2
,IDC_12,IDC_12,CAATCCCGAC,TACTACTCGG,preliminary_1_rna,SI-TT-C2
,NAT_4,NAT_4,TTAATACGCG,ACCCGAGGTG,preliminary_1_rna,SI-TT-D2""" > ${run_dir}/prelim_1_rna.samplesheet.csv

cd $run_dir
bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --sample-sheet=${run_dir}/prelim_1_rna.samplesheet.csv
```


Phase 1 rna workflow b 
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
run_dir="${proj_dir}/sequencing_data/EXP220412HM/220413_A01058_0228_AHCVLCDRX2-raw"
echo """
[Header]
EMFileVersion,4
 
[Reads]
28
90
 
[Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,Original_Sample_ID
,DCIS_1,DCIS_1,GTAGCCCTGT,ATAGATGCTC,phase_1_rna,SI-TT-A5
,IDC_1,IDC_1,TCGGCTCTAC,AGACCATCGG,phase_1_rna,SI-TT-B5
,IDC_2,IDC_2,TCCGTTGGAT,GCGAGAACGT,phase_1_rna,SI-TT-C5
,IDC_6,IDC_6,TGGTTCGGGT,CTCCTGCCAC,phase_1_rna,SI-TT-D5
,IDC_7,IDC_7,CGCGGTAGGT,CAACATCCTG,phase_1_rna,SI-TT-E5
,DCIS_2,DCIS_2,CGGCTGGATG,GTGCTTATCA,phase_1_rna,SI-TT-F5
,IDC_8,IDC_8,ATAGGGCGAG,ACTCGATGCA,phase_1_rna,SI-TT-G5
,IDC_9,IDC_9,AGCAAGAAGC,AGAAACACAA,phase_1_rna,SI-TT-H5
,IDC_3,IDC_3,TAACGCGTGA,GAAGTTAGGG,phase_1_rna,SI-TT-A6
,IDC_4,IDC_4,AATGCCATGA,GCATTACGTA,phase_1_rna,SI-TT-B6
,IDC_10,IDC_10,ACGACTACCA,TTAGGGTCGT,phase_1_rna,SI-TT-C6""" > ${run_dir}/phase_1_rna.samplesheet.csv


cd $run_dir
bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --sample-sheet=${run_dir}/phase_1_rna.samplesheet.csv
```
Phase 2 rna workflow b 
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
run_dir="${proj_dir}/sequencing_data/EXP220628HM/220713_A01058_0246_BHFMTNDRX2"
echo """
[Header]
EMFileVersion,4
 
[Reads]
28
90
 
[Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,Original_Sample_ID
1,NAT_11,NAT_11,ATGACGTCGC,ATCCTGACCT,phase_2_rna,SI-TT-G3
1,DCIS_3,DCIS_3,CCCGTTCTCG,CCAATCCGTC,phase_2_rna,SI-TT-H3
1,NAT_14,NAT_14,GCGCTTATGG,CTAGCCAGGC,phase_2_rna,SI-TT-G4
1,IDC_11,IDC_11,AGTTTCCTGG,CTGTGTGGCA,phase_2_rna,SI-TT-H4""" > ${run_dir}/phase_2_rna.samplesheet.csv

cd $run_dir
bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --sample-sheet=${run_dir}/phase_2_rna.samplesheet.csv
```

## Generation of ATAC libraries bcl2fastq
Preliminary 1 atac data
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
run_dir="${proj_dir}/sequencing_data/EXP220412HM/211229_A01058_0202_AHT33MDRXY"
echo """
[Header]
EMFileVersion,4

[Reads]
150
150

[Data]
Lane,Sample_ID,Sample_Name,index,Sample_Project
,SI-NA-A2_1,ILC_1,AGCCCTTT,preliminary_1_atac
,SI-NA-A2_2,ILC_1,CAAGTCCA,preliminary_1_atac
,SI-NA-A2_3,ILC_1,GTGAGAAG,preliminary_1_atac
,SI-NA-A2_4,ILC_1,TCTTAGGC,preliminary_1_atac
,SI-NA-B2_1,IDC_5,AAGTTGAT,preliminary_1_atac
,SI-NA-B2_2,IDC_5,CCCACCCA,preliminary_1_atac
,SI-NA-B2_3,IDC_5,GGTCGAGC,preliminary_1_atac
,SI-NA-B2_4,IDC_5,TTAGATTG,preliminary_1_atac
,SI-NA-C2_1,IDC_12,AATCACTA,preliminary_1_atac
,SI-NA-C2_2,IDC_12,CCGAGAAC,preliminary_1_atac
,SI-NA-C2_3,IDC_12,GTAGTGCG,preliminary_1_atac
,SI-NA-C2_4,IDC_12,TGCTCTGT,preliminary_1_atac
,SI-NA-D2_1,NAT_4,ACATTCCG,preliminary_1_atac
,SI-NA-D2_2,NAT_4,CTGCGGTA,preliminary_1_atac
,SI-NA-D2_3,NAT_4,GACACAAT,preliminary_1_atac
,SI-NA-D2_4,NAT_4,TGTGATGC,preliminary_1_atac""" > ${run_dir}/preliminary_1_atac.samplesheet.csv


cd $run_dir
bcl2fastq \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --sample-sheet=${run_dir}/preliminary_1_atac.samplesheet.csv
```

Phase 1 atac data
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
run_dir="${proj_dir}/sequencing_data/EXP220411HM/220412_A01058_0226_AHCVKLDRX2-raw"
echo """
[Header]
EMFileVersion,4

[Reads]
150
150

[Data]
Lane,Sample_ID,Sample_Name,index,Sample_Project
,SI-NA-A3_1,DCIS_1,AAAGCATA,phase_1_atac
,SI-NA-A3_2,DCIS_1,CTGCAGCC,phase_1_atac
,SI-NA-A3_3,DCIS_1,GCCTTTAT,phase_1_atac
,SI-NA-A3_4,DCIS_1,TGTAGCGG,phase_1_atac
,SI-NA-B3_1,IDC_1,ATTGGACG,phase_1_atac
,SI-NA-B3_2,IDC_1,CAGCTTAC,phase_1_atac
,SI-NA-B3_3,IDC_1,GGCAAGGA,phase_1_atac
,SI-NA-B3_4,IDC_1,TCATCCTT,phase_1_atac
,SI-NA-C3_1,IDC_2,ACGTTACA,phase_1_atac
,SI-NA-C3_2,IDC_2,CGTAGGTT,phase_1_atac
,SI-NA-C3_3,IDC_2,GACGACGG,phase_1_atac
,SI-NA-C3_4,IDC_2,TTACCTAC,phase_1_atac
,SI-NA-D3_1,IDC_6,ACTTCACT,phase_1_atac
,SI-NA-D3_2,IDC_6,CGAAGTTG,phase_1_atac
,SI-NA-D3_3,IDC_6,GAGCACGC,phase_1_atac
,SI-NA-D3_4,IDC_6,TTCGTGAA,phase_1_atac
,SI-NA-E3_1,IDC_7,AACAAGTC,phase_1_atac
,SI-NA-E3_2,IDC_7,CGGCTCCA,phase_1_atac
,SI-NA-E3_3,IDC_7,GTATGTAT,phase_1_atac
,SI-NA-E3_4,IDC_7,TCTGCAGG,phase_1_atac
,SI-NA-F3_1,DCIS_2,AGTCTGTA,phase_1_atac
,SI-NA-F3_2,DCIS_2,CAGAATAG,phase_1_atac
,SI-NA-F3_3,DCIS_2,GCCTCCGT,phase_1_atac
,SI-NA-F3_4,DCIS_2,TTAGGACC,phase_1_atac
,SI-NA-G3_1,IDC_8,ATGTCCAG,phase_1_atac
,SI-NA-G3_2,IDC_8,CGACGTCA,phase_1_atac
,SI-NA-G3_3,IDC_8,GCTATAGC,phase_1_atac
,SI-NA-G3_4,IDC_8,TACGAGTT,phase_1_atac
,SI-NA-H3_1,IDC_9,ACACCTAA,phase_1_atac
,SI-NA-H3_2,IDC_9,CGTTTGGG,phase_1_atac
,SI-NA-H3_3,IDC_9,GACAAACC,phase_1_atac
,SI-NA-H3_4,IDC_9,TTGGGCTT,phase_1_atac
,SI-NA-A4_1,IDC_3,AGAACGCC,phase_1_atac
,SI-NA-A4_2,IDC_3,CATGGCAG,phase_1_atac
,SI-NA-A4_3,IDC_3,GTCTTTGA,phase_1_atac
,SI-NA-A4_4,IDC_3,TCGCAATT,phase_1_atac
,SI-NA-B4_1,IDC_4,AGGGACTG,phase_1_atac
,SI-NA-B4_2,IDC_4,CCTCTAAC,phase_1_atac
,SI-NA-B4_3,IDC_4,GACAGGCT,phase_1_atac
,SI-NA-B4_4,IDC_4,TTATCTGA,phase_1_atac
,SI-NA-C4_1,IDC_10,ACATTGGC,phase_1_atac
,SI-NA-C4_2,IDC_10,CTTAGTCA,phase_1_atac
,SI-NA-C4_3,IDC_10,GAGCCCAT,phase_1_atac
,SI-NA-C4_4,IDC_10,TGCGAATG,phase_1_atac""" > ${run_dir}/phase_1_atac.samplesheet.csv


cd $run_dir
bcl2fastq \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --sample-sheet=${run_dir}/phase_1_atac.samplesheet.csv
```

Phase 2 atac data
```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
run_dir="${proj_dir}/sequencing_data/EXP220629HM/220713_A01058_0247_AHFJY3DRX"

echo """
[Header]
EMFileVersion,4

[Reads]
150
150

[Data]
Lane,Sample_ID,Sample_Name,index,Sample_Project
,SI-NA-B6_1,NAT_11,AACGCGAA,phase_2_atac
,SI-NA-B6_2,NAT_11,CTATTTGG,phase_2_atac
,SI-NA-B6_3,NAT_11,GCGCACCT,phase_2_atac
,SI-NA-B6_4,NAT_11,TGTAGATC,phase_2_atac
,SI-NA-B5_1,DCIS_3,ATCGTACT,phase_2_atac
,SI-NA-B5_2,DCIS_3,CATCAGTG,phase_2_atac
,SI-NA-B5_3,DCIS_3,GGGACTAC,phase_2_atac
,SI-NA-B5_4,DCIS_3,TCATGCGA,phase_2_atac
,SI-NA-B2_1,NAT_14,AAGTTGAT,phase_2_atac
,SI-NA-B2_2,NAT_14,CCCACCCA,phase_2_atac
,SI-NA-B2_3,NAT_14,GGTCGAGC,phase_2_atac
,SI-NA-B2_4,NAT_14,TTAGATTG,phase_2_atac
,SI-NA-B1_1,IDC_11,AGGCTACC,phase_2_atac
,SI-NA-B1_2,IDC_11,CTAGCTGT,phase_2_atac
,SI-NA-B1_3,IDC_11,GCCAACAA,phase_2_atac
,SI-NA-B1_4,IDC_11,TATTGGTG,phase_2_atac""" > ${run_dir}/phase_2_atac.samplesheet.csv


cd $run_dir
bcl2fastq \
		  --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ${run_dir} \
            --sample-sheet=${run_dir}/phase_2_atac.samplesheet.csv
```

# Initial Processing

## Specify File Location
Generate libraries csv file specifying fastq locations for cellranger-arc. Note that preliminary batch of data is different format because the core split data rather than supplied bcl files.

```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
cd ${proj_dir}/sequencing_data

#preliminary batch
prelim_atac_dir="${proj_dir}/sequencing_data/EXP211227HM/211229_A01058_0202_AHT33MDRXY/EXP211227HM"
prelim_rna_dir="${proj_dir}/sequencing_data/EXP211228HM/211229_A01058_0201_BHT55FDRXY/EXP211228HM"

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-A2_1,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-A2_2,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-A2_3,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-A2_4,EXP211227HM_RM_1,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_1,Gene Expression""" > ${proj_dir}/sequencing_data/ILC_1.csv

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-B2_2,EXP211227HM_RM_2,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-B2_3,EXP211227HM_RM_2,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-B2_4,EXP211227HM_RM_2,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_2,Gene Expression""" > ${proj_dir}/sequencing_data/IDC_5.csv

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-C2_1,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-C2_2,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-C2_3,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-C2_4,EXP211227HM_RM_3,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_3,Gene Expression""" > ${proj_dir}/sequencing_data/IDC_12.csv

echo """fastqs,sample,library_type
${prelim_atac_dir}/SI-NA-D2_1,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-D2_2,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-D2_3,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_atac_dir}/SI-NA-D2_4,EXP211227HM_RM_4,Chromatin Accessibility
${prelim_rna_dir},EXP211228HM_RM_4,Gene Expression""" > ${proj_dir}/sequencing_data/NAT_4.csv

#phase 1
for i in "DCIS_1" "IDC_1" "IDC_2" "IDC_6" "IDC_7" "DCIS_2" "IDC_8" "IDC_9" "IDC_3" "IDC_4" "IDC_10"; do
echo """fastqs,sample,library_type
${proj_dir}/sequencing_data/fq/HFJY3DRX2/${i},${i},Chromatin Accessibility
${proj_dir}/sequencing_data/fq/HCVLCDRX2/${i},${i},Gene Expression""" > ${proj_dir}/sequencing_data/${i}.csv; done

#phase 2
phase2_atac_dir="${proj_dir}/sequencing_data/EXP220628HM/220713_A01058_0246_BHFMTNDRX2/Data/Intensities/BaseCalls/phase_2_atac"
phase2_rna_dir="${proj_dir}/sequencing_data/EXP220629HM/220713_A01058_0247_AHFJY3DRX/Data/Intensities/BaseCalls/phase_2_rna"
for i in "NAT_11" "DCIS_3" "NAT_14" "IDC_11"; do
echo """fastqs,sample,library_type
${phase2_atac_dir}/sequencing_data/fq/HCVKLDRX2/${i},${i},Chromatin Accessibility
${phase2_rna_dir}/sequencing_data/fq/HFMTNDRX2/${i},${i},Gene Expression""" > ${proj_dir}/sequencing_data/${i}.csv ; done

```


## Deeper sequencing performed for all samples, final cellranger run for manuscript.

Sample processing followed same pipeline as before. Rerunning cellranger to keep sample names consistent.
Initial bcl2fastq processing performed by Hugo Cros and the OHSU Sequencing Core.
New data moved to /home/groups/CEDAR/mulqueen/bc_multiome/sequencing_data/second_round_sequencing_combined
New cellranger output to /home/groups/CEDAR/mulqueen/bc_multiome/cellranger_data/second_round

```bash
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
old_dir="/home/groups/CEDAR/cros/projects/Ryan_multiome_newdata/data"
new_dir="/home/groups/CEDAR/mulqueen/bc_multiome/sequencing_data/second_round_sequencing_combined"
sample_list=($(ls ${proj_dir}/sequencing_data/second_round_sequencing_combined/["s|R"]*.csv))
list_length=$(echo ${#sample_list[@]})

#hardcoded rename of samples, sorted alphabetically by original name
#sample 15 (NAT 11) is excluded for reseq, so just added from first round
sample_names=("ILC_1" "IDC_5" "IDC_12" "NAT_4" "IDC_3" "IDC_4" "IDC_10" "DCIS_3" "NAT_14" "DCIS_1" "IDC_11" "IDC_1" "IDC_2" "IDC_6" "IDC_7" "DCIS_2" "IDC_8" "IDC_9")
names_length=$(echo ${#sample_names[@]})
echo $list_length
echo $names_length
#rewrite csv files to current directory after moved and rename output
for i in $(eval echo {0..$list_length}); do
  outname=${sample_names[i]}
  csv_in=${sample_list[i]}
  cat $csv_in | sed 's|/home/groups/CEDAR/cros/projects/Ryan_multiome_newdata/data|/home/groups/CEDAR/mulqueen/bc_multiome/sequencing_data/second_round_sequencing_combined|g' > ${outname}.csv
done

```
Then use cellranger-arc count for processing.
Using slurm submission for multinodes

cellranger_processing.slurm
```bash
#!/bin/bash
#SBATCH --nodes=1 #request 19 nodes
#SBATCH --array=1-19 #19 samples from csv files
#SBATCH --tasks-per-node=2 ##we want our node to do N tasks at the same time
#SBATCH --cpus-per-task=20 ##ask for CPUs per task (5 * 8 = 40 total requested CPUs)
#SBATCH --mem-per-cpu=3gb ## request gigabyte per cpu
#SBATCH --time=10:00:00 ## ask for 5 hour on the node
#SBATCH --partition="exacloud"
#SBATCH --

proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
file_in=$(ls ${proj_dir}/sequencing_data/second_round_sequencing_combined/*.csv | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}')

outname=${file_in::-4};
prefix=`basename $outname`
cellranger-arc count --id=${prefix} \
 --reference=${proj_dir}/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
 --libraries=${file_in} \
 --localcores=20 \
 --localmem=60 
```
```bash
sbatch cellranger_processing.slurm
```

error: FASTQ path doesn't exist: "/home/groups/CEDAR/mulqueen/bc_multiome/sequencing_data/first_round_sequencing/EXP220628HM/220713_A01058_0246_BHFMTNDRX2/Data/Intensities/BaseCalls/phase_2_rna/"
