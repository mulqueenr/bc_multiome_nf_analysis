
Download cistarget data bases.

```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/src/
git clone https://github.com/aertslab/create_cisTarget_databases
wget https://resources.aertslab.org/cistarget/programs/cbust
chmod a+x cbust

#prepare cistarget database
cd /home/groups/CEDAR/mulqueen/bc_multiome/ref/
mkdir -p aertslab_motif_collection
wget -O aertslab_motif_collection/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
cd aertslab_motif_collection; unzip -q v10nr_clust_public.zip
```



#convert output peaks to bed format
#run this is SIF with bedtools
cd /home/groups/CEDAR/mulqueen/bc_multiome/ref/aertslab_motif_collection

awk 'OFS="\t" {split($1,a,",");split(a[2],b,"-"); print b[1],b[2],b[3]}' /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/atac_counts.peaks.csv | tail -n +2 | tr -d '"' | awk 'OFS="\t" {print $1,$2,$3}' > /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/atac_peaks.bed

REGION_BED="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/atac_peaks.bed"
GENOME_FASTA="/home/groups/CEDAR/mulqueen/bc_multiome/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
CHROMSIZES="/home/groups/CEDAR/mulqueen/bc_multiome/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/chrNameLength.txt"
DATABASE_PREFIX="bc_multiome"
SCRIPT_DIR="/home/groups/CEDAR/mulqueen/bc_multiome/src/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${DATABASE_PREFIX}.fa \
        1000 \
        yes

```

```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects
```

```bash
#just run with jaspar motifs
ls /home/groups/CEDAR/mulqueen/bc_multiome/ref/aertslab_motif_collection/v10nr_clust_public/singletons/jaspar* > motifs.txt

OUT_DIR=""${PWD}""
CBDIR="/home/groups/CEDAR/mulqueen/bc_multiome/ref/aertslab_motif_collection/v10nr_clust_public/singletons"
DATABASE_PREFIX="bc_multiome"
FASTA_FILE="${OUT_DIR}/${DATABASE_PREFIX}.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
SCRIPT_DIR="/home/groups/CEDAR/mulqueen/bc_multiome/src/create_cisTarget_databases"

export PATH="$PATH:/home/groups/CEDAR/mulqueen/bc_multiome/src/create_cisTarget_databases" #cbust is in here and needs to be on PATH

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20
```