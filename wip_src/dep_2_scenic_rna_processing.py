#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/scenicplus.sif
#cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects

import os
import pandas as pd
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(
    description="Function to run single cistopic model.")

parser.add_argument("-g", "--rna_genes", default = "epi_rna_counts.genes.csv", help = "Genes list csv format")
parser.add_argument("-c", "--rna_cells", default = "epi_rna_counts.cells.csv", help = "Cells list csv format")
parser.add_argument("-n", "--rna_matrix", default = "epi_rna_counts.mtx", help = "Raw counts, mtx")
parser.add_argument("-m", "--meta", default = 'epi_metadata.rna.csv', help = "Metadata table, csv")
parser.add_argument("-o", "--outDir", default ="./", help = "Output Directory")

args = parser.parse_args()

# Project directory and files
outDir = args.outDir


rna_genes =  pd.read_csv(args.rna_genes)
rna_cells =  pd.read_csv(args.rna_cells)
cell_data =  pd.read_csv(args.meta)
adata = sc.read_mtx(args.rna_matrix)

#standard QC processing
adata = adata.transpose()
adata.obs = cell_data
adata.obs_names = rna_cells['x']
adata.var_names = rna_genes['x']
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write(os.path.join(outDir, "scenicplus_"+args.rna_genes.split('.')[0]+"_rna.h5ad"))

