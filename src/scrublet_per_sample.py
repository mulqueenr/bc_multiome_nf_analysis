#!/usr/bin/env python
#NOTE THIS IS FROM https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Scrublet.html
#SLIGHT MODIFICATIONS TO WORK FOR MY PROCESSING
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif" 
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome --bind /home/groups/CEDAR/mulqueen/bc_multiome/bc_multiome_nf_analysis/src/:/src $sif
# python /src/scrublet_per_sample.py -m IDC_2/outs/filtered_feature_bc_matrix.h5 -o IDC_2/outs
import argparse
import sys
import os
parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = False, help = "cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory")
args = parser.parse_args()
#args.counts_matrix="IDC_2/outs/filtered_feature_bc_matrix.h5"
#args.outdir="IDC_2/outs"
import collections
import scipy.sparse as sp_sparse
import tables
import scrublet as scr
import scipy.io
import numpy as np
import pandas as pd

CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])

def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            feature_ref[key] = getattr(feature_group, key.decode()).read()
        return CountMatrix(feature_ref, barcodes, matrix)



filtered_matrix_h5 = get_matrix_from_h5(args.counts_matrix)
counts_matrix = filtered_matrix_h5[2]
counts_matrix = counts_matrix.T

## Basic run with scrublet
dbl_rate = counts_matrix.shape[0]/1000 * 0.008
print('Counts matrix shape: {} rows (cells), {} columns (features)'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
dataframe = pd.concat([barcodes_df, results, scores], axis=1)
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")
print("Writing output.\n")
dataframe.to_csv(os.path.join(args.outdir,'scrublet_results.tsv'), sep = "\t", index = False)
