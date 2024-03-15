#!/usr/bin/env python
#NOTE THIS IS FROM https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Scrublet.html
#SLIGHT MODIFICATIONS TO WORK FOR MY PROCESSING
# python /src/scrublet_per_sample.py -m IDC_2/outs/filtered_feature_bc_matrix.h5 -o IDC_2/outs

import argparse
import sys
import os
import scanpy
parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
# parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix directory containing matrix files or full path to matrix.mtx. Can also also provide the 10x h5.")
# parser.add_argument("-t", "--scrublet_doublet_threshold", required = False, default = None, type = float, help = "Manually Set the scrublet doublet threshold location. For running a second time if scrublet incorrectly places the threshold the first time")
#parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory")
args = parser.parse_args()

import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
import numba
import numba.typed

# Get path of mods directory from current script directory
mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods"
sys.path.append(mods_path)
import read10x

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)


plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

## Basic run with scrublet
counts_matrix = scanpy.read_10x_h5(args.counts_matrix)
dbl_rate = counts_matrix.shape[0]/1000 * 0.008
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix.X, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)


### Plotting and saving
scrub.plot_histogram();
plt.savefig(os.path.join(args.outdir,'doublet_score_histogram.png'))
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True);
plt.savefig(os.path.join(args.outdir,'UMAP.png'))
results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
dataframe = pd.concat([barcodes_df, results, scores], axis=1)
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")
print("Writing output.\n")
dataframe.to_csv(os.path.join(args.outdir,'scrublet_results.tsv'), sep = "\t", index = False)

### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

print("Writing summary.\n")

print("Writing summary to {}.".format(os.path.join(args.outdir,'scrublet_summary.tsv')))
summary.to_csv(os.path.join(args.outdir,'scrublet_summary.tsv'), sep = "\t", index = False)

print("Done!")
