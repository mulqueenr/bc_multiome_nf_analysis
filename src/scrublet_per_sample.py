import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import sys
import pandas as pd
np.random.seed(0)
x=sys.argv[1]
input_dir=x+"/outs"
outname=x.split("/")[-1]

#Load the raw counts matrix as a scipy sparse matrix with cells as rows and genes as columns.
counts_matrix = scipy.io.mmread(input_dir + '/filtered_feature_bc_matrix/matrix.mtx.gz').T.tocsc()
cellIDs=gzip.open(input_dir + '/filtered_feature_bc_matrix/barcodes.tsv.gz',"rb").read().split()

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
#Run scrublet
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

df = pd.DataFrame({'cellid':cellIDs, 'doublet_scores':doublet_scores,'predicted_doublets':predicted_doublets})
df.to_csv(outname+'.scrublet.tsv', index=False, sep="\t")
print("Done with sample: "+outname)
print("Saved output to: "+outname+'.scrublet.tsv')