```bash
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/scenicplus.sif
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects
```

```python
#also ran topics #[15,20,25,30,35,40, into same tmpDir
import os
import pickle
import scanpy as sc
import pandas as pd
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import evaluate_models
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

tmpDir="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/tmp/"
outDir="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/output/"

models = []
for file in os.listdir(tmpDir):
    if file.endswith(".pkl") and file.startswith("Topic"):
        model = pickle.load(open(os.path.join(tmpDir, file), "rb"))
        models.append(model)

#output all together as model list
pickle.dump(
    models,
    open(os.path.join(outDir, "models.pkl"), "wb")
)

cistopic_models=pickle.load(open(os.path.join(outDir, "models.pkl"),"rb"))
cistopic_obj=pickle.load(open(os.path.join(outDir, "cistopic_obj.pkl"),"rb"))

#i think 40 is right
cistopic_model = evaluate_models(
    models,
    select_model = 40,
    plot_metrics = True,
    return_model = True ,
    save = "LDA_topic_select.pdf"
)


cistopic_obj.add_LDA_model(cistopic_model)
pickle.dump(
    cistopic_obj,
    open(os.path.join(outDir, "cistopic_obj.pkl"), "wb")
)

#go through this https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html
#then

#cluster sanity check
find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)
run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

run_tsne(
    cistopic_obj,
    target  = 'cell', scale=True)

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['assigned_celltype', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=4,
    text_size=10,
    dot_size=5,
    save="cistopic_umap.pdf")

#heatmap of cell type per topic
cell_topic_heatmap(
    cistopic_obj,
    variables = ['assigned_celltype'],
    scale = False,
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 10),
    save="cistopic_celltype_heatmap.pdf"
)


#heatmap of cell type per topic
cell_topic_heatmap(
    cistopic_obj,
    variables = ['Diag_MolDiag'],
    scale = False,
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 10),
    save="cistopic_diagmoldiat_heatmap.pdf"
)

#binarize topics 
#by region
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)
#by cell
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)

#from this, determine topic metrics (like specificity)
topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1

plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.savefig('topic_metrics.pdf')

#now annotate by cell type and by diagnoses
topic_annot_celltype = topic_annotation(
    cistopic_obj,
    annot_var='assigned_celltype',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)

topic_annot_celltype.to_csv("topic_specificity_celltype.csv")

topic_annot_diagmoldiag = topic_annotation(
    cistopic_obj,
    annot_var='Diag_MolDiag',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)
topic_annot_diagmoldiag.to_csv("topic_specificity_diag_moldiag.csv")


#find DARs
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np

imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

normalized_imputed_acc_obj = normalize_scores(
    imputed_acc_obj, 
    scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)

markers_dict_celltype = find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='assigned_celltype',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir="/home/users/mulqueen/tmp",
    split_pattern = '-'
)

#note this is across all cells, so probably not the best, should subset to cancer
markers_dict_diag_moldiag = find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='Diag_MolDiag',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir="/home/users/mulqueen/tmp",
    split_pattern = '-'
)

#save final cistopic obj
pickle.dump(
    cistopic_obj,
    open(os.path.join(outDir, "cistopic_obj.pkl"), "wb")
)

from pycisTopic.utils import region_names_to_coordinates
import re

#save bed files
os.makedirs(os.path.join(outDir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "DARs_cell_type"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "DARs_Diag_MolDiag"), exist_ok = True)

#topic lists
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

#cell type DAR
for cell_type in markers_dict_celltype:
    region_names_to_coordinates(
        markers_dict_celltype[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

#Diag MolDiag DAR
for Diag_MolDiag in markers_dict_diag_moldiag:
    regions=region_names_to_coordinates(markers_dict_diag_moldiag[Diag_MolDiag].index).sort_values(["Chromosome", "Start", "End"])
    outname=re.sub('\+', 'plus', Diag_MolDiag)
    outname=re.sub("\-","minus",outname)
    outname=re.sub("\/","_",outname)
    outname=re.sub("\ ","_",outname)
    regions.to_csv(
        os.path.join(outDir, "region_sets", "DARs_Diag_MolDiag", f"{outname}.bed"),
        sep = "\t",
        header = False, index = False
    )

#removed tag in initial cistopic generation too
cistopic_obj.cell_names=[re.sub("___cisTopic","",x) for x in cistopic_obj.cell_names]
#save final cistopic obj
pickle.dump(
    cistopic_obj,
    open(os.path.join(outDir, "cistopic_obj.pkl"), "wb")
)


```

```bash
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/output

mkdir -p scplus_pipeline

scenicplus init_snakemake \
--out_dir scplus_pipeline

mkdir -p outs
mkdir -p tmp

cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/output


cistopic_obj="$PWD/cistopic_obj.pkl"
GEX_anndata="$PWD/scenicplus_rna.h5ad"
region_set_folder="$PWD/region_sets"
ctx_db="$PWD/bc_multiome.regions_vs_motifs.rankings.feather"
dem_db="$PWD/bc_multiome.motifs_vs_regions.scores.feather"
path_to_motif_ann="/home/groups/CEDAR/mulqueen/bc_multiome/ref/aertslab_motif_collection/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl"
#configure yaml
#using # as delim
sed -i "s#cisTopic_obj_fname: \"\"#cisTopic_obj_fname: \"$cistopic_obj\"#" scplus_pipeline/Snakemake/config/config.yaml

sed -i "s#GEX_anndata_fname: \"\"#GEX_anndata_fname: \"$GEX_anndata\"#" scplus_pipeline/Snakemake/config/config.yaml

sed -i "s#region_set_folder: \"\"#region_set_folder: \"$region_set_folder\"#" scplus_pipeline/Snakemake/config/config.yaml

sed -i "s#ctx_db_fname: \"\"#ctx_db_fname: \"$ctx_db\"#" scplus_pipeline/Snakemake/config/config.yaml

sed -i "s#dem_db_fname: \"\"#dem_db_fname: \"$dem_db\"#" scplus_pipeline/Snakemake/config/config.yaml

sed -i "s#path_to_motif_annotations: \"\"#path_to_motif_annotations: \"$path_to_motif_ann\"#" scplus_pipeline/Snakemake/config/config.yaml

cd scplus_pipeline/Snakemake/
snakemake --cores 20


#additional functions to maybe prerun?
scenicplus prepare_data download_genome_annotations \
--species hsapiens \
--biomart_host http://www.ensembl.org \
--genome_annotation_out_fname /home/jupyter-ayang/notebooks/scenic_plus/pbmc_development/output/genome_annotation.tsv\
--chromsizes_out_fname /home/jupyter-ayang/notebooks/scenic_plus/pbmc_development/output/chromsizes.tsv

scenicplus prepare_data prepare_GEX_ACC \
--cisTopic_obj_fname $cistopic_obj \
--GEX_anndata_fname $GEX_anndata \
--out_file "${PWD}/ACC_GEX.h5mu" \
--bc_transform_func "lambda x: f'{x}'"




#add output files genome_annotation.tsv and 
```

