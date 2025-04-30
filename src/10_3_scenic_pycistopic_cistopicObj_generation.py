#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/scenicplus.sif
#cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects
import os
from pycisTopic.cistopic_class import *
from scipy.io import mmread
import pickle

parser.add_argument("-f", "--frag_path", default = 'frag_paths.csv', help = "List of fragment locations, csv format")
parser.add_argument("-n", "--atac_counts", default = 'atac_counts.mtx', help = "Raw counts, mtx")
parser.add_argument("-c", "--atac_cells", default = 'atac_counts.cells.csv', help = "List of cells, csv")
parser.add_argument("-p", "--atac_peaks", default = 'atac_counts.peaks.csv', help = "List of peaks, csv")
parser.add_argument("-m", "--meta", default = 'metadata.csv', help = "Metadata table, csv")
parser.add_argument("-o", "--outDir", default ="./", help = "Output Directory")

# Project directory and files
frag_path=args.frag_path
meta=args.meta
atac_counts=args.atac_counts
atac_cells=args.atac_cells
atac_peaks=args.atac_peaks

# Output directory
if not os.path.exists(outDir):
    os.makedirs(outDir)

# Create cisTopic object
atac_counts = mmread(atac_counts)
atac_peaks =  pd.read_csv(atac_peaks)
atac_peaks = [peak.replace("-",":",1) for peak in atac_peaks['x']]  #reformat name
atac_cells =  pd.read_csv(atac_cells)
cell_data =  pd.read_csv(meta,dtype="string")

cistopic_obj = create_cistopic_object(fragment_matrix=atac_counts.tocsr(),
    cell_names=atac_cells['x'].tolist(),
    region_names=atac_peaks,
    tag_cells=False)

# Adding cell information
cistopic_obj.add_cell_data(cell_data)
pickle.dump(
    cistopic_obj,
    open(os.path.join(outDir, "scenicplus"+args.atac_counts.split('_')[0]+"_cistopic_obj.pkl"), "wb")
)
