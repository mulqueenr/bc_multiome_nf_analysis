#container
"""
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome /home/groups/CEDAR/mulqueen/bc_multiome/scenicplus.sif
cd /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects
"""
#example run:
"""
python pycistopic_run_topic.py \
--topicCount 5 \
--cistopicObj epi_cistopic_obj.pkl \
--outDir /home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/output_epi/
"""

import os
from pycisTopic.cistopic_class import *
import scanpy as sc
import pickle
from pycisTopic.lda_models import run_cgs_models_mallet
import argparse


parser = argparse.ArgumentParser(
    description="Function to run single cistopic model.")

parser.add_argument("-c", "--topicCount", required = False, default = 5, help = "Integer of topic count to be used in run_cgs_models_mallet.",dtype=int)
parser.add_argument("-d", "--outDir", required = False, default = '/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round4/seurat_objects/output_epi/', help = "The project directory directory containing the cistopic object.")
parser.add_argument("-o", "--cistopicObj", required = False, default = "epi_cistopic_obj.pkl", help = "Cistopic object pkl.")
parser.add_argument("-m", "--memory", required = False, default = '750G', help = "Memory for MALLET_MEMORY")
parser.add_argument("-t", "--taskCpus", required = False, default = 10, help = "CPUS to run",dtype=int)
parser.add_argument("-M", "--mallet", required = False, default ="/container_mallet/bin/mallet", help = "Mallet path, built into container")

args = parser.parse_args()

# Project directory and files
os.chdir(args.outDir)

#tmp dir
tmpDir = projDir + '/tmp_epi/'
if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)


cistopic_obj=pickle.load(open(os.path.join(args.outDir, args.cistopicObj), "rb"))

# Run models with mallet
os.environ['MALLET_MEMORY'] = args.memory
mallet_path= args.mallet

# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[args.topicCount], 
    n_cpu=args.taskCpus,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=tmpDir,
    save_path=tmpDir,
    mallet_path=mallet_path,
)