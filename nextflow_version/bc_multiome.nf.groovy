// Declare syntax version
//https://github.com/danrlu/Nextflow_cheatsheet/blob/main/nextflow_cheatsheet.pdf
//module load nextflow/21.10.1 (on exacloud)
nextflow.enable.dsl=2
// Script parameters
params.proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
params.outdir = "${params.proj_dir}/nf_analysis"
params.sample_dir="${params.proj_dir}/cellranger_data"
params.ref = "${params.proj_dir}/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
params.src_dir="${params.proj_dir}/src"

log.info """

		================================================
		    Breast Cancer Multiome NF PIPELINE v1.0
		================================================
		NF Working Directory : ${workflow.launchDir}
		Project Directory : ${params.proj_dir}
		Reference Directory: ${params.ref}
		Looking for samples in: ${params.sample_dir}
		================================================

""".stripIndent()


//Data correction.
process SCRUBLET_RNA {
	//Perform scrublet on raw RNA count.

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		tuple val(sample_name), path(sample_dir)
	script:
		"""
		python ${params.src_dir}/scrublet_per_sample.py ${sample_dir}
		"""

}

process SOUPX_RNA {
	//Perform soupX on raw RNA counts.

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		tuple val(sample_name), path(sample_dir)
	script:
		"""
		Rscript ${params.src_dir}/soupx_per_sample.R ${sample_name} ${sample_dir}
		"""
}


process SEURAT_GENERATION {
	//Initialize Seurat Object per sample.
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		path("${sample_dir.simpleName}.SeuratObject.rds")
	script:
		"""
		Rscript ${params.src_dir}/initialize_seurat.R ${sample_name} ${sample_dir}
		"""
}

//process MERGE_SAMPLES_CALLPEAKS {}
	//Initialize Seurat Object per sample.
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(seurat_objects)
	output:
		tuple path("merged.SeuratObject.rds"), path("combined.peakset.rds")
	script:
		"""
		Rscript ${params.src_dir}/seurat_merge_and_callpeaks.R ${seurat_objects}
		"""
}

//process DIM_REDUCTION_PER_SAMPLE {}

//process CISTOPIC_PER_SAMPLE {}

//process PUBLIC_DATA_LABEL_TRANSFER {}

//process MERGE_SEURAT_OBJECT {}

//process MERGED_CLUSTER {}

//process MERGED_CELLTYPE_BARPLOTS_AND_ALLUVIAL {}

//process MERGED_CHROMVAR {}

//process MERGED_GENE_ACTIVITY {}

// CNV CALLING
//process INFERCNV_RNA_PER_SAMPLE {}

//process CASPER_RNA_PER_SAMPLE {}

//process COPYKAT_RNA_PER_SAMPLE {}

//process COPYSCAT_ATAC_PER_SAMPLE

//CELL TYPE ANALYSIS


workflow {
	/* SETTING UP VARIABLES */
		def fasta_ref = Channel.value(params.ref)
		sample_dir = Channel
     	.fromPath("${params.sample_dir}/*/" , type: 'dir')
     	.map { [it.name, it ] }

	/* CELLRANGER TO SEURAT OBJECTS */
		seurat_objects=SCRUBLET_RNA(sample_dir) \ 
		| SOUPX_RNA \
		| SEURAT_GENERATION \
		| collect \
		| MERGE_SAMPLES_CALLPEAKS



}

