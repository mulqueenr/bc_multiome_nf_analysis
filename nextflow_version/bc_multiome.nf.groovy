// Declare syntax version
//module load nextflow/21.10.1 (on exacloud)
nextflow.enable.dsl=2

  //////////////////////////////
 ///	Script Parameters	///
//////////////////////////////
params.proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
params.outdir = "${params.proj_dir}/nf_analysis"
params.sample_dir="${params.proj_dir}/cellranger_data"
params.ref = "${params.proj_dir}/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
params.src_dir="${params.proj_dir}/src"
params.force_rewrite="false"

log.info """

		================================================
		    Breast Cancer Multiome NF PIPELINE v1.0
		================================================
		NF Working Directory : ${workflow.launchDir}
		Project Directory : ${params.proj_dir}
		Reference Directory: ${params.ref}
		Looking for samples in: ${params.sample_dir}
		Force Rewriting Data: ${params.force_rewrite}
		================================================

""".stripIndent()

  //////////////////////////
 ///	Data correction	///
//////////////////////////
process SCRUBLET_RNA {
	//Perform scrublet on raw RNA count.

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		tuple val(sample_name), path(sample_dir)

	script:
		"""
		if (test -e ${sample_dir}/outs/${sample_name}.scrublet.tsv) && (! ${params.force_rewrite}); then
		echo "Using stored scrublet output: ${sample_dir}/outs/${sample_name}.scrublet.tsv"
		else
		python ${params.src_dir}/scrublet_per_sample.py \\
		${sample_dir}
		fi
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
		if (test -e ${sample_dir}/outs/soupx_corrected_counts.rds) && (! ${params.force_rewrite}); then
		echo "Using stored scrublet output: ${sample_dir}/outs/soupx_corrected_counts.rds"
		else
		Rscript ${params.src_dir}/soupx_per_sample.R \\
		${sample_name} \\
		${sample_dir}
		fi
		"""
}

  //////////////////////////////////////
 ///	Seurat Sample Processing	///
//////////////////////////////////////
process SEURAT_GENERATION {
	//Initialize Seurat Object per sample.
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		path("${sample_dir.simpleName}.SeuratObject.rds")

	script:
		"""
		mkdir -p ${params.outdir}/plots
		Rscript ${params.src_dir}/initialize_seurat.R \\
		${sample_name} \\
		${sample_dir} \\
		${params.outdir}/plots
		"""
}

process MERGE_SAMPLES_CALLPEAKS {
	//Initialize Seurat Object per sample.
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(seurat_objects)
	output:
		tuple path("merged.SeuratObject.rds"), path("combined.peakset.rds")

	script:
		"""
		Rscript ${params.src_dir}/seurat_merge_and_callpeaks.R \\
		${seurat_objects}
		"""
}

process DIM_REDUCTION_PER_SAMPLE {
	//Reanalyze ATAC data with combined peak set and perform dim reduction.

	input:
		path(obj_in)
		tuple path(merged_object), path(combined_peaks)
	output:
		path("${obj_in}")
	script:
	"""
	${params.src_dir}/seurat_dim_reduction_per_sample.R \\
	${combined_peaks} \\
	${obj_in} \\
	${params.outdir}/plots 
	"""
}

process CISTOPIC_PER_SAMPLE {
	//Run cisTopic on sample ATAC data
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(obj_in)
	output:
		tuple path("${obj_in}"),path("*.CisTopicObject.rds")
	script:
	"""
	${params.src_dir}/seurat_cistopic_per_sample.R \\
	${obj_in} \\
	${params.outdir}/plots
	"""
}

//process PUBLIC_DATA_LABEL_TRANSFER {}
//${params.src_dir}/seurat_public_data_label_transfer.R


  //////////////////////////////////////////////////////
 ///	Sample Integration and Merged Processing	///
//////////////////////////////////////////////////////
//process MERGE_SEURAT_OBJECT {}

//process MERGED_CLUSTER {}

//process MERGED_CELLTYPE_BARPLOTS_AND_ALLUVIAL {}

//process MERGED_CHROMVAR {}

//process MERGED_GENE_ACTIVITY {}

  //////////////////////////////////
 ///	CNV Calling Per Sample	///
//////////////////////////////////
//process INFERCNV_RNA_PER_SAMPLE {}

//process CASPER_RNA_PER_SAMPLE {}

//process COPYKAT_RNA_PER_SAMPLE {}

//process COPYSCAT_ATAC_PER_SAMPLE

  //////////////////////////////
 ///	Cell Type Analysis	///
//////////////////////////////


workflow {
	/* SETTING UP VARIABLES */
		def fasta_ref = Channel.value(params.ref)
		sample_dir = Channel
     	.fromPath("${params.sample_dir}/*/" , type: 'dir')
     	.map { [it.name, it ] }

	/* Data correction */
		sample_seurat_objects=SCRUBLET_RNA(sample_dir) \
		| SOUPX_RNA \
		| SEURAT_GENERATION
		
	/* Seurat Sample Processing */
		merged_seurat_object =
		collect(sample_seurat_objects) \
		| MERGE_SAMPLES_CALLPEAKS

		DIM_REDUCTION_PER_SAMPLE(sample_seurat_objects,merged_seurat_object) \
		| CISTOPIC_PER_SAMPLE

}

