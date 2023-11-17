// Declare syntax version
//module load nextflow/21.10.1 (on exacloud)
nextflow.enable.dsl=2

  //////////////////////////////
 ///	Script Parameters	///////
//////////////////////////////
params.proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
params.outdir = "${params.proj_dir}/nf_analysis"
params.sample_dir="${params.proj_dir}/cellranger_data"
params.ref = "${params.proj_dir}/ref"
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
		path(sample_dir)

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

process MERGE_SAMPLES_CALLPEAKS {
	//REQUIRES MACS3 AND SAMTOOLS INSTALL IN PATH
	//Initialize Seurat Object per sample.

	input:
		path(sample_dir)
	output:
		tuple path(sample_dir),path("merged*bed")
	script:
		"""
		#merge all ATAC bam files
		samples_arr=(${sample_dir})
		for i in "\${samples_arr[@]}"; do echo \${i}"/outs/atac_possorted_bam.bam"; done > fofn.txt
		samtools cat -@ 10 -b fofn.txt | samtools sort -@ 10 - > out.bam

		#run macs2 to call atac peaks
		/home/groups/CEDAR/mulqueen/src/miniconda3/bin/macs2 \\
		callpeak -f BAMPE \\
		-t out.bam \\
		-g hs \\
		-n merged \\
		-q 0.01
		"""

}

process DIM_REDUCTION_PER_SAMPLE {
	//Generate per sample seurat object and perform dim reduction.
	//Reanalyze ATAC data with combined peak set and perform dim reduction.

	input:
		tuple path(sample_dir), path(combined_peaks)
	output:
		path("*SeuratObject.rds")
	script:
	"""
	Rscript ${params.src_dir}/seurat_dim_reduction_per_sample.R \\
	${combined_peaks} \\
	${sample_dir} \\
	${params.outdir}/plots \\
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
	Rscript ${params.src_dir}/seurat_cistopic_per_sample.R \\
	${obj_in} \\
	${params.outdir}/plots
	"""
}

process PUBLIC_DATA_LABEL_TRANSFER_PER_SAMPLE {
	//Run single-cell label trasfer using available RNA data
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		tuple path(obj_in), path(cistopic_in)
	output:
		path("${obj_in}")

	script:
	"""
	Rscript ${params.src_dir}/seurat_public_data_label_transfer.R \\
	${obj_in} \\
	${params.outdir}/plots
	"""
}


  //////////////////////////////////////////////////////
 ///	Sample Integration and Merged Processing	///
//////////////////////////////////////////////////////

process MERGED_CLUSTER {
	//Run merge seurat objects again and run LIGER on merged seurat object.
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(seurat_objects)
	output:
		tuple path("${merged_in}"),path("*.genecounts.rds")

	script:
	"""
	Rscript ${params.src_dir}/merged_cluster_liger.R \\
	"${seurat_objects}"
	"""
}

process MERGED_CHROMVAR {	
	//Run chromVAR on merged seurat object.
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(merged_in)
	output:
		path("${merged_in}")

	script:
	"""
	Rscript ${params.src_dir}/seurat_public_data_label_transfer.R \\
	${merged_in}
	"""
}

process MERGED_GENE_ACTIVITY {
	//Run Signac Gene activity function on seurat object.
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(merged_in)
	output:
		tuple path("${merged_in}"), path("*.GeneActivity.rds")

	script:
	"""
	Rscript ${params.src_dir}/seurat_public_data_label_transfer.R \\
	${merged_in}
	"""
}

//process MERGED_CELLTYPE_BARPLOTS_AND_ALLUVIAL {}

//process MOLECULAR_CHARACTERIZATION {}

//process NMF_METAPROGRAMS {}
  //////////////////////////////////
 ///	CNV Calling Per Sample	///
//////////////////////////////////
process INFERCNV_RNA_PER_SAMPLE {
	//Run InferCNV per sample.

	input:
		path(merged_in)
	output:
		path("${merged_in}")
	script:
	"""
	mkdir -p ${params.outdir}/cnv

	Rscript ${params.src_dir}/infercnv_per_sample.R \\
	${merged_in} \\
	${params.outdir}/cnv
	"""
}

process CASPER_RNA_PER_SAMPLE {
	//Run CASPER per sample.

	input:
		path(merged_in)
	output:
		path("${merged_in}")
	script:
	"""
	Rscript ${params.src_dir}/casper_per_sample.R \\
	${merged_in} \\
	${ref} \\
	${params.sample_dir} \\
	${params.outdir}/cnv
	"""
}

process COPYKAT_RNA_PER_SAMPLE {
	//Run CopyKAT per sample.

	input:
		path(merged_in)
	output:
		path("${merged_in}")
	script:
	"""
	Rscript ${params.src_dir}/copykat_per_sample.R \\
	${merged_in} \\
	${params.outdir}/cnv
	"""
}

/*
process COPYSCAT_ATAC_PER_SAMPLE {
	//Run CopyscAT per sample.

	input:
		path(merged_in)
	output:
		path("${merged_in}")
	script:
	"""
	Rscript ${params.src_dir}/copykat_per_sample.R \\
	${merged_in} \\
	${params.outdir}/cnv \\
	${ref} \\
	${params.sample_dir}

	"""
}
*/
  //////////////////////////////
 ///	Cell Type Analysis	/////
//////////////////////////////


workflow {
	/* SETTING UP VARIABLES */
		sample_dir = Channel.fromPath("${params.sample_dir}/*/" , type: 'dir').map { [it.name, it ] }

	// DATA CORRECTION
		for_generating_seurat_objects=
		SCRUBLET_RNA(sample_dir) \
		| SOUPX_RNA \
		| collect \
		| MERGE_SAMPLES_CALLPEAKS

	// DATA PROCESSING 
		merged_seurat_object =
		for_generating_seurat_objects \
		| DIM_REDUCTION_PER_SAMPLE \
		| CISTOPIC_PER_SAMPLE \
		| PUBLIC_DATA_LABEL_TRANSFER_PER_SAMPLE \
		| collect \
		| MERGED_CLUSTER \
		| MERGED_CHROMVAR \
		| MERGED_GENE_ACTIVITY
		
	// CNV CALLING 
		merged_seurat_object =
		merged_seurat_object \
		| INFERCNV_RNA_PER_SAMPLE \
		| CASPER_RNA_PER_SAMPLE \
		| COPYKAT_RNA_PER_SAMPLE 
		// COPYSCAT_ATAC_PER_SAMPLE

}

