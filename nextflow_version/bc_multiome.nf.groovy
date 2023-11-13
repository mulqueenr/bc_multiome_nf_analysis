// Declare syntax version
//https://github.com/danrlu/Nextflow_cheatsheet/blob/main/nextflow_cheatsheet.pdf
//module load nextflow/21.10.1 (on exacloud)
nextflow.enable.dsl=2
// Script parameters
params.proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
params.outdir = "${params.proj_dir}/nf_analysis"
params.ref = "${params.proj_dir}/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
params.sample_dir="${params.proj_dir}/sequencing_data"

log.info """

		================================================
		    Breast Cancer Multiome NF PIPELINE v1.0
		================================================
		NF Working Directory : ${workflow.launchDir}
		Reference Directory: ${params.ref}
		Project Directory : ${params.proj_dir}
		Looking for samples in: ${params.sample_dir}
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE 
process CELLRANGER_ARC { 
	//Generate Undetermined Fastq Files from BCL Files.
	//Based on https://github.com/nf-core/modules/tree/master/modules/nf-core/cellranger/count
	publishDir "${params.outdir}/cellranger_out/", mode: 'copy', overwrite: true

	input:
		path sample
	output:
		tuple val("${sample.simpleName}"), path("**/outs/**")
	script:
		"""
		cellranger-arc count --id=${sample.simpleName} \\
		--reference=${params.ref} \\
		--libraries=${sample} \\
		--localcores=${task.cpus} \\
		--localmem=90
		"""
}

/* //Data correction.
process SCRUBLET_RNA {
	//Perform scrublet on raw RNA count.
	publishDir "${params.outdir}/cellranger_out/", mode: 'copy', overwrite: true


	input:
		path val(sample),path(cellranger_output)
	script:
		"""
		scrublet_per_sample.py ${cellranger_output}
		"""

}

process SOUPX_RNA {
	//Perform soupX on raw RNA counts.

	publishDir "${params.outdir}/cellranger_out/", mode: 'copy', overwrite: true


	input:
		path val(sample),path(cellranger_output)
	script:
		"""
		soupx_per_sample.py ${sample} ${cellranger_output}
		"""
}

*/
//process SEURAT_GENERATION {}

//process MERGE_SAMPLES_CALLPEAKS {}

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
		def sample = Channel.fromPath( "${params.sample_dir}/NAT_11.csv" )

	/* CELLRANGER PIPELINE */
		cellranger_output=CELLRANGER_ARC(sample) //\
		//| SCRUBLET_RNA \
		//| SOUPX_RNA



}

