// Declare syntax version
//module load nextflow/21.10.1 (on exacloud)
nextflow.enable.dsl=2

  //////////////////////////////
 ///	Script Parameters	///////
//////////////////////////////
params.proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"
params.outdir = "${params.proj_dir}/nf_analysis"
params.sample_dir="${params.proj_dir}/cellranger_data/second_round" 
params.ref = "${params.proj_dir}/ref"
params.src_dir="${params.proj_dir}/bc_multiome_nf_analysis/src"
params.force_rewrite="false"
//params.merged_bed="${params.proj_dir}/merged_500bp.bed"
params.sample_metadata="${params.proj_dir}/sample_metadata.csv"

log.info """

		================================================
		    Breast Cancer Multiome NF PIPELINE v1.0
		================================================
		NF Working Directory : ${workflow.launchDir}
		Project Directory : ${params.proj_dir}
		Script Directory: ${params.src_dir}
		Reference Directory: ${params.ref}
		Singularity Container: /home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif
		Looking for samples in: ${params.sample_dir}
		Force Rewriting Data: ${params.force_rewrite}
		Supplied Sample Metadata: ${params.sample_metadata}
		Supplied Bed Files of Peaks: ${params.merged_bed}
		If null, will generate peaks.
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
	//TODO: Update with R getopts library

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
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	cpus 20

	input:
		path(sample_dir)
	output:
		path("*.nf.bed")
	script:
		"""
		#merge all ATAC bam files
		samples_arr=(${sample_dir})
		for i in "\${samples_arr[@]}"; do echo \${i}"/outs/atac_possorted_bam.bam"; done > fofn.txt
		samtools cat -@ ${task.cpus} -b fofn.txt \\
		| samtools sort -@ ${task.cpus} - > out.bam

		#run macs2 to call atac peaks
		macs2 \\
		callpeak -f BAMPE \\
		-t out.bam \\
		-g hs \\
		-n merged \\
		-q 0.01

		#format as bam and filter chr
		awk 'OFS="\\t" {print \$1,\$2,\$3}' merged_peaks.narrowPeak | grep "chr" | grep -v "chrY" > merged.nf.bed
		
		#take summits and then expand to 250bp in either direction
		#awk 'OFS="\\t" {print \$1,int(\$2)-250,int(\$2)+250}' merged_summits.bed > peaks_500bp.bed
		#awk 'OFS="\\t" sub(/-./,\"1\")1' peaks_500bp.bed > peaks_500bp.nonneg.bed #set negative values to 1
		#bed will be filtered to only main chroms and no negative values.
		#merge peaks that overlap
		#bedtools merge -i peaks_500bp.nonneg.bed > mergedpeaks_500bp.nf.bed
		"""

}


process SUPPLIED_MERGED_PEAKS {
		//Copy supplied bed file. If one is given to the --merged_peaks argument on initialization of pipeline.
		publishDir "${params.outdir}", mode: 'copy', overwrite: true
		input:
			path(merged_bed)
		output:
			path("${merged_bed.simpleName}.nf.bed")
		script:
		"""
		cp -f ${merged_bed} ${merged_bed.simpleName}.nf.bed
		"""

}

process DIM_REDUCTION_PER_SAMPLE {
	//Generate per sample seurat object and perform dim reduction.
	//Reanalyze ATAC data with combined peak set and perform dim reduction.
	//Note at this stage the seurat object is unfiltered.

	input:
		tuple val(sample_name), path(sample_dir)
		path(combined_peaks)
	output:
		path("*SeuratObject.rds")
	script:
	"""
	Rscript ${params.src_dir}/seurat_dim_reduction_per_sample.R \\
	-p ${combined_peaks} \\
	-s ${sample_dir} \\
	-o ${params.outdir}
	"""
}

process MERGED_PUBLIC_DATA_LABEL_TRANSFER {
	//Run single-cell label trasfer using available RNA data
    publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true

	input:
		path(seurat_objects)
		path(metadata)
	output:
		path("merged.public_transfer.SeuratObject.rds")

	script:
	"""
	Rscript ${params.src_dir}/seurat_public_data_label_transfer.R \\
	-s "${seurat_objects}" \\
	-r ${params.ref} \\
	-m ${metadata} \\
	-o ${params.outdir}/plots
	"""
}


process CISTOPIC_PER_SAMPLE {
	//Run cisTopic on sample ATAC data

  cpus 3
	publishDir "${params.outdir}/seurat_objects/cistopic", mode: 'copy', overwrite: true

	input:
		val(merged_in)
		tuple val(sample_name), path(sample_dir)
	output:
		path("${sample_name}.cistopic.SeuratObject.rds")
	script:
	"""
	Rscript ${params.src_dir}/seurat_cistopic_per_sample_onmerged.R \\
	-i ${merged_in} \\
	-s ${sample_name} \\
	-o ${params.outdir}/cistopic
	"""
}


process TITAN_PER_SAMPLE {
	//Run TITAN on sample RNA data

	publishDir "${params.outdir}/seurat_objects/titan", mode: 'copy', overwrite: true

	cpus 3
	input:
		val(merged_in)
		tuple val(sample_name), path(sample_dir)
	output:
		path("${sample_name}.titan.SeuratObject.rds")
	script:
	"""
	Rscript ${params.src_dir}/seurat_titan_per_sample_onmerged.R \\
	-i ${merged_in} \\
	-s ${sample_name} \\
	-o ${params.outdir}/titan
	"""
}


//process INTEGRATE_TITAN_CISTOPIC_FACTORS {
	//Combine TITAN and cisTOPIC output factors
//}


  //////////////////////////////////////////////////////
 ///	Sample Integration and Merged Processing	///
//////////////////////////////////////////////////////

process MERGED_CLUSTER {
	//Run merge seurat objects again and run LIGER on merged seurat object.
  cpus 20

	input:
		path(merged_in)
	output:
		path("*.liger.SeuratObject.rds")

	script:
	"""
	Rscript ${params.src_dir}/merged_cluster_liger.R \\
	-i "${merged_in}" \\
	-o ${params.outdir}/plots
	"""
}

process MERGED_CHROMVAR {	
	//Run chromVAR on merged seurat object.

	input:
		path(merged_in)
	output:
		path("*.chromvar.SeuratObject.rds")

	script:
	"""
	Rscript ${params.src_dir}/chromvar_merged_samples.R \\
	-i ${merged_in}
	"""
}

process MERGED_GENE_ACTIVITY {
	//Run Signac Gene activity function on seurat object.
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*.rds"

	input:
		path(merged_in)

	output:
		path("*.geneactivity.SeuratObject.rds")

	script:
	"""
	Rscript ${params.src_dir}/geneactivity_merged_sample.R \\
	-i ${merged_in}
	"""
}

//process MERGED_CELLTYPE_BARPLOTS_AND_ALLUVIAL {}

//process MOLECULAR_CHARACTERIZATION {}

//process NMF_METAPROGRAMS {}


  //////////////////////////////
 ///	Cell Type Analysis	/////
//////////////////////////////

workflow {
	/* SETTING UP VARIABLES */
		sample_dir = Channel.fromPath("${params.sample_dir}/[I|N|D]*/" , type: 'dir').map { [it.name, it ] }
		sample_metadata = Channel.fromPath("${params.sample_metadata}")
		//Sample_dir finds all folders that start with I (IDC/ILC), N (NAT), or D (DCIS), but that regex filter can be removed//

	// DATA CORRECTION
		merged_peaks_input=
		SCRUBLET_RNA(sample_dir) \
		| SOUPX_RNA 
		
		if ( params.merged_bed ) {
			//Merged bed file supplied
			merged_peaks = 
			Channel.fromPath("${params.merged_bed}") \
			| SUPPLIED_MERGED_PEAKS \
			| toList
  	} else {
  		//Make merged bed file of peaks
	  	merged_peaks = 
	  	merged_peaks_input \
	  	| collect \
	    | MERGE_SAMPLES_CALLPEAKS
  	}


	// DATA PROCESSING 
		seurat_object_list =
		DIM_REDUCTION_PER_SAMPLE(sample_dir,merged_peaks) \
		| collect

		merged_seurat_object=
		MERGED_PUBLIC_DATA_LABEL_TRANSFER(seurat_object_list,sample_metadata)

		cistopic_object_list=CISTOPIC_PER_SAMPLE(merged_seurat_object,sample_dir) 
		titan_object_list=TITAN_PER_SAMPLE(merged_seurat_object,sample_dir)

		merged_seurat_object \
		| MERGED_CLUSTER \
		| MERGED_CHROMVAR \
		| MERGED_GENE_ACTIVITY
		
}
/*
#Example running
cd /home/groups/CEDAR/mulqueen/bc_multiome #move to project directory
git clone https://github.com/mulqueenr/bc_multiome_nf_analysis.git #pull github repo

module load singularity #load singularity
module load nextflow #load nextflow
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"

#run nextflow with defaults
nextflow run bc_multiome_nf_analysis/nextflow_version/bc_multiome.nf.groovy \
-with-dag bc_multiome.flowchart.png \
-with-report bc_multiome.report.html \
-with-singularity $sif \
-resume


#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

#Error in GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) : 
  Please install biovizBase
https://www.bioconductor.org/packages/biovizBase/

*/


