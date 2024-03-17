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
	 cpus 5
	 label 'scrub'
	 containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		tuple val(sample_name), path(sample_dir)

	script:
		"""
		python /src/scrublet_per_sample.py \\
		-m ${sample_dir}/outs/filtered_feature_bc_matrix.h5 \\
		-o ${sample_dir}/outs
		"""

}

process SOUPX_RNA {
	//Perform soupX on raw RNA counts.
	//TODO: Update with R getopts library
  cpus 5
  label 'inhouse'
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		path(sample_dir)

	script:
		"""
		Rscript /src/soupx_per_sample.R \\
		${sample_name} \\
		${sample_dir}
		"""
}

  //////////////////////////////////////
 ///	Seurat Sample Processing	///
//////////////////////////////////////

process MERGE_SAMPLES_CALLPEAKS {
	//REQUIRES MACS3 AND SAMTOOLS INSTALL IN PATH
	//Initialize Seurat Object per sample.
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	cpus 30
	label 'inhouse'
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

		#run macs3 to call atac peaks
		macs3 \\
		callpeak -f BAMPE \\
		-t out.bam \\
		-g hs \\
		-n merged \\
		-B \\
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
		label 'inhouse'
		input:
			path(merged_bed)
			path(sample_dir)
		output:
			path("${merged_bed}")
		script:
		"""
		touch ${merged_bed}
		"""

}

process DIM_REDUCTION_PER_SAMPLE {
	//Generate per sample seurat object and perform dim reduction.
	//Reanalyze ATAC data with combined peak set and perform dim reduction.
	//Note at this stage the seurat object is unfiltered.
  cpus 5
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		tuple val(sample_name), path(sample_dir)
		path(combined_peaks)
	output:
		path("${sample_dir.simpleName}.SeuratObject.rds")
	script:
	"""
	Rscript /src/seurat_dim_reduction_per_sample.R \\
	-p ${combined_peaks} \\
	-s ${sample_dir} \\
	-o ${params.outdir}
	"""
}

process MERGED_PUBLIC_DATA_LABEL_TRANSFER {
	//Run single-cell label trasfer using available RNA data
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true
  cpus 5
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		path(seurat_objects)
		path(metadata)
	output:
		path("merged.public_transfer.SeuratObject.rds")

	script:
	"""
	Rscript /src/seurat_public_data_label_transfer.R \\
	-s "${seurat_objects}" \\
	-r ${params.ref} \\
	-m ${metadata} \\
	-o ${params.outdir}/plots
	"""
}


process CISTOPIC_PER_SAMPLE {
	//Run cisTopic on sample ATAC data
	publishDir "${params.outdir}/seurat_objects/cistopic", mode: 'copy', overwrite: true
  maxForks 1
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		tuple val(sample_name), path(merged_in)
		path(final_object)

	output:
		path("${sample_name}.cistopic.SeuratObject.rds")
	script:
	"""
	Rscript /src/seurat_cistopic_per_sample_onmerged.R \\
	-i ${merged_in} \\
	-s ${sample_name} \\
	-o ${params.outdir}/cistopic
	"""
}


process TITAN_PER_SAMPLE {
	//Run TITAN on sample RNA data
	publishDir "${params.outdir}/seurat_objects/titan", mode: 'copy', overwrite: true
  maxForks 1
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		tuple val(sample_name), path(merged_in)
		path(final_object)

	output:
		path("${sample_name}.titan.SeuratObject.rds")
	script:
	"""
	Rscript /src/seurat_titan_per_sample_onmerged.R \\
	-i ${merged_in} \\
	-s ${sample_name} \\
	-o ${params.outdir}/titan
	"""
}


  //////////////////////////////////////////////////////
 ///	Sample Integration and Merged Processing	///
//////////////////////////////////////////////////////

process MERGED_CLUSTER {
	//Run merge seurat objects again and run LIGER on merged seurat object.
  cpus 5
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		path(merged_in)
	output:
		path("*.liger.SeuratObject.rds")

	script:
	"""
	Rscript /src/merged_cluster_liger.R \\
	-i ${merged_in} \\
	-o ${params.outdir}/plots
	"""
}

process MERGED_CHROMVAR {	
	//Run chromVAR on merged seurat object.
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		path(merged_in)
	output:
		path("*.chromvar.SeuratObject.rds")

	script:
	"""
	Rscript /src/chromvar_merged_samples.R \\
	-i ${merged_in}
	"""
}

process MERGED_GENE_ACTIVITY {
	//Run Signac Gene activity function on seurat object.
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*.rds"
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	input:
		path(merged_in)

	output:
		path("*.geneactivity.SeuratObject.rds")

	script:
	"""
	Rscript /src/geneactivity_merged_sample.R \\
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
			//Merged bed file supplied, still take merged_peaks_input to ensure QC is done on all samples before proceeding
	  	merged_peaks_input | collect
			merged_peaks = Channel.fromPath("${params.merged_bed}") 
			SUPPLIED_MERGED_PEAKS(merged_peaks,merged_peaks_input) \
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

		merged_seurat_object =
		MERGED_PUBLIC_DATA_LABEL_TRANSFER(seurat_object_list,sample_metadata)

		merged_out =
		MERGED_CLUSTER(merged_seurat_object)
		| MERGED_CHROMVAR \
		| MERGED_GENE_ACTIVITY \
		| collect

		//generate tuple of sample names with merged object 
		//for splitting merged seurat object in parallel
		merged_object_sample_split=
		seurat_object_list
		.flatten()
		.map{[it.name.toString().split("[.]")[0]]}
		.combine(merged_seurat_object)

		//merged out is to serialize these		
		cistopic_object_list=
		CISTOPIC_PER_SAMPLE(merged_object_sample_split,merged_out) \
		| collect
		
		titan_object_list =
		TITAN_PER_SAMPLE(merged_object_sample_split,merged_out) \
		| collect

}
/*
#Example running
cd /home/groups/CEDAR/mulqueen/bc_multiome #move to project directory
git clone https://github.com/mulqueenr/bc_multiome_nf_analysis.git #pull github repo

module load singularity/3.8.0 #load singularity
module load nextflow/21.10.1 #load nextflow
proj_dir="/home/groups/CEDAR/mulqueen/bc_multiome"

cd /home/groups/CEDAR/mulqueen/bc_multiome
nextflow run bc_multiome_nf_analysis/nextflow_version/bc_multiome.nf.groovy \
--force_rewrite true \
--outdir ${proj_dir}/nf_analysis_round3 \
--sample_dir ${proj_dir}/cellranger_data/third_round \
-resume


module load singularity/3.8.0 #load singularity
module load nextflow/21.10.1 #load nextflow
sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
bed="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis/merged.nf.bed" #using established bed file

*/


