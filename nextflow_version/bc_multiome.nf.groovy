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
params.merged_bed=null
params.sample_metadata="${params.proj_dir}/sample_metadata.csv"
log.info """

		================================================
		    Breast Cancer Multiome NF PIPELINE v1.2
		================================================
		NF Working Directory : ${workflow.launchDir}
		Project Directory : ${params.proj_dir}
		Script Directory: ${params.src_dir}
		Reference Directory: ${params.ref}
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
	//Output back into cellranger out directory
	 label 'scrub'
	 containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		tuple val(sample_name), path(sample_dir)

	script:
		"""
		python /src/1_preprocessing_scrublet_per_sample.py \\
		-m ${sample_dir}/outs/filtered_feature_bc_matrix.h5 \\
		-o ${sample_dir}/outs
		"""

}

process SOUPX_RNA {
	//Perform soupX on raw RNA counts.
	//Output back into cellranger out directory
  label 'inhouse'
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"

	input:
		tuple val(sample_name), path(sample_dir)
	output:
		path(sample_dir)

	script:
		"""
		Rscript /src/2_preprocessing_soupx_per_sample.R \\
		-o ${sample_name} \\
		-d ${sample_dir}
		"""
}

  //////////////////////////////////////
 ///		Peak Calling			///
//////////////////////////////////////

process MERGE_SAMPLES_CALLPEAKS {
	//TODO update SIF to include macs3 and samtools 
	//REQUIRES MACS3 AND SAMTOOLS INSTALL IN PATH 
	//Initialize Seurat Object per sample.
	//NOTE merged_fragments are not properly cellIDs for merging across all samples. 
	//Use passqc.atac_fragments.tsv.gz for each cell ranger object for actual data upload
	publishDir "${params.outdir}/peaks", mode: 'copy', overwrite: true, pattern: "merged.nf.bed"
	publishDir "${params.outdir}/peaks", mode: 'copy', overwrite: true, pattern: "merged_fragments.sorted.tsv.gz"

	cpus 30
	//label 'inhouse'
	input:
		path(sample_dir)
	output:	
		path("merged.nf.bed"), emit: merged_bed
		path("merged_fragments.sorted.tsv.gz"), emit: merged_fragments
	script:
		"""
		#merge all ATAC fragment files
		samples_arr=(${sample_dir})

		for i in "\${samples_arr[@]}"; 
		do zcat \${i}"/outs/atac_fragments.tsv.gz" | \\
		awk -v sample=\${i} 'OFS="\\t" {print \$1,\$2,\$3,\$4,\$5,sample}'; done >> merged_fragments.tsv 

		pigz -p ${task.cpus} merged_fragments.tsv

		zcat merged_fragments.tsv.gz | \\
		grep -v "^#" | \\
		sort --parallel=${task.cpus} -T . -k1,1 -k2,2n -k3,3n - > merged_fragments.sorted.tsv

		pigz -p ${task.cpus} merged_fragments.sorted.tsv

		#run macs3 to call atac peaks
		#note because this is using bed, it does not account for cell id specifics, or counts in extra columns,
		#so basically calling on deduplicated fragments
				
		macs3 callpeak -f BEDPE \\
		--tempdir . \\
		-t merged_fragments.sorted.tsv.gz \\
		-g hs -n merged \\
		-B -q 0.01	

		#format as bam and filter chr
		awk 'OFS="\\t" {print \$1,\$2,\$3}' merged_peaks.narrowPeak | grep "chr" | grep -v "chrY" > merged.nf.bed
		"""

}

  //////////////////////////////////////
 ///	Seurat Sample Processing	///
//////////////////////////////////////

process MERGE_SAMPLES_AND_FILTER {
	//Run single-cell label trasfer using available RNA data
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*rds"
    publishDir "${params.outdir}/plots/sample_qc", mode: 'copy', overwrite: true, pattern: "*pdf"

	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir},${params.sample_dir}"
	label 'inhouse'
	input:
		path(sample_dir)
		path(merged_peaks)
		path(metadata)
	output:
		path("1_merged.unfiltered.SeuratObject.rds"), emit: unfiltered_obj
		path("2_merged.scrublet_count_filtered.SeuratObject.rds"), emit: obj
		path("*pdf"), emit: plots

	script:
	"""
	Rscript /src/3_preprocessing_seurat_sample_filtering.R \\
	-s . \\
	-p ${merged_peaks} \\
	-m ${metadata}
	"""
	
}

process MERGED_PUBLIC_DATA_LABEL_TRANSFER {
	//Run single-cell label trasfer using available RNA data
	//All reference data must have a metadata column of "celltype" to label transfer
	//Data processing conducted in the SIF is included in /public_data/ folder
  containerOptions "--bind ${params.src_dir}:/src/,${params.outdir},${params.ref}:/ref/"
  label 'inhouse'
  publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*rds"
  publishDir "${params.outdir}/plots/label_transfer", mode: 'copy', overwrite: true, pattern: "*pdf"

	
	input:
		path(obj)
	output:
		path("3_merged.public_transfer.SeuratObject.rds"), emit: seurat_obj
		path("*pdf"), emit: public_transfer_plots

	script:
	"""
	Rscript /src/4_preprocessing_seurat_public_data_label_transfer.R \\
	-i ${obj} \\
	-r /ref/
	"""
}


process MERGED_CHROMVAR {	
	//Run chromVAR on merged seurat object.
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'inhouse'
	publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*rds"

	input:
		path(merged_in)
	output:
		path("4_merged.chromvar.SeuratObject.rds")

	script:
	"""
	Rscript /src/5_preprocessing_chromvar_merged_samples.R \\
	-i ${merged_in}
	"""
}

process MERGED_GENE_ACTIVITY {
	//Run Signac Gene activity function on seurat object.
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.outdir}"
	publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*.rds"
	label 'inhouse'

	input:
		path(merged_in)

	output:
		path("5_merged.geneactivity.SeuratObject.rds")

	script:
	"""
	Rscript /src/6_preprocessing_geneactivity_merged_sample.R  \\
	-i ${merged_in}
	"""
}


  //////////////////////////////////////////////////////
 ///	Sample Clustering and Merged Processing		///
//////////////////////////////////////////////////////

process FIG1_MERGED_CLUSTER {
	//Run merge seurat objects again and run LIGER on merged seurat object.
	//This script is not generalizable. It includes hard coded metadata and a set seet for assigning clusters.
	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*.rds"
	publishDir "${params.outdir}/plots/fig1", mode: 'copy', overwrite: true, pattern: "*.pdf"

	label 'inhouse'
	input:
		path(merged_in)
	output:
		path("6_merged.celltyping.SeuratObject.rds"), emit: cluster_obj
		path("*pdf"), emit: fig_pdf

	script:
	"""
	Rscript /src/7_cluster_data.R \\
	-i ${merged_in}
	"""
}

process FIG2_PSEUDOBULK_ANALYSIS {
	//currently uses some local libraries, and so I bound that directory. to be fixed in the SIF
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.ref}:/ref/,${params.outdir},/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3"
	publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*.rds"
	publishDir "${params.outdir}/plots/fig2", mode: 'copy', overwrite: true, pattern: "*.pdf"
	publishDir "${params.outdir}/plots/fig2", mode: 'copy', overwrite: true, pattern: "*.tsv", optional: true

	label 'inhouse'
	input:
		path(merged_in)
	output:
		path("6_merged.celltyping.SeuratObject.rds"), emit: tf_obj, includeInputs: true
		path("*pdf"), emit: fig_pdf
		path("*tsv"), emit: fig_tsv, optional: true

	script:
	"""
	Rscript /src/8_tf_markers_coverage_plots.R \\
	-i ${merged_in} \\
	-r /ref/nakshatri/nakshatri_multiome.rds
	"""
}

process FIG2_SCSUBTYPE {
	//currently uses some local libraries, and so I bound that directory. to be fixed in the SIF
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.outdir},/home/users/mulqueen/R/x86_64-conda-linux-gnu-library/4.3"
	publishDir "${params.outdir}/seurat_objects", mode: 'copy', overwrite: true, pattern: "*.rds"
	publishDir "${params.outdir}/plots/fig2", mode: 'copy', overwrite: true, pattern: "*.pdf"

	label 'inhouse'
	input:
		path(merged_in)
	output:
		path("7_merged.scsubtype.SeuratObject.rds"), emit: scsubtype_obj
		path("*pdf"), emit: fig_pdf

	script:
	"""
	Rscript /src/9_scsubtype_merged_sample.R \\
	-i ${merged_in}
	"""
}


process SCENICPLUS_FORMATTING_FROM_SIGNAC {
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.outdir}"
	publishDir "${params.outdir}/scenic_input", mode: 'copy', overwrite: true, pattern: "*.{csv,mtx}"

	label 'inhouse'
	input:
		path(merged_in)
	output:
		tuple path("metadata.atac.csv"),path("frag_paths.csv"),
				path("atac_counts.cells.csv"),path("atac_counts.peaks.csv"),
				path("atac_counts.mtx"), emit: scenic_all_cells_atac

		tuple path("metadata.rna.csv"),path("rna_counts.cells.csv"),
				path("rna_counts.genes.csv"),path("rna_counts.mtx"), emit: scenic_all_cells_rna

		tuple path("epi_metadata.atac.csv"),path("epi_frag_paths.csv"),
				path("epi_atac_counts.cells.csv"),path("epi_atac_counts.peaks.csv"),
				path("epi_atac_counts.mtx"), emit: scenic_epi_atac

		tuple path("epi_metadata.rna.csv"),path("epi_rna_counts.cells.csv"),
				path("epi_rna_counts.genes.csv"),path("epi_rna_counts.mtx"), emit: scenic_epi_rna
	script:
	"""
	Rscript /src/10_1_scenic_input_from_signac.R \\
	-i ${merged_in} 
	"""
}


process SCENICPLUS_RNA_PREPROCESSING {
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.outdir}"
	publishDir "${params.outdir}/scenic_input/rna", mode: 'copy', overwrite: true, pattern: "*.{h5ad}"

	label 'scenic'
	input:
		tuple path(meta),path(cells),path(genes),path(counts)

		
	output:
		path("*.h5ad")
	script:
	"""
	#make tmp directories for numba processing
	export NUMBA_CACHE_DIR=\${PWD}/.tmp
	export MPLCONFIGDIR=\${PWD}/.tmp
	mkdir -p \$NUMBA_CACHE_DIR

	python /src/10_2_scenic_rna_processing.py \\
	--rna_genes $genes \\
	--rna_cells $cells \\
	--rna_matrix $counts \\
	--meta $meta
	"""
}

process SCENICPLUS_ATAC_PREPROCESSING {
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.outdir}"
	publishDir "${params.outdir}/scenic_input/atac", mode: 'copy', overwrite: true, pattern: "*.{pkl}"
	label 'scenic'
	input:
		tuple path(meta),path(frag),path(cells),path(peaks),path(counts)
	output:
		path("*cistopic_obj.pkl"), emit: cistopic_obj
	script:
	"""
	export NUMBA_CACHE_DIR=\${PWD}/.tmp
	export MPLCONFIGDIR=\${PWD}/.tmp
	mkdir -p \$NUMBA_CACHE_DIR

	python /src/10_3_scenic_pycistopic_cistopicObj_generation.py \\
	--frag_path $frag  \\
	--atac_cells $cells \\
	--atac_peaks $peaks \\
	--atac_counts $counts \\
	--meta $meta
	"""
}


process SCENICPLUS_CISTARGET_ON_PEAKS {
	//based on https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html
	containerOptions "--bind ${params.proj_dir},${params.src_dir}:/src/,${params.outdir},${params.ref}:/ref/"
	publishDir "${params.outdir}/scenic_input/atac", mode: 'copy', overwrite: true, pattern: "*.{pkl}"

	label 'scenic'
	input:
		path(peaks)
		path(scriptDir), stageAs: "/create_cisTarget_databases"
	output:
		path("bc_multiome/"), emit: cistarget_db
	script:
	"""
	GENOME_FASTA="/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
	CHROMSIZES="/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/chrNameLength.txt"
	CBDIR="/ref/aertslab_motif_collection/v10nr_clust_public/singletons"

	#create peaks fasta with background
	create_cisTarget_databases/create_fasta_with_padded_bg_from_bed.sh \
			\${GENOME_FASTA} \
			\${CHROMSIZES} \
			${peaks} \
			bc_multiome.fa \
			1000 \
			yes
	
	ls \${CBDIR}/jaspar*cb > motifs.txt #list jaspar motifs

	create_cisTarget_databases/create_cistarget_motif_databases.py \
		-f bc_multiome.fa \
		-M \${CBDIR} \
		-m motifs.txt \
		-o bc_multiome \
		--bgpadding 1000 \
		-t ${task.cpus}
	"""
}

process SCENICPLUS_RUN_CISTOPIC {
	//run pycistopic per topic 
	publishDir "${params.outdir}/scenic_output/cistopic/topics", mode: 'copy', overwrite: true
	executor 'slurm'
	cpus 8
	time 16.hour
	memory 200.GB
	errorStrategy 'retry'
	maxRetries 3
	queue 'batch'

	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'scenic'

	input:
		tuple val(topic_count), path(cistopic)
	output:
		path("Topic*")
	script:
	"""
	export NUMBA_CACHE_DIR=\${PWD}/tmp
	export MPLCONFIGDIR=\${PWD}/tmp
	mkdir -p \$NUMBA_CACHE_DIR

	python /src/pycistopic_run_topic.py \\
	--topicCount ${topic_count} \\
	--cistopicObj ${cistopic} \\
	--outDir "." \\
	--taskCpus 5 \\
	--memory 200G
	"""
}


process SCENICPLUS_MERGE_CISTOPICS {
	//merge all topic models and select best. output merged cistopic object
	publishDir "${params.outdir}/scenic_output/cistopic/plots", mode: 'copy', overwrite: true, pattern="*pdf"
	publishDir "${params.outdir}/scenic_output/cistopic/topic_specificity", mode: 'copy', overwrite: true, pattern="*csv"
	publishDir "${params.outdir}/scenic_output/cistopic/region_sets/topic_bed", mode: 'copy', overwrite: true, pattern="*bed"
	publishDir "${params.outdir}/scenic_output/cistopic/region_sets/celltype_dar_bed", mode: 'copy', overwrite: true, pattern="*bed"
	publishDir "${params.outdir}/scenic_output/cistopic", mode: 'copy', overwrite: true, pattern="cistopic_obj.pkl"

	containerOptions "--bind ${params.src_dir}:/src/,${params.outdir}"
	label 'scenic'

	input:
		path(topics)
		path(cistopic_obj)
	output:
		path("*pdf"), emit: cistopic_plots
		path("*csv"), emit: cistopic_specificity
		path("cistopic_obj.pkl"), emit: cistopic_obj, includeInputs: true
		path("./region_sets/Topics_otsu/*bed"), emit: topic_bed
		path("./region_sets/DARs_cell_type/*bed"), emit: celltype_bed


	script:
	"""
	export NUMBA_CACHE_DIR=\${PWD}/tmp
	export MPLCONFIGDIR=\${PWD}/tmp
	mkdir -p \$NUMBA_CACHE_DIR

	python /src/10_4_scenic_cistopic_combine_topics.py \\
	--cistopicObj ${cistopic_obj} \\
	--outDir "./" \\
	--tmpDir "./"
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

//process  {}

//process MOLECULAR_CHARACTERIZATION {}

//process NMF_METAPROGRAMS {}


  //////////////////////////////
 ///	Cell Type Analysis	/////
//////////////////////////////

workflow {
	// SETTING UP VARIABLES
		sample_dir = \
		Channel.fromPath("${params.sample_dir}/[I|N|D]*/" , type: 'dir').map { [it.name, it ] }
		//Sample_dir finds all folders that start with I (IDC/ILC), N (NAT), or D (DCIS), but that regex filter can be removed//

		sample_metadata = \
		Channel.fromPath("${params.sample_metadata}")

	// DATA CORRECTION
		SCRUBLET_RNA(sample_dir) \
		| SOUPX_RNA \
		| set { merged_peaks_input }

	//ATAC PEAK GENERATION
	// supplied bed file, or call them from merged fragments file using macs3
		if ( params.merged_bed ) {
			merged_peaks = Channel.fromPath("${params.merged_bed}")
			}
		else {
			merged_peaks_input \
			| collect \
			| MERGE_SAMPLES_CALLPEAKS

			merged_peaks = \
			MERGE_SAMPLES_CALLPEAKS.out.merged_bed
		}

	// DATA PREPROCESSING 
		//Merge filtered seurat objects, add sample metadata
		cellranger_out = \
		merged_peaks_input \
		| collect

		merged_obj = \
		MERGE_SAMPLES_AND_FILTER(cellranger_out, merged_peaks, sample_metadata)
		
		//Merge sample, public data label transfers
		merged_seurat_object = \
		merged_obj.obj \
		| MERGED_PUBLIC_DATA_LABEL_TRANSFER

		merged_seurat_object.seurat_obj \
		| MERGED_CHROMVAR \
		| MERGED_GENE_ACTIVITY \
		| FIG1_MERGED_CLUSTER \
		
		FIG2_PSEUDOBULK_ANALYSIS(FIG1_MERGED_CLUSTER.out.cluster_obj)
		FIG2_SCSUBTYPE(FIG2_PSEUDOBULK_ANALYSIS.out.tf_obj)

	//SCENIC BLOCK
		scriptDir=Channel.fromPath("${params.proj_dir}/src/create_cisTarget_databases/")
		topicList=Channel.of(15,20,25,30,35,40,45,50,55,60)

		//just run epi only for now
		//Prepare for scanpy/cistopic
		SCENICPLUS_FORMATTING_FROM_SIGNAC(FIG2_SCSUBTYPE.out.scsubtype_obj)

		//Run scanpy
		SCENICPLUS_RNA_PREPROCESSING(SCENICPLUS_FORMATTING_FROM_SIGNAC.out.scenic_epi_rna)
		
		//Run cistopic
		SCENICPLUS_ATAC_PREPROCESSING(SCENICPLUS_FORMATTING_FROM_SIGNAC.out.scenic_epi_atac)
		
		//topics =
		topics = 
		topicList.combine(SCENICPLUS_ATAC_PREPROCESSING.out.cistopic_obj) \
		| map { topic, obj -> tuple(topic, obj)}
		
		topics = SCENICPLUS_RUN_CISTOPIC(topics) \
		| collect

		SCENICPLUS_MERGE_CISTOPICS(topics,
			SCENICPLUS_ATAC_PREPROCESSING.out.cistopic_obj)

		//Make cistarget db
		//SCENICPLUS_CISTARGET_ON_PEAKS(merged_peaks,scriptDir)





}

/*
		

		//generate tuple of sample names with merged object for splitting merged seurat object in parallel
		seurat_object_list
		| flatten \
		| map{[it.name.toString().split("[.]")[0]]} \
		| combine(merged_seurat_object) \
		| set { merged_object_sample_split }

		//cistopic per sample for GRN via chromatin		
		CISTOPIC_PER_SAMPLE(merged_object_sample_split,merged_out) \
		| collect \
		| set { cistopic_object_list }
		
		//titan per sample for GRN via rna
		TITAN_PER_SAMPLE(merged_object_sample_split,merged_out) \
		| collect \
		| set { titan_object_list }
*/

/*
Example run is in README.md file.
*/



