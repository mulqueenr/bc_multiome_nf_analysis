#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(reshape2)
library(optparse)
library(circlize)
library(RColorBrewer)

option_list = list(
  make_option(c("-s", "--object_input"), type="character", default=NULL, 
              help="Seurat Object", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#cistopic=readRDS(opt$cistopic)
#titan=readRDS(opt$titan)
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3")
opt$object_input="/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects/merged.chromvar.SeuratObject.rds"
dat=readRDS(opt$object_input)

met<-dat@meta.data
met<-met[!duplicated(met$sample_ID),]

#met$sample_weight
#met$Diagnosis
#met$Mol_Diagnosis
#met$sampled_site

#transferring these from sample metadata sheet.
mass=c("DCIS_01"="0.24", "DCIS_02"="0.45", "DCIS_03"="0.19", "IDC_01"="0.90", "IDC_02"="0.70", "IDC_03"="0.27", "IDC_04"="0.18", "IDC_14"="0.89", "IDC_16"="1.06", "ILC_04"="0.92", "IDC_05"="0.27", "IDC_06"="0.16", "IDC_07"="0.19", "IDC_08"="0.18", "IDC_09"="0.21", "IDC_10"="0.93", "IDC_11"="0.08", "IDC_12"="0.11", "IDC_13"="0.98", "IDC_15"="1.38", "ILC_01"="0.23", "ILC_02"="1.76", "ILC_03"="1.08", "ILC_05"="0.77", "NAT_04"="0.39", "NAT_11"="0.31", "NAT_14"="0.67")
age=c('DCIS_01'='31', 'DCIS_02'='49', 'DCIS_03'='61', 'IDC_01'='75', 'IDC_10'='68', 'IDC_11'='NA', 'IDC_12'='NA', 'IDC_02'='51', 'IDC_03'='74', 'IDC_04'='67', 'IDC_05'='34', 'IDC_06'='76', 'IDC_07'='44', 'IDC_08'='63', 'IDC_09'='63', 'ILC_01'='57', 'ILC_02'='64','NAT_11'='37', 'NAT_14'='50', 'NAT_04'='67', 'IDC_13'='68', 'IDC_14'='40', 'IDC_15'='43', 'IDC_16'='75', 'ILC_03'='71', 'ILC_04'='65', 'ILC_05'='34')
sampled_site=c('DCIS_01'='Primary', 'DCIS_02'='Primary', 'DCIS_03'='Primary', 'IDC_01'='Primary', 'IDC_02'='Primary', 'IDC_03'='Primary', 'IDC_04'='Primary', 'IDC_14'='Primary', 'IDC_16'='Primary', 'ILC_04'='Primary', 'IDC_05'='Primary', 'IDC_06'='Primary', 'IDC_07'='Primary', 'IDC_08'='Primary', 'IDC_09'='Primary', 'IDC_10'='Primary', 'IDC_11'='Primary', 'IDC_12'='Primary', 'IDC_13'='Primary', 'IDC_15'='Primary', 'ILC_01'='Primary', 'ILC_02'='Metastasis of Lymph Node', 'ILC_03'='Primary', 'ILC_05'='Primary', 'NAT_04'='NAT', 'NAT_11'='NAT', 'NAT_14'='NAT')
plot_order=c('DCIS_01'='1', 'DCIS_02'='2', 'DCIS_03'='3', 'IDC_01'='4', 'IDC_02'='5', 'IDC_03'='6', 'IDC_04'='7', 'IDC_14'='8', 'IDC_16'='9', 'ILC_04'='10', 'IDC_05'='11', 'IDC_06'='12', 'IDC_07'='13', 'IDC_08'='14', 'IDC_09'='15', 'IDC_10'='16', 'IDC_11'='17', 'IDC_12'='18', 'IDC_13'='19', 'IDC_15'='20', 'ILC_01'='21', 'ILC_02'='22', 'ILC_03'='23', 'ILC_05'='24', 'NAT_04'='25', 'NAT_11'='26', 'NAT_14'='27') 
met$wgs=c('DCIS_01'='0', 'DCIS_02'='1', 'DCIS_03'='1', 'IDC_01'='1', 'IDC_02'='1', 'IDC_03'='1', 'IDC_04'='1', 'IDC_14'='0', 'IDC_16'='0', 'ILC_04'='0', 'IDC_05'='0', 'IDC_06'='1', 'IDC_07'='1', 'IDC_08'='1', 'IDC_09'='1', 'IDC_10'='1', 'IDC_11'='1', 'IDC_12'='0', 'IDC_13'='0', 'IDC_15'='0', 'ILC_01'='0', 'ILC_02'='0', 'ILC_03'='0', 'ILC_05'='0', 'NAT_04'='0', 'NAT_11'='1', 'NAT_14'='1')
met$spatial_atac=c('DCIS_01'='0', 'DCIS_02'='0', 'DCIS_03'='0', 'IDC_01'='1', 'IDC_02'='1', 'IDC_03'='0', 'IDC_04'='0', 'IDC_14'='0', 'IDC_16'='0', 'ILC_04'='0', 'IDC_05'='0', 'IDC_06'='1', 'IDC_07'='0', 'IDC_08'='0', 'IDC_09'='1', 'IDC_10'='0', 'IDC_11'='0', 'IDC_12'='0', 'IDC_13'='0', 'IDC_15'='0', 'ILC_01'='0', 'ILC_02'='0', 'ILC_03'='0', 'ILC_05'='0', 'NAT_04'='0', 'NAT_11'='0', 'NAT_14'='0')
met$mass<-as.numeric(mass[met$sample])
met$age<-as.numeric(age[met$sample])
met$sampled_site<-sampled_site[met$sample]
met$plot_order<-as.numeric(plot_order[met$sample])
met$multiome<-1
met$wgs<-0
met$spatial_atac<-0

sample_heatmap<-met[c("sample","age","mass","Diagnosis","Mol_Diagnosis","sampled_site","plot_order","multiome","wgs","spatial_atac")]
row.names(sample_heatmap)<-sample_heatmap$sample
sample_heatmap<-sample_heatmap[order(sample_heatmap$plot_order),]
sample_heatmap<-sample_heatmap[c("age","mass","Diagnosis","Mol_Diagnosis","sampled_site","multiome","wgs","spatial_atac")]


hist_col=c("NAT"="#99CCFF","DCIS"="#CCCCCC","IDC"="#FF9966","ILC"="#006633")
clin_col=c("IDC ER+/PR−/HER2+"="#FFCC66", "DCIS"="#CCCCCC", "IDC ER+/PR−/HER2−"="#FFCCCC", "IDC ER+/PR+/HER2−"="#FF6699", "ILC ER+/PR−/HER2−"="#66CC99", "ILC ER+/PR+/HER2−"="#006633")
sampled_col=c("Primary"="#8A4C80","Metastasis of Lymph Node"="#4c9173","NAT"="#99CCFF")
age_col=colorRamp2(c(min(sample_heatmap$age),max(sample_heatmap$age)),c("#d789d7","#2a3d66"))
mass_col=colorRamp2(c(min(sample_heatmap$mass),max(sample_heatmap$mass)),c("#f2fc9f","#b05977"))

assay_col=c("0"="white","1"="black")

ha = rowAnnotation(age=sample_heatmap$age,
                      mass=sample_heatmap$mass,
                      histological_type=sample_heatmap$Diagnosis,
                      molecular_type=sample_heatmap$Mol_Diagnosis,
                      sampled_site=sample_heatmap$sampled_site,
                      col = list(age=age_col,
                                  mass=mass_col,
                                    histological_type =hist_col,
                                    clinical_subtype=clin_col,
                                    sampled_site=sampled_col))

pdf("sample_metadata.heatmap.pdf")
plt<-Heatmap(sample_heatmap[c("multiome","wgs","spatial_atac")],
 cluster_columns=F,cluster_rows=F,
 left_annotation=ha,
 col=assay_col)
print(plt)
dev.off()
