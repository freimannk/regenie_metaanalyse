# Meta-analysis of datasets analysed by Regenie-Nextflow-pipeline 


**Mandatory input parameters in nextflow.config file**
- params.regenie_result_directories : QTL-mapping result directories run by Regenie QTL-mapping workflow (e.g '~path1/MRCA/eQTL_analysis/regenie/,~path2/MRCE/eQTL_analysis/regenie/')
- params.meta_study_name : study name
- params.outdir : directory of the output
- params.filter_logp_threshold : -logp10 threshold
- params.phenotype_Ensembl_metadata : Ensemble metadata file in /data direcory