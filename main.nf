#!/bin/bash nextflow
nextflow.enable.dsl=2

include { METAANALYSIS } from './modules/metaanalysis'
include { FIND_SIG_RESULTS } from './modules/find_sig_results'




workflow { 
    regenie_directories = params.regenie_result_directories?.split(',') as List
    regenie_directories_regexes = regenie_directories.collect { it + "*.regenie.gz" }
    all_datasets_regenie_files = Channel.fromPath(regenie_directories_regexes )
    all_datasets_regenie_files_map = all_datasets_regenie_files.map{regenie_file -> tuple(regenie_file.name.split('\\+++_')[1].split('.regenie.gz')[0], regenie_file)}
    all_regenie_files_phenotype_grouped = all_datasets_regenie_files_map.groupTuple()
    METAANALYSIS(all_regenie_files_phenotype_grouped)
    FIND_SIG_RESULTS(METAANALYSIS.out.meta_parquet_ch)
    FIND_SIG_RESULTS
        .out
        .collectFile(storeDir: "${params.outdir}/${params.meta_study_name}/significant_meta_results", name: "${params.meta_study_name}_significant_meta.tsv", keepHeader: true) 
} 
