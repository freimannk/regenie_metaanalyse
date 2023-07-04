#!/bin/bash nextflow
//nextflow.enable.dsl=2

process METAANALYSIS {
   
    publishDir "${params.outdir}/${params.meta_study_name}/results", mode:'copy'
    container = 'quay.io/fhcrc-microbiome/python-pandas:4a6179f'


    input:
    val phenotype_id_reagenie_files
   
    output:
    path '*.parquet', emit: meta_parquet_ch

    script:
    phenotype_id = phenotype_id_reagenie_files[0]
    reagenie_files = phenotype_id_reagenie_files[1]
    
    """
        metaanalysis.py -r '${reagenie_files}' -f ${phenotype_id} -n ${params.meta_study_name}

    """
}