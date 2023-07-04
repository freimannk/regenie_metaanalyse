#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process FIND_SIG_RESULTS{  
    container = 'quay.io/fhcrc-microbiome/python-pandas:4a6179f'

    input:
    file parquet_file

    output:
    file '*_sig_results.txt' optional true

    script:
        """
            find_sig_results.py -p ${parquet_file} -s ${params.filter_logp_threshold} -n ${parquet_file.simpleName}


        """
}