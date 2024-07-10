#!/usr/bin/env python3

import argparse
import math
import time
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.stats as st
from scipy.stats import norm

parser = argparse.ArgumentParser(description="Input file list and variants pattern.")

parser.add_argument('-r', '--regenie_results', required=True, type=str,
                    help="String of list of regenie files.")
parser.add_argument('-f', '--phenotype_id', required=True, type=str,
                    help="Phenotype id.")
parser.add_argument('-n', '--meta_study_name', required=True, type=str,
                    help="Meta study name.")


args = parser.parse_args()

def format_table(df):
    df = df.reset_index()
    df['chromosome'] = df.ID.str.extract(r'chr(\d+|X)_', expand=True).astype(str)
    df['chromosome'] = df['chromosome'].replace(['X'], '23')
    df['chromosome'] = df['chromosome'].astype(int)
    df['position'] = df.ID.str.extract(r'_(\d+)_', expand=True).astype(int)
    df['ref'] = df.ID.str.extract(r'_([ACGTN]+)_', expand=True).astype(str)
    df['alt'] = df.ID.str.extract(r'_([ACGTN]+)$', expand=True).astype(str)
    df = df.rename(columns={'phenotype_id': 'molecular_trait_id', 'pooled_effects': 'beta',"pooled_se":"se"})
    df = df [["chromosome", "position","ref","alt","molecular_trait_id","beta","se","nlog10p","n_datasets","sample_size","ac"]]  # TODO:change!
    return df

def write_pq_file(df, file_name):
    fields = [pa.field('chromosome', pa.int32()),
    pa.field('position', pa.int32()),
    pa.field('ref', pa.string()),
    pa.field('alt', pa.string()),
    pa.field('molecular_trait_id', pa.string()),
    pa.field('beta', pa.float64()),
    pa.field('se', pa.float64()),
    pa.field('nlog10p', pa.float64()),
    pa.field('n_datasets', pa.int32()),
    pa.field('sample_size', pa.int32()),
    pa.field('ac', pa.float64())]
    my_schema = pa.schema(fields)
    table = pa.Table.from_pandas(df, schema=my_schema)
    pq.write_table(table, f'{file_name}_{phenotype_id}.parquet',compression='snappy')
    

def get_weight_array(b_se_array):
    get_w = lambda x: 1 / (x ** 2)
    return get_w(b_se_array)


def replace_zeros_w_nan(array):
    array[array == 0] = np.nan


def replace_nan_with_zero(array):
    array[np.isnan(array)] = 0


def check_if_contains_zero(array):
    if (array == 0).any():
        return True
    return False


def prepare_studies(paths):
    datasets = []
    for index, path in enumerate(paths):
        dataset = MetaStudy(str(index), path.strip())
        datasets.append(dataset)
    return datasets


class MetaStudy:
    def __init__(self, dataset_name: str, regenie_file_path: str):
        self._dataset_name = dataset_name
        self._regenie_file_path = regenie_file_path
        self._df = None
        self._b_array = None
        self._b_se_array = None
        self._w_array = None
        self._study_w_es_array = None
        self._set_dataset_df_format()

    def _get_df_from_regenie_file(self):
        self._df = pd.read_csv(self._regenie_file_path, sep=' ', compression="gzip",
                            index_col='ID')

    def _rename_columns(self):
        self._df.rename(columns={'BETA': f'{self._dataset_name}_b', 'SE': f'{self._dataset_name}_b_se',
                                'N': f'{self._dataset_name}_N', 'AC': f'{self._dataset_name}_ac'}, inplace=True)

    def _select_needed_columns(self):
        if "CoLaus_tensorQTL" not in (self._regenie_file_path):
            columns = ['BETA', "SE", "N", "A1FREQ"]
            self._df = self._df[columns]
            self._df["AC"] = self._df["A1FREQ"] * self._df["N"] * 2
            self._df.drop(["A1FREQ"], axis=1, inplace=True)

    def _set_dataset_df_format(self):
        self._get_df_from_regenie_file()
        self._select_needed_columns()
        self._rename_columns()

    def get_study_df(self):
        return self._df

    def get_study_name(self):
        return self._dataset_name

    def get_weight_array(self):
        return self._w_array

    def get_w_es_array(self):
        return self._study_w_es_array

    def set_b_array(self, b_array):
        self._b_array = b_array

    def set_b_se_array(self, b_se_array):
        self._b_se_array = b_se_array

    def _calculate_w_array(self):
        study_w_array = get_weight_array(self._b_se_array)
        replace_nan_with_zero(study_w_array)
        self._w_array = study_w_array

    def _calculate_w_es_array(self):
        study_w_es_array = np.multiply(self._w_array, self._b_array)
        replace_nan_with_zero(study_w_es_array)
        self._study_w_es_array = study_w_es_array

    def calculate(self):
        self._calculate_w_array()
        self._calculate_w_es_array()


class MetaAnalysis:
    def __init__(self, study_name: str):
        self.studies = []
        self._merged_dataset_df = None
        self._meta_study_name = study_name

    def add_studies_to_analyze(self, *meta_studies: MetaStudy):
        for study in meta_studies:
            self.studies.append(study)

    def _merge_studies(self):
        self._merged_dataset_df = self.studies[0].get_study_df()
        for study in self.studies[1:]:
            self._merged_dataset_df = self._merged_dataset_df.join(study.get_study_df(),
                                                                   how='outer')


    def _select_needed_columns(self):
        selected_columns = ["id"]
        for study in self.studies:
            selected_columns.append(f"{study.get_study_name()}_b")
            selected_columns.append(f"{study.get_study_name()}_b_se")
        self._merged_dataset_df = self._merged_dataset_df[selected_columns]

    def _calculate_studies_summed_wES(self):
        studies_w_es_arrays = []
        for study in self.studies:
            studies_w_es_arrays.append(study.get_w_es_array())
        return (np.array(studies_w_es_arrays)).sum(axis=0)

    def _set_studies_working_arrays(self):
        for study in self.studies:
            study_b_array = self._merged_dataset_df[f"{study.get_study_name()}_b"].to_numpy()
            study.set_b_array(study_b_array)
            study_b_se_array = self._merged_dataset_df[f"{study.get_study_name()}_b_se"].to_numpy()
            study.set_b_se_array(study_b_se_array)
            study.calculate()

    def _calculate_studies_summed_weight(self):
        studies_weight_arrays = []
        for study in self.studies:
            studies_weight_arrays.append(study.get_weight_array())
        studies_weight_arrays_sum = (np.array(studies_weight_arrays)).sum(axis=0)
        if not check_if_contains_zero(studies_weight_arrays_sum):
            return studies_weight_arrays_sum
        else:
            replace_zeros_w_nan(studies_weight_arrays_sum)  # to avoid zerodivision
            return studies_weight_arrays_sum

    def _calculate_pooled_effects(self):
        return np.divide(self._calculate_studies_summed_wES(), self._calculate_studies_summed_weight())

    def _calculate_pooled_se(self):
        find_es_se = np.vectorize(lambda x: math.sqrt((1 / x)))
        return find_es_se(self._calculate_studies_summed_weight())

    def _calculate_p_values(self):
        ci_low, ci_upper = st.norm.interval(alpha=0.95, loc=self._merged_dataset_df['pooled_effects'].to_numpy(),
                                            scale=self._merged_dataset_df['pooled_se'].to_numpy())
        SE = (ci_upper - ci_low) / (2 * 1.96)
        Z = np.divide(np.absolute((self._merged_dataset_df['pooled_effects']).to_numpy()), SE)
        pval = 2 * (norm.sf(abs(Z)))
        self._merged_dataset_df["p_value"] = pval
        self._merged_dataset_df['nlog10p'] = - np.log10(self._merged_dataset_df['p_value'])


    def _calculate_nlog10p(self):
        effects = self._merged_dataset_df['pooled_effects'].to_numpy()
        se = self._merged_dataset_df['pooled_se'].to_numpy()
        Z = effects / se
        log_p = norm.logsf(np.abs(Z))
        log10p = -log_p / np.log(10)
        self._merged_dataset_df["nlog10p"] = log10p


    def _get_statistic_df(self):
        self._merge_studies()
        self._set_studies_working_arrays()
        self._merged_dataset_df["pooled_effects"] = self._calculate_pooled_effects()
        self._merged_dataset_df["pooled_se"] = self._calculate_pooled_se()
        self._merged_dataset_df["phenotype_id"] = phenotype_id
        #self._calculate_p_values()
        self._calculate_nlog10p()  
        self._find_missing_variants_count()
        self._calculate_sample_size()
        self._calculate_sample_AC()
        return self._merged_dataset_df[
            ["phenotype_id", 'pooled_effects', 'pooled_se', 'nlog10p', 'n_datasets', "sample_size", "ac"]]  # TODO: change!

    def _find_missing_variants_count(self):
        self._merged_dataset_df["n_datasets"] = self._merged_dataset_df.loc[:,
                                                           self._merged_dataset_df.columns.str.endswith(
                                                               '_b_se')].notnull().sum(axis=1)

    def analyze(self):
        calculated_df = self._get_statistic_df()
        calculated_table_df = format_table(calculated_df)
        write_pq_file(calculated_table_df,self._meta_study_name)
        #table = pa.Table.from_pandas(calculated_df)
        #pq.write_table(table, f'{self._meta_study_name}_{phenotype_id}.parquet')

    def _calculate_sample_size(self):
        self._merged_dataset_df["sample_size"] = self._merged_dataset_df.loc[:,
                                       self._merged_dataset_df.columns.str.endswith(
                                           '_N')].sum(axis=1)

    def _calculate_sample_AC(self):
        self._merged_dataset_df["ac"] = self._merged_dataset_df.loc[:,
                                        self._merged_dataset_df.columns.str.endswith(
                                            '_ac')].sum(axis=1)

    def _generate_variant_info(self):
        self._merged_dataset_df = self._merged_dataset_df.reset_index()
        self._merged_dataset_df['variant_position'] = self._merged_dataset_df.ID.str.extract(r'_(\d+)_').astype(int)
        self._merged_dataset_df['variant_chromosome'] = self._merged_dataset_df.ID.str.extract(r'chr(\d+|X)_').astype(str)


start = time.time()
phenotype_id = args.phenotype_id

regenie_results = args.regenie_results
paths_list = list(map(str, regenie_results.strip('[]').split(',')))
studies = prepare_studies(paths_list)
meta_study_name = args.meta_study_name
metaanalysis = MetaAnalysis(meta_study_name)
metaanalysis.add_studies_to_analyze(*studies)
metaanalysis.analyze()
end = time.time()
print("DONE!")
print(f"Time to finish workflow: {end - start}")