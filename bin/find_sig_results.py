#!/usr/bin/env python3

import argparse
import csv
import pyarrow.parquet as pq



parser = argparse.ArgumentParser(description="")

parser.add_argument('-p', '--parquet', required=True, type=str,
                    help="Parquet file path.")
parser.add_argument('-s', '--logp_threshold', required=True, type=float,
                    help="p value threshold.")
parser.add_argument('-n', '--file_name', required=True, type=str,
                    help="File name.")


args = parser.parse_args()
pq_table = pq.read_table(args.parquet , filters=[("nlog10p", ">=", args.logp_threshold)])
df = pq_table.to_pandas()

if not df.empty:
    df.to_csv(f"{args.file_name}_sig_results.txt", sep="\t", quoting=csv.QUOTE_NONE,
                index=True)