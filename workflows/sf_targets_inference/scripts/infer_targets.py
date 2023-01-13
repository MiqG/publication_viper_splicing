#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# For each drug, measure the shortest path length from each of its targets
# to each of the genes significantly associated to the drug sensitivity profile

import argparse
import pandas as pd
import gc
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

# variables
N_JOBS = 1
SAVE_PARAMS = {"sep":"\t", "compression":"gzip", "index":False}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_viper_splicing'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
SUPPORT_DIR = os.path.join(ROOT,"support")
splicing_file = os.path.join(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
genexpr_file = os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
gene_sources_file = os.path.join(SUPPORT_DIR,"GOBP_RNA_SPLICING-ensembl.txt")
n_jobs = 10
method = "correlation"
"""
##### FUNCTIONS #####
def load_data(splicing_file, genexpr_file, gene_sources_file):
    # read
    splicing = pd.read_table(splicing_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)
    gene_sources = list(pd.read_table(gene_sources_file, header=None)[0])

    # subset
    genes_oi = set(gene_sources).intersection(genexpr.index)
    genexpr = genexpr.loc[genes_oi].copy()

    # check order
    common_samples = set(splicing.columns).intersection(genexpr.columns)
    splicing = splicing[common_samples].copy()
    genexpr = genexpr[common_samples].copy()

    gc.collect()

    return splicing, genexpr


def compute_correlation_single(splicing, genexpr_single):
    correl = splicing.T.corrwith(genexpr_single.T, method="spearman")
    correl.name = genexpr_single.name

    return correl


def compute_correlations(splicing, genexpr, n_jobs):
    correls = Parallel(n_jobs=n_jobs)(
        delayed(compute_correlation_single)(splicing, genexpr.loc[gene_oi])
        for gene_oi in tqdm(genexpr.index)
    )
    result = pd.concat(correls, axis=1)

    # add n observations
    result["n_obs"] = splicing.notna().sum(axis=1)

    return result


def infer_targets(splicing, genexpr, method, n_jobs):

    if method == "correlation":
        result = compute_correlations(splicing, genexpr, n_jobs)

    result = result.reset_index()

    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--splicing_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--gene_sources_file", type=str)
    parser.add_argument("--method", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    splicing_file = args.splicing_file
    genexpr_file = args.genexpr_file
    gene_sources_file = args.gene_sources_file
    method = args.method
    n_jobs = args.n_jobs
    output_file = args.output_file

    print(args)

    # load
    print("Loading data...")
    splicing, genexpr = load_data(splicing_file, genexpr_file, gene_sources_file)

    print("Inferring targets...")
    result = infer_targets(splicing, genexpr, method, n_jobs)

    # save
    print("Saving data...")
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
