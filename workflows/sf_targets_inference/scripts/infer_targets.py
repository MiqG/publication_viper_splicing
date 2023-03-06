#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# For each splicing-related gene, infer the exons sensitive to its activity

import argparse
import pandas as pd
import gc
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
from scipy import stats
from statsmodels.stats import multitest
import tempfile
import os

# variables
N_JOBS = 1
THRESH_FDR = 0.05
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}
ARACNE_BIN = "~/repositories/ARACNe-AP/dist/aracne.jar"
N_ARACNE_BOOTSTRAPS = 10 ##########

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_viper_splicing'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
SUPPORT_DIR = os.path.join(ROOT,"support")
splicing_file = os.path.join(PREP_DIR,"event_psi_imputed","CCLE-EX.tsv.gz")
genexpr_file = os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
gene_sources_file = os.path.join(SUPPORT_DIR,"GOBP_RNA_SPLICING-ensembl.txt")
correlation_spearman_file = None
n_jobs = 10
method = "correlation_spearman"
"""
##### FUNCTIONS #####
def load_data(
    splicing_file, genexpr_file, gene_sources_file, correlation_spearman_file
):
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

    # drop events and genes with no variation
    splicing = splicing.loc[splicing.std(1) > 0]
    genexpr = genexpr.loc[genexpr.std(1) > 0]
    
    # log-transform gene expression
    genexpr = np.log2(genexpr + 1)
    
    # normalize
    #splicing = (splicing - splicing.mean(1).values.reshape(-1, 1)) / splicing.std(
    #    1
    #).values.reshape(-1, 1)
    #genexpr = (genexpr - genexpr.mean(1).values.reshape(-1, 1)) / genexpr.std(
    #    1
    #).values.reshape(-1, 1)

    # load correlation spearman
    if correlation_spearman_file is not None:
        correlation_spearman = pd.read_table(correlation_spearman_file)
    else:
        correlation_spearman = None

    gc.collect()

    return splicing, genexpr, correlation_spearman


def corr(a, b, method):
    idx = np.isfinite(a) & np.isfinite(b)
    a = a[idx]
    b = b[idx]

    corr_funcs = {
        "correlation_spearman": stats.spearmanr,
        "correlation_pearson": stats.pearsonr,
    }
    corr_func = corr_funcs[method]

    try:
        statistic, pvalue = corr_func(a, b)
        result = pd.Series(
            {
                "statistic": statistic,
                "pvalue": pvalue,
                "n_samples": len(a),
                "method": method,
            }
        )
    except:
        result = pd.Series(
            {
                "statistic": np.nan,
                "pvalue": np.nan,
                "n_samples": np.nan,
                "method": method,
            }
        )

    return result


def compute_correlation_single(splicing, genexpr_single, method):
    correl = splicing.apply(lambda x: corr(x, genexpr_single, method), axis=1).dropna()
    correl["GENE"] = genexpr_single.name
    correl = correl.reset_index()

    # prepare regulon
    ## likelihood
    correl["likelihood"] = np.abs(correl["statistic"])
    ## tfmode
    correl["tfmode"] = correl["statistic"]
    ## "splicing_factor" and "target" columns
    correl["splicing_factor"] = correl["GENE"]
    correl["target"] = correl["EVENT"]
    ## keep only significant correlations
    correl["padj"] = np.nan
    _, fdr, _, _ = multitest.multipletests(
        correl.loc[np.isfinite(correl["pvalue"]), "pvalue"], method="fdr_bh"
    )
    correl.loc[np.isfinite(correl["pvalue"]), "padj"] = fdr
    ## filter out irrelevant splicing changes
    correl = correl.loc[correl["padj"] < THRESH_FDR].copy()

    return correl


def compute_correlations(splicing, genexpr, n_jobs, method):
    correls = Parallel(n_jobs=n_jobs)(
        delayed(compute_correlation_single)(
            splicing, genexpr.loc[gene_oi], method
        )
        for gene_oi in tqdm(genexpr.index)
    )
    result = pd.concat(correls)

    return result


def execute_aracne(
    splicing, genexpr, aracne_bin, n_jobs, n_bootstraps, correlation_spearman
):
    with tempfile.TemporaryDirectory() as tmpdirname:
        print("created temporary directory", tmpdirname)

        # write aracne inputs
        print("Writing inputs in", tmpdirname, "...")
        ## SF gene expression
        genexpr.reset_index().to_csv(
            os.path.join(tmpdirname, "genexpr.tsv"), sep="\t", index=None
        )
        ## list of SFs
        pd.DataFrame(genexpr.index).to_csv(
            os.path.join(tmpdirname,"splicing_factors.txt"), sep="\t", index=None, header=False
        )
        ## event PSI (imputed)
        splicing.dropna().reset_index().to_csv(
            os.path.join(tmpdirname, "splicing.tsv"), sep="\t", index=None
        )

        # run aracne
        ## find MI threshold
        print("Finding MI threshold...")
        cmd = """
        java -Xmx5G -jar {aracne_bin} \
                --expfile_upstream {genexpr} \
                --tfs {splicing_factors} \
                --expfile_downstream {splicing} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --seed 1 \
                --threads {n_jobs} \
                --calculateThreshold
        """.format(
            aracne_bin=aracne_bin,
            genexpr=os.path.join(tmpdirname, "genexpr.tsv"),
            splicing_factors=os.path.join(tmpdirname, "splicing_factors.txt"),
            splicing=os.path.join(tmpdirname, "splicing.tsv"),
            output_dir=tmpdirname,
            n_jobs=n_jobs,
        )
        exit = os.system(cmd)
        assert exit == 0
        ## run aracne with boostraps
        print("Running ARACNE bootstraps..")
        cmd = """
        i=1
        while [ $i -le {n_bootstraps} ]
        do
        java -Xmx5G -jar {aracne_bin} \
                --expfile_upstream {genexpr} \
                --tfs {splicing_factors} \
                --expfile_downstream {splicing} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --seed $i \
                --threads {n_jobs}
        i=$(( $i + 1 ))
        done
        """.format(
            aracne_bin=aracne_bin,
            genexpr=os.path.join(tmpdirname, "genexpr.tsv"),
            splicing_factors=os.path.join(tmpdirname, "splicing_factors.txt"),
            splicing=os.path.join(tmpdirname, "splicing.tsv"),
            output_dir=tmpdirname,
            n_jobs=n_jobs,
            n_bootstraps=n_bootstraps
        )
        print(cmd)
        exit = os.system(cmd)
        assert exit == 0
        ## consolidate bootstraps
        print("Consolidating bootstraps...")
        cmd = """
        java -Xmx5G -jar {aracne_bin} \
                --output {output_dir} \
                --threads {n_jobs} \
                --consolidate
         """.format(
            aracne_bin=aracne_bin, output_dir=tmpdirname, n_jobs=n_jobs
        )
        exit = os.system(cmd)
        assert exit == 0

        # load result
        result = pd.read_table(os.path.join(tmpdirname, "network.txt"))
        
        print("Shape result:", result.shape)

    # prepare regulon
    ## "splicing_factor" and "target" columns
    result["splicing_factor"] = result["Regulator"]
    result["target"] = result["Target"]
    ## likelihood
    result["likelihood"] = result["MI"]
    ## tfmode
    result = pd.merge(
        result,
        correlation_spearman[["splicing_factor", "target", "tfmode"]],
        on=["splicing_factor", "target"],
        how="left",
    )
    result.loc[result["tfmode"].isnull(), "tfmode"] = 0

    return result


def infer_targets(splicing, genexpr, method, n_jobs, aracne_bin, correlation_spearman):

    if method == "correlation_spearman":
        result = compute_correlations(splicing, genexpr, n_jobs, "correlation_spearman")

    elif method == "correlation_pearson":
        result = compute_correlations(splicing, genexpr, n_jobs, "correlation_pearson")

    elif method == "aracne":
        assert correlation_spearman is not None

        result = execute_aracne(
            splicing,
            genexpr,
            aracne_bin,
            n_jobs,
            N_ARACNE_BOOTSTRAPS,
            correlation_spearman,
        )

    result = result.reset_index()

    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--splicing_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--gene_sources_file", type=str)
    parser.add_argument(
        "--correlation_spearman_file",
        type=str,
        default=None,
        help="only for method=aracne",
    )
    parser.add_argument("--method", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--aracne_bin", type=str, default=ARACNE_BIN)

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
    aracne_bin = args.aracne_bin
    correlation_spearman_file = args.correlation_spearman_file

    print(args)

    # load
    print("Loading data...")
    splicing, genexpr, correlation_spearman = load_data(
        splicing_file, genexpr_file, gene_sources_file, correlation_spearman_file
    )

    print("Inferring targets...")
    result = infer_targets(
        splicing, genexpr, method, n_jobs, aracne_bin, correlation_spearman
    )

    # save
    print("Saving data...")
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
