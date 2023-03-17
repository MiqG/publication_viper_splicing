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
from aracne import ARACNe, FindThreshold

# variables
N_JOBS = 1
THRESH_FDR = 0.05
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}
ARACNE_BIN = "~/repositories/ARACNe-AP/dist/aracne.jar"
N_ARACNE_BOOTSTRAPS = 5  # 100

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_viper_splicing'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
SUPPORT_DIR = os.path.join(ROOT,"support")
targets_file = os.path.join(PREP_DIR,"event_psi_imputed","CCLE-EX.tsv.gz")
upstream_regs_file = os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
upstream_regs_oi_file = os.path.join(SUPPORT_DIR,"GOBP_RNA_SPLICING-ensembl.txt")
correlation_spearman_file = None
n_jobs = 10
method = "aracne_py"
"""
##### FUNCTIONS #####
def load_data(
    targets_file, upstream_regs_file, upstream_regs_oi_file, correlation_spearman_file
):
    # read
    targets = pd.read_table(targets_file, index_col=0)
    upstream_regs = pd.read_table(upstream_regs_file, index_col=0)
    upstream_regs_oi = list(pd.read_table(upstream_regs_oi_file, header=None)[0])

    # index names
    targets.index.name = "target"
    upstream_regs.index.name = "upstream_regulator"

    # subset
    genes_oi = set(upstream_regs_oi).intersection(upstream_regs.index)
    upstream_regs = upstream_regs.loc[genes_oi].copy()

    # check order
    common_samples = set(targets.columns).intersection(upstream_regs.columns)
    targets = targets[common_samples].copy()
    upstream_regs = upstream_regs[common_samples].copy()

    # drop events and genes with no variation
    targets = targets.loc[targets.std(1) > 1]
    upstream_regs = upstream_regs.loc[upstream_regs.std(1) > 1]

    # load correlation spearman
    if correlation_spearman_file is not None:
        correlation_spearman = pd.read_table(correlation_spearman_file)
    else:
        correlation_spearman = None

    gc.collect()

    return targets, upstream_regs, correlation_spearman


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


def compute_correlation_single(targets, upstream_reg_single, method):
    correl = targets.apply(
        lambda x: corr(x, upstream_reg_single, method), axis=1
    ).dropna()
    correl["upstream_regulator"] = upstream_reg_single.name
    correl = correl.reset_index()

    # prepare regulon
    ## likelihood
    correl["likelihood"] = np.abs(correl["statistic"])
    ## tfmode
    correl["tfmode"] = correl["statistic"]
    ## keep only significant correlations
    correl["padj"] = np.nan
    _, fdr, _, _ = multitest.multipletests(
        correl.loc[np.isfinite(correl["pvalue"]), "pvalue"], method="fdr_bh"
    )
    correl.loc[np.isfinite(correl["pvalue"]), "padj"] = fdr
    ## filter out irrelevant associations
    correl = correl.loc[correl["padj"] < THRESH_FDR].copy()

    return correl


def compute_correlations(targets, upstream_regs, n_jobs, method):
    correls = Parallel(n_jobs=n_jobs)(
        delayed(compute_correlation_single)(targets, upstream_regs.loc[reg_oi], method)
        for reg_oi in tqdm(upstream_regs.index)
    )
    result = pd.concat(correls)
    result = result.reset_index()

    return result


def execute_aracne(targets, upstream_regs, aracne_bin, n_jobs, output_dir, step, random_seed):
    # write aracne inputs
    upstream_regs_file = os.path.join(output_dir, "upstream_regs.tsv")
    upstream_regulators_file = os.path.join(output_dir, "upstream_regulators.txt")
    targets_file = os.path.join(output_dir, "targets.tsv")

    if (
        (not os.path.isfile(upstream_regs_file))
        | (not os.path.isfile(upstream_regulators_file))
        | (not os.path.isfile(targets_file))
    ):
        print("Writing inputs in", output_dir, "...")
        ## SF gene expression
        upstream_regs.reset_index().to_csv(upstream_regs_file, sep="\t", index=None)
        ## list of SFs
        pd.DataFrame(upstream_regs.index).to_csv(
            upstream_regulators_file, sep="\t", index=None, header=False
        )
        ## event PSI (imputed)
        targets.dropna().reset_index().to_csv(targets_file, sep="\t", index=None)

    # run aracne
    if step=="aracne_threshold":
        ## find MI threshold
        print("Finding MI threshold...")
        cmd = """
        java -Xmx5G -jar {aracne_bin} \
                --expfile_upstream {upstream_regs} \
                --tfs {upstream_regulators} \
                --expfile_downstream {targets} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --seed {random_seed} \
                --threads {n_jobs} \
                --calculateThreshold
        """.format(
            aracne_bin=aracne_bin,
            upstream_regs=os.path.join(output_dir, "upstream_regs.tsv"),
            upstream_regulators=os.path.join(output_dir, "upstream_regulators.txt"),
            targets=os.path.join(output_dir, "targets.tsv"),
            output_dir=output_dir,
            n_jobs=n_jobs,
            random_seed=random_seed
        )
        exit = os.system(cmd)
        assert exit == 0
    
    elif step=="aracne_bootstrap":
        ## run aracne with boostraps
        print("Running ARACNE bootstraps..")
        cmd = """
        java -Xmx5G -jar {aracne_bin} \
                --expfile_upstream {upstream_regs} \
                --tfs {upstream_regulators} \
                --expfile_downstream {targets} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --seed {random_seed} \
                --threads {n_jobs}
        """.format(
            aracne_bin=aracne_bin,
            upstream_regs=os.path.join(output_dir, "upstream_regs.tsv"),
            upstream_regulators=os.path.join(output_dir, "upstream_regulators.txt"),
            targets=os.path.join(output_dir, "targets.tsv"),
            output_dir=output_dir,
            n_jobs=n_jobs,
            random_seed=random_seed,
        )
        print(cmd)
        exit = os.system(cmd)
        assert exit == 0
        
    
    return None


def execute_aracne_full(
    targets, upstream_regs, aracne_bin, n_jobs, n_bootstraps, correlation_spearman
):
    with tempfile.TemporaryDirectory() as tmpdirname:
        print("created temporary directory", tmpdirname)

        # write aracne inputs
        print("Writing inputs in", tmpdirname, "...")
        ## SF gene expression
        upstream_regs.reset_index().to_csv(
            os.path.join(tmpdirname, "upstream_regs.tsv"), sep="\t", index=None
        )
        ## list of SFs
        pd.DataFrame(upstream_regs.index).to_csv(
            os.path.join(tmpdirname, "upstream_regulators.txt"),
            sep="\t",
            index=None,
            header=False,
        )
        ## event PSI (imputed)
        targets.dropna().reset_index().to_csv(
            os.path.join(tmpdirname, "targets.tsv"), sep="\t", index=None
        )

        # run aracne
        ## find MI threshold
        print("Finding MI threshold...")
        cmd = """
        java -Xmx5G -jar {aracne_bin} \
                --expfile_upstream {upstream_regs} \
                --tfs {upstream_regulators} \
                --expfile_downstream {targets} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --seed 1 \
                --threads {n_jobs} \
                --calculateThreshold
        """.format(
            aracne_bin=aracne_bin,
            upstream_regs=os.path.join(tmpdirname, "upstream_regs.tsv"),
            upstream_regulators=os.path.join(tmpdirname, "upstream_regulators.txt"),
            targets=os.path.join(tmpdirname, "targets.tsv"),
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
                --expfile_upstream {upstream_regs} \
                --tfs {upstream_regulators} \
                --expfile_downstream {targets} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --seed $i \
                --threads {n_jobs}
        i=$(( $i + 1 ))
        done
        """.format(
            aracne_bin=aracne_bin,
            upstream_regs=os.path.join(tmpdirname, "upstream_regs.tsv"),
            upstream_regulators=os.path.join(tmpdirname, "upstream_regulators.txt"),
            targets=os.path.join(tmpdirname, "targets.tsv"),
            output_dir=tmpdirname,
            n_jobs=n_jobs,
            n_bootstraps=n_bootstraps,
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
    ## "upstream_regulator" and "target" columns
    result["upstream_regulator"] = result["Regulator"]
    result["target"] = result["Target"]
    ## likelihood
    result["likelihood"] = result["MI"]
    ## tfmode
    result = pd.merge(
        result,
        correlation_spearman[["upstream_regulator", "target", "tfmode"]],
        on=["upstream_regulator", "target"],
        how="left",
    )
    result.loc[result["tfmode"].isnull(), "tfmode"] = 0

    return result


def run_aracne_py(targets, upstream_regs, n_jobs):

    arac = ARACNe(
        n_jobs=n_jobs,
        seed=RANDOM_SEED,
        n_iterations=10,  ######
        consolidate_threshold=1.1,
        # threshold_class = FindThreshold
    )
    regulons = arac.fit_transform(upstream_regs, targets)

    # prepare regulon
    ## likelihood
    result["likelihood"] = result["mutual_information_ap"]
    ## tfmode
    result["tfmode"] = result["spearman_correlation"]

    return regulons


def infer_targets(
    targets, upstream_regs, method, n_jobs, aracne_bin, correlation_spearman, output_dir, random_seed
):

    if method == "correlation_spearman":
        result = compute_correlations(
            targets, upstream_regs, n_jobs, "correlation_spearman"
        )

    elif method == "correlation_pearson":
        result = compute_correlations(
            targets, upstream_regs, n_jobs, "correlation_pearson"
        )

    elif method == "aracne_ap":
        assert correlation_spearman is not None

        result = execute_aracne_full(
            targets,
            upstream_regs,
            aracne_bin,
            n_jobs,
            N_ARACNE_BOOTSTRAPS,
            correlation_spearman,
        )

    elif method == "aracne_py":
        result = run_aracne_py(targets, upstream_regs, n_jobs)

    elif method == "aracne_threshold":
        result = execute_aracne(targets, upstream_regs, aracne_bin, n_jobs, output_dir, method, random_seed)
        
    elif method == "aracne_bootstrap":
        result = execute_aracne(targets, upstream_regs, aracne_bin, n_jobs, output_dir, method, random_seed)

    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets_file", type=str)
    parser.add_argument("--upstream_regs_file", type=str)
    parser.add_argument("--upstream_regs_oi_file", type=str)
    parser.add_argument(
        "--correlation_spearman_file",
        type=str,
        default=None,
        help="only for method=aracne",
    )
    parser.add_argument("--method", type=str)
    parser.add_argument("--random_seed", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)
    parser.add_argument("--output_file", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default=None)
    parser.add_argument("--aracne_bin", type=str, default=ARACNE_BIN)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    targets_file = args.targets_file
    upstream_regs_file = args.upstream_regs_file
    upstream_regs_oi_file = args.upstream_regs_oi_file
    method = args.method
    n_jobs = args.n_jobs
    output_file = args.output_file
    output_dir = args.output_dir
    aracne_bin = args.aracne_bin
    random_seed = args.random_seed
    correlation_spearman_file = args.correlation_spearman_file
    
    print(args)

    # load
    print("Loading data...")
    targets, upstream_regs, correlation_spearman = load_data(
        targets_file,
        upstream_regs_file,
        upstream_regs_oi_file,
        correlation_spearman_file,
    )
    
    print("Inferring targets...")
    result = infer_targets(
        targets.iloc[:100],
        upstream_regs.iloc[:100],
        method,
        n_jobs,
        aracne_bin,
        correlation_spearman,
        output_dir = output_dir,
        random_seed = random_seed
    )

    # save
    if result is not None:
        print("Saving data...")
        result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
