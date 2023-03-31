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
import statsmodels.api as sm
from statsmodels.stats import multitest
from sklearn.model_selection import train_test_split
import tempfile
import os

# variables
N_JOBS = 1
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}
ARACNE_BIN = "~/repositories/ARACNe-AP/dist/aracne.jar"
TEST_SIZE = 0.15
N_ITER = 1

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_viper_splicing'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
SUPPORT_DIR = os.path.join(ROOT,"support")
targets_file = os.path.join(PREP_DIR,"event_psi_imputed","CCLE-EX.tsv.gz")
regulators_file = os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
regulators_oi_file = os.path.join(SUPPORT_DIR,"GOBP_RNA_SPLICING-ensembl.txt")
correlation_spearman_file = None
n_jobs = 10
method = "aracne_py"
"""
##### FUNCTIONS #####
def pairwise(lst1, lst2, lim=np.inf):
    i = 1
    for x in lst1:
        for y in lst2:
            if x != y:
                if i <= lim:
                    i = i + 1
                    yield (x, y)
                    
                    
def load_data(regulators_file, targets_file):
    print("Loading data...")

    # read
    regulators = pd.read_table(regulators_file, index_col=0)
    targets = pd.read_table(targets_file, index_col=0)

    # index names
    regulators.index.name = "regulator"
    targets.index.name = "target"

    # check order
    common_samples = set(targets.columns).intersection(regulators.columns)
    regulators = regulators[common_samples].copy()
    targets = targets[common_samples].copy()

    gc.collect()

    return regulators, targets


def compute_correlation_single(regulator_single, target_single, method):
    
    regulator, target = regulator_single.name, target_single.name
    
    idx = np.isfinite(regulator_single) & np.isfinite(target_single)
    a = regulator_single[idx]
    b = target_single[idx]

    corr_funcs = {
        "correlation_spearman": stats.spearmanr,
        "correlation_pearson": stats.pearsonr,
    }
    corr_func = corr_funcs[method]

    try:
        statistic, pvalue = corr_func(a, b)
        correl = pd.Series(
            {
                "regulator": regulator,
                "target": target,
                "statistic": statistic,
                "pvalue": pvalue,
                "n_samples": len(a),
            }
        )
    except:
        correl = pd.Series(
            {
                "regulator": regulator,
                "target": target,
                "statistic": np.nan,
                "pvalue": np.nan,
                "n_samples": np.nan,
            }
        )

    return correl


def compute_correlations(regulators, targets, n_jobs, method):
    pairs = pairwise(regulators.index, targets.index)
    
    correls = Parallel(n_jobs=n_jobs)(
        delayed(compute_correlation_single)(regulators.loc[reg], targets.loc[tar], method)
        for reg, tar in tqdm(pairs)
    )
    result = pd.DataFrame(correls)
    
    # adjust pvalues
    result["padj"] = np.nan
    _, fdr, _, _ = multitest.multipletests(
        result.loc[np.isfinite(result["pvalue"]), "pvalue"], method="fdr_bh"
    )
    result.loc[np.isfinite(result["pvalue"]), "padj"] = fdr
    
    # add association
    result["association"] = result["statistic"]
    
    return result


def execute_aracne(regulators_file, targets_file, aracne_bin, n_jobs):
    
    regulators = pd.read_table(regulators_file, index_col=0)
    regulators_oi = "\n".join(list(regulators.index)) 
    n_samples = regulators.shape[1]

    with tempfile.TemporaryDirectory() as tmpdirname:

        ## run aracne with boostraps
        print("Running ARACNE...")
        cmd = """
        # hack p-value threshold
        echo 0.0 > {output_dir}/miThreshold_p1E-8_samples{n_samples}.txt

        # create TFs file
        echo '{regulators_oi}' > {output_dir}/tfs.txt

        # compute mutual informations
        java -Xmx5G -jar {aracne_bin} \
                --expfile_upstream {regulators} \
                --tfs {output_dir}/tfs.txt \
                --expfile_downstream {targets} \
                --output {output_dir} \
                --pvalue 1E-8 \
                --threads {n_jobs} \
                --nobootstrap \
                --nodpi
        """.format(
            aracne_bin=aracne_bin,
            regulators=regulators_file,
            regulators_oi=regulators_oi,
            targets=targets_file,
            output_dir=tmpdirname,
            n_jobs=n_jobs,
            n_samples=n_samples
        )
        print(cmd)
        exit = os.system(cmd)
        assert exit == 0

        # load
        result = pd.read_table(os.path.join(tmpdirname, "nobootstrap_network.txt"))

    # prepare output
    result = result.rename(columns={"Regulator": "regulator", "Target": "target", "MI":"association"})

    return result


def get_summary_stats(df, col_oi):
    summary_stats = {
        col_oi + "_mean": np.mean(df[col_oi]),
        col_oi + "_median": np.median(df[col_oi]),
        col_oi + "_std": np.std(df[col_oi]),
        col_oi + "_q25": np.quantile(df[col_oi], 0.25),
        col_oi + "_q75": np.quantile(df[col_oi], 0.75),
    }
    return summary_stats


def fit_olsmodel(y, X, n_iterations):
    regulator = y.name
    target = X.columns[0]
    
    summaries = []
    for i in range(n_iterations):
        # split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE, random_state=i)

        # fit linear model to training data
        model = sm.OLS(y_train, X_train).fit()

        # log-likelihood test
        model_null = sm.OLS(y_train, X_train[["intercept"]]).fit()
        lr_stat, lr_pvalue, lr_df = model.compare_lr_test(model_null)

        # score using test data
        prediction = model.predict(X_test)
        pearson_coef, pearson_pvalue = stats.pearsonr(prediction, y_test)
        spearman_coef, spearman_pvalue = stats.spearmanr(prediction, y_test)

        # prepare output
        summary_it = {
            "iteration": i,
            "target_coefficient": model.params[target],
            "target_stderr": model.bse[target],
            "target_zscore": model.params[target] / model.bse[target],
            "target_pvalue": model.pvalues[target],
            "intercept_coefficient": model.params["intercept"],
            "intercept_stderr": model.bse["intercept"],
            "intercept_zscore": model.params["intercept"] / model.bse["intercept"],
            "intercept_pvalue": model.pvalues["intercept"],
            "n_obs": model.nobs,
            "rsquared": model.rsquared,
            "pearson_correlation": pearson_coef,
            "pearson_pvalue": pearson_pvalue,
            "spearman_correlation": spearman_coef,
            "spearman_pvalue": spearman_pvalue,
            "lr_stat": lr_stat,
            "lr_pvalue": lr_pvalue,
            "lr_df": lr_df,
        }
        summaries.append(summary_it)
        
    summaries = pd.DataFrame(summaries)

    # compute average likelihood-ratio test
    avg_lr_stat = np.mean(summaries["lr_stat"])
    avg_lr_df = np.round(summaries["lr_df"].mean())
    lr_pvalue = stats.chi2.sf(avg_lr_stat, avg_lr_df)

    # prepare output
    ## summary
    summary = {"regulator": regulator, "target": target, "n_obs": model.nobs}
    summary.update(get_summary_stats(summaries, "target_coefficient"))
    summary.update(get_summary_stats(summaries, "intercept_coefficient"))
    summary.update(get_summary_stats(summaries, "rsquared"))
    summary.update(get_summary_stats(summaries, "pearson_correlation"))
    summary.update(get_summary_stats(summaries, "spearman_correlation"))
    summary.update(get_summary_stats(summaries, "lr_stat"))
    summary.update(
        {"lr_df": lr_df, "lr_pvalue": lr_pvalue,}
    )
    summary = pd.Series(summary)

    return summary


def compute_lm_single(regulator_single, target_single, n_iterations):

    # model variables
    X = pd.DataFrame([target_single]).T
    y = regulator_single
    
    # drop missing values
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan].copy()
    y = y[~is_nan].copy()
    
    # fit
    try:
        # add intercept
        X["intercept"] = 1.0

        # fit full model
        summary = fit_olsmodel(y, X, n_iterations)

    except:
        X["intercept"] = np.nan

        # create empy summary
        summary = pd.Series(
            np.nan,
            index=[
                "regulator",
                "target",
                "n_obs",
                "target_coefficient_mean",
                "target_coefficient_median",
                "target_coefficient_std",
                "target_coefficient_q25",
                "target_coefficient_q75",
                "intercept_coefficient_mean",
                "intercept_coefficient_median",
                "intercept_coefficient_std",
                "intercept_coefficient_q25",
                "intercept_coefficient_q75",
                "rsquared_mean",
                "rsquared_median",
                "rsquared_std",
                "rsquared_q25",
                "rsquared_q75",
                "pearson_correlation_mean",
                "pearson_correlation_median",
                "pearson_correlation_std",
                "pearson_correlation_q25",
                "pearson_correlation_q75",
                "spearman_correlation_mean",
                "spearman_correlation_median",
                "spearman_correlation_std",
                "spearman_correlation_q25",
                "spearman_correlation_q75",
                "lr_stat_mean",
                "lr_stat_median",
                "lr_stat_std",
                "lr_stat_q25",
                "lr_stat_q75",
                "lr_df",
                "lr_pvalue",
            ],
        )
        summary["regulator"] = regulator_single.name
        summary["target"] = target_single.name

    # add some more info to the summary
    summary["target_mean"] = target_single.mean()
    summary["target_std"] = target_single.std()

    return summary

    
def compute_lm(regulators, targets, n_jobs, n_iterations):
    # standardize regulators
    regulators = regulators - regulators.mean(axis=1).values.reshape(-1,1)
    regulators = regulators / regulators.std(axis=1).values.reshape(-1,1)
    
    # fit linear models
    pairs = pairwise(regulators.index, targets.index)
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_lm_single)(regulators.loc[reg], targets.loc[tar], n_iterations)
        for reg, tar in tqdm(pairs)
    )
    result = pd.DataFrame(results)
    
    # adjust pvalues
    result["lr_padj"] = np.nan
    _, fdr, _, _ = multitest.multipletests(
        result.loc[np.isfinite(result["lr_pvalue"]), "lr_pvalue"], method="fdr_bh"
    )
    result.loc[np.isfinite(result["lr_pvalue"]), "lr_padj"] = fdr
    
    # add association
    result["association"] = result["target_coefficient_mean"]
    
    return result
    
    
def compute_associations(regulators, targets, method, n_jobs, aracne_bin, n_iterations):
    print("Computing associations...")

    if method == "correlation_spearman":
        result = compute_correlations(
            regulators, targets, n_jobs, "correlation_spearman"
        )

    elif method == "correlation_pearson":
        result = compute_correlations(
            regulators, targets, n_jobs, "correlation_pearson"
        )

    elif method == "aracne":
        result = execute_aracne(regulators, targets, aracne_bin, n_jobs)
        
    elif method == "linear_model":
        result = compute_lm(regulators, targets, n_jobs, n_iterations)
        
    result["method"] = method

    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets_file", type=str)
    parser.add_argument("--regulators_file", type=str)
    parser.add_argument("--method", type=str)
    parser.add_argument("--n_jobs", type=int)
    parser.add_argument("--n_iterations", type=int, default=N_ITER)
    parser.add_argument("--aracne_bin", type=str, default=ARACNE_BIN)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    targets_file = args.targets_file
    regulators_file = args.regulators_file
    method = args.method
    n_jobs = args.n_jobs
    output_file = args.output_file
    aracne_bin = args.aracne_bin
    n_iterations = args.n_iterations

    print(args)

    if method == "aracne":
        result = compute_associations(
            regulators_file,
            targets_file,
            method,
            n_jobs,
            aracne_bin,
            n_iterations
        )
    else:
        regulators, targets = load_data(regulators_file, targets_file)

        result = compute_associations(
            regulators,
            targets,
            method,
            n_jobs,
            aracne_bin,
            n_iterations
        )
    
    # save
    print("Saving data...")
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
