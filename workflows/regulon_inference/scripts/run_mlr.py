import os
import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import gc
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats

# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["OMP_NUM_THREADS"] = "1"

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

##### FUNCTIONS #####
def fit_single(x_genexpr, y_splicing, thresh_pvalue, thresh_fdr):
    X = x_genexpr
    y = y_splicing

    # drop observations with missing values
    is_missing = X.isnull().any(axis=1) | y.isnull()
    X = X.loc[~is_missing]
    y = y[~is_missing]

    X = X - X.median(axis=0).values.reshape(1,-1)
    X["intercept"] = 1

    try:
        # Fit the linear regression model
        model_full = sm.OLS(y, X).fit()

        # save
        result = pd.DataFrame({
            "EVENT": y_splicing.name,
            "splicing_factor": model_full.params.index,
            "assoc_coef": model_full.params.values,
        })
        result = result.loc[~(result["splicing_factor"]=="intercept")].copy()
        result["zscore_coef"] = (result["assoc_coef"] - result["assoc_coef"].mean()) / result["assoc_coef"].std()
        result["zscore_pvalue"] = stats.norm.sf(abs(result["zscore_coef"]))*2 
        result["zscore_fdr"] = np.nan
        _, fdr, _, _ = multipletests(result.loc[~result["zscore_pvalue"].isnull(), "zscore_pvalue"], method="fdr_bh")
        result.loc[~result["zscore_pvalue"].isnull(),"zscore_fdr"] = fdr

    except:
        result = pd.Series({
            "EVENT": y_splicing.name,
            "splicing_factor": x_genexpr.columns,
            "assoc_coef": np.nan,
            "zscore_coef": np.nan,
            "zscore_pvalue": np.nan,
            "zscore_fdr": np.nan
        })
        result = pd.DataFrame(result).T
    
    # keep only significant interactions
    if thresh_fdr is not None:
        result = result.loc[result["zscore_fdr"] <= thresh_fdr].copy() 
    elif thresh_pvalue is not None:
        result = result.loc[result["zscore_pvalue"] <= thresh_pvalue].copy() 
    
    gc.collect()

    return result


def fit_all(regulators, targets, thresh_pvalue, thresh_fdr, n_jobs, batch_size=10000):
    
    n_targets = targets.shape[0]
    n_batches = int(np.ceil(n_targets / batch_size))
    
    results = []
    for i in range(n_batches):
        # subset batch
        start, end = i*batch_size, i*batch_size+batch_size
        
        # fit
        results_batch = Parallel(n_jobs=n_jobs)(
            delayed(fit_single)(
                x_genexpr = regulators.T, # regulators as features
                y_splicing = targets.loc[event_oi],
                thresh_pvalue = thresh_pvalue,
                thresh_fdr = thresh_fdr
            )
            for event_oi in tqdm(targets.index[start:end])
        )
        results_batch = pd.concat(results_batch)
        results.append(results_batch)
        
    results = pd.concat(results)
    
    return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regulators_file', type=str)
    parser.add_argument('--regulators_oi_file', type=str)
    parser.add_argument('--targets_file', type=str, default=None)
    parser.add_argument('--n_jobs', type=int)
    parser.add_argument('--thresh_pvalue', type=float, default=None)
    parser.add_argument('--thresh_fdr', type=float, default=None)
    parser.add_argument('--output_file', type=str)
    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    regulators_file = args.regulators_file
    regulators_oi_file = args.regulators_oi_file
    targets_file = args.targets_file
    n_jobs = args.n_jobs
    thresh_pvalue = args.thresh_pvalue
    thresh_fdr = args.thresh_fdr
    output_file = args.output_file
    
    # load
    regulators = pd.read_table(regulators_file, index_col=0)
    regulators_oi = pd.read_table(regulators_oi_file, header=None)[0].tolist()
    targets = pd.read_table(targets_file, index_col=0)
    gc.collect()

    # subset
    ## regulators of interest
    regulators = regulators.loc[regulators.index.isin(regulators_oi)].copy()
    ## samples
    common_samples = set(regulators.columns).intersection(targets.columns)
    regulators = regulators[common_samples].dropna().copy()
    targets = targets[common_samples].copy()
    
    # drop events and genes with too many missing values
    targets = targets.loc[~targets.isnull().all(axis=1)].copy()
    regulators = regulators.loc[~regulators.isnull().all(axis=1)].copy()
    
    # drop events and genes with no variation
    targets = targets.loc[targets.std(axis=1) > 0].copy()
    regulators = regulators.loc[regulators.std(axis=1) > 0].copy()
    
    # fit models
    results = fit_all(regulators, targets, thresh_pvalue, thresh_fdr, n_jobs)
    
    # prep outputs
    results["regulator"] = results["splicing_factor"]
    results["target"] = results["EVENT"]
    results["likelihood"] = np.abs(results["assoc_coef"])
    results["tfmode"] = np.sign(results["assoc_coef"])

    # save
    results.to_csv(output_file, **SAVE_PARAMS)    
    
    
##### SCRIPT #####
if __name__ == '__main__':
    main()
    print('Done!')    
