#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Evaluate how well an algorithm recapitulates the exons sensitive to changes in the activity of
# splicing factor
#
# Script outline
# --------------
# - define true positives and true negatives
# - compute TPR, FPR, Precision, Recall
#

import argparse
import pandas as pd
import gc
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
from scipy import stats

# variables
N_JOBS = 1
THRESH_DPSI = 5
THRESH_DPSI_REL = 50 # 50% possible PSI change
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_viper_splicing'
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,"results","sf_targets_inference")
inference_file = os.path.join(RESULTS_DIR,"files","target_inference","correlation.tsv.gz")
ground_truth_kd_dpsi_file = os.path.join(PREP_DIR,'ground_truth_kd','ENCORE','delta_psi-EX.tsv.gz')
ground_truth_kd_rel_file = os.path.join(PREP_DIR,'ground_truth_kd','ENCORE','delta_psi_rel-EX.tsv.gz')
n_jobs=10
"""

##### FUNCTIONS #####
def load_data(inference_file, ground_truth_kd_dpsi_file, ground_truth_kd_rel_file):
    # read
    inference = pd.read_table(inference_file, index_col=0)
    ground_truth_kd_dpsi = pd.read_table(ground_truth_kd_dpsi_file, index_col=0)
    ground_truth_kd_rel = pd.read_table(ground_truth_kd_rel_file, index_col=0)

    # subset
    common_events = set(inference.index).intersection(ground_truth_kd_dpsi.index)
    common_splicing_factors = set(inference.columns).intersection(
        ground_truth_kd_dpsi.columns
    )
    inference = inference.loc[common_events, common_splicing_factors].copy()
    ground_truth_kd_dpsi = ground_truth_kd_dpsi.loc[
        common_events, common_splicing_factors
    ].copy()
    ground_truth_kd_rel = ground_truth_kd_rel.loc[
        common_events, common_splicing_factors
    ].copy()

    gc.collect()

    return inference, ground_truth_kd_dpsi, ground_truth_kd_rel


def evaluate_prediction_kd(targets_true, targets_pred):
    """
    evaluate how well the algorithm recapitulates the interactions
    """

    threshs_rel = [THRESH_DPSI_REL] # only interested in changes higher than 50% possible
    threshs_class = np.arange(0, 1, 0.01)
    evaluation = []
    for thresh_rel in threshs_rel:
        # define ground truth
        y_true = -np.sign(targets_true)
        y_true[np.abs(targets_true) < thresh_rel] = 0

        for thresh_class in threshs_class:
            # define classes
            y_pred = np.sign(targets_pred)
            y_pred[np.abs(targets_pred) < thresh_class] = 0

            # compute frequencies of
            n_correct = (y_pred == y_true).sum()
            n_incorrect = (y_pred != y_true).sum()
            n_total = len(y_pred)
            prop_correct = n_correct / n_total
            prop_incorrect = n_incorrect / n_total

            # compute rate scores
            ## TPR: proportion of targets correctly classified
            is_true_target = y_true != 0
            tpr = (
                y_pred[is_true_target] == y_true[is_true_target]
            ).sum() / is_true_target.sum()
            ## TNR: proportion of not targets correctly identified
            is_true_nottarget = y_true == 0
            tnr = (
                y_pred[is_true_nottarget] == y_true[is_true_nottarget]
            ).sum() / is_true_nottarget.sum()
            ## FPR: proportion of targets incorrectly identified
            is_pred_target = y_pred != 0
            fpr = (
                y_pred[is_pred_target] != y_true[is_pred_target]
            ).sum() / is_pred_target.sum()
            ## Precision: from the predicted targets, how many are true targets
            is_pred_target = y_pred != 0
            precision = (
                y_pred[is_pred_target] == y_true[is_pred_target]
            ).sum() / is_pred_target.sum()
            ## F1 score
            recall = tpr
            f1 = 2 * precision * recall / (precision + recall)
            
            # compute chi squared test
            tab = pd.crosstab(y_true, y_pred)
            try:
                chisq_stat, chisq_pvalue, _, _ = stats.chi2_contingency(tab)
            except:
                chisq_stat, chisq_pvalue = np.nan, np.nan
            
            # save
            eval_i = pd.Series(
                {
                    "n_correct": n_correct,
                    "n_incorrect": n_incorrect,
                    "n_total": n_total,
                    "prop_correct": prop_correct,
                    "prop_incorrect": prop_incorrect,
                    "n_tp_pred": (
                        y_pred[is_true_target] == y_true[is_true_target]
                    ).sum(),
                    "total_tp": is_true_target.sum(),
                    "tpr": tpr,
                    "n_tn_pred": (
                        y_pred[is_true_nottarget] == y_true[is_true_nottarget]
                    ).sum(),
                    "total_tn": is_true_nottarget.sum(),
                    "tnr": tnr,
                    "n_pred_not_tp": (
                        y_pred[is_pred_target] != y_true[is_pred_target]
                    ).sum(),
                    "n_pred_tp": (
                        y_pred[is_pred_target] != y_true[is_pred_target]
                    ).sum(),
                    "total_pred": is_pred_target.sum(),
                    "fpr": fpr,
                    "precision": precision,
                    "recall": recall,
                    "f1": f1,
                    "threshold_classification": thresh_class,
                    "threshold_dpsi_rel": str(thresh_rel),
                    "chisq_stat": chisq_stat,
                    "chisq_pvalue": chisq_pvalue,
                    "chisq_log10_pvalue": -np.log10(chisq_pvalue)
                }
            )
            evaluation.append(eval_i)

    evaluation = pd.DataFrame(evaluation)

    return evaluation


def evaluate_targets_single(targets_pred, targets_true, targets_dpsi, ensembl, thresh_dpsi=THRESH_DPSI):
    # drop missing values
    to_keep = targets_pred.notna() & targets_true.notna() & targets_dpsi.notna()
    targets_pred = targets_pred[to_keep]
    targets_true = targets_true[to_keep]
    targets_dpsi = targets_dpsi[to_keep]

    # we don't trust changes smaller than |dPSI|<5, consider them 0
    targets_true[np.abs(targets_dpsi) < thresh_dpsi] = 0

    # evaluate target predictions
    evaluation = evaluate_prediction_kd(targets_true, targets_pred)
    evaluation["KD_ENSEMBL"] = ensembl
    
    return evaluation

    
def evaluate_targets(
    inference, ground_truth_kd_dpsi, ground_truth_kd_rel, n_jobs
):
    
    evaluations = Parallel(n_jobs=n_jobs)(
        delayed(evaluate_targets_single)(
            inference[sf_oi], ground_truth_kd_rel[sf_oi], ground_truth_kd_dpsi[sf_oi], sf_oi
        )
        for sf_oi in tqdm(inference.columns)
    )
    evaluations = pd.concat(evaluations)
    
    return evaluations


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inference_file", type=str)
    parser.add_argument("--ground_truth_kd_dpsi_file", type=str)
    parser.add_argument("--ground_truth_kd_rel_file", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    inference_file = args.inference_file
    ground_truth_kd_dpsi_file = args.ground_truth_kd_dpsi_file
    ground_truth_kd_rel_file = args.ground_truth_kd_rel_file
    n_jobs = args.n_jobs
    output_file = args.output_file

    print(args)

    # load
    print("Loading data...")
    inference, ground_truth_kd_dpsi, ground_truth_kd_rel = load_data(
        inference_file, ground_truth_kd_dpsi_file, ground_truth_kd_rel_file
    )
    
    print("Evaluating targets...")
    evaluations = evaluate_targets(inference, ground_truth_kd_dpsi, ground_truth_kd_rel, n_jobs)

    # save
    print("Saving data...")
    evaluations.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
