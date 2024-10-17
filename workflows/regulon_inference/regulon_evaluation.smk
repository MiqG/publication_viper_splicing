import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SRC_DIR = os.path.join(ROOT,"src")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = EVENT_TYPES

PERT_SPLICING_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{omic_type}.tsv.gz')
}

EVAL_DATASETS = list(PERT_SPLICING_FILES.keys())

PERT_FILES = {
    "EX": PERT_SPLICING_FILES,
}

METADATA_FILES = [
    os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz")
]

REGULON_SETS = [
    # general comparison
    "experimentally_derived_regulons_pruned",
    "aracne_regulons_CardosoMoreira2020",
    "aracne_regulons_PANCAN_STN",
    "aracne_regulons_PANCAN_PT",
    "aracne_regulons_combined",
    "mlr_regulons_CardosoMoreira2020",
    "mlr_regulons_PANCAN_STN",
    "mlr_regulons_PANCAN_PT",
    "mlr_regulons_combined",
    "postar3_clip_regulons",
    # clip vs empirical networks
    "postar3_and_experimental_regulons",
    "experimental_without_postar3_regulons",
    "splicinglore_regulons",
    # empirical vs computational networks
    "aracne_and_experimental_regulons",
    "mlr_and_experimental_regulons",
    "aracne_and_mlr_regulons"
]

TOP_N = [100, 90, 80, 70, 60, 50, 40]
ROBUSTNESS_EVAL_SETS = ["top{N}_experimentally_derived_regulons_pruned".format(N=n) for n in TOP_N]
THRESH_DPSI = [5,10,15,20,25,30,35,40,45]
THRESHOLDS_EVAL_SETS = ["dPSIthresh{thresh}_experimentally_derived_regulons_pruned".format(thresh=t) for t in THRESH_DPSI]
REGULON_SETS = REGULON_SETS + ROBUSTNESS_EVAL_SETS + THRESHOLDS_EVAL_SETS

METHODS_ACTIVITY = ["viper","correlation_pearson","correlation_spearman","gsea"]
METHODS_ACTIVITY = {r: METHODS_ACTIVITY for r in REGULON_SETS}
METHODS_ACTIVITY["postar3_clip_regulons"] = ["viper","gsea"]

# prepare lists for zip expand
OMIC_TYPE_ZIPLIST = ["EX" for r in REGULON_SETS for m in METHODS_ACTIVITY[r] for e in EVAL_DATASETS]
REGULON_SETS_ZIPLIST = [r for r in REGULON_SETS for m in METHODS_ACTIVITY[r] for e in EVAL_DATASETS]
EVAL_DATASETS_ZIPLIST = [e for r in REGULON_SETS for m in METHODS_ACTIVITY[r] for e in EVAL_DATASETS]
METHODS_ACTIVITY_ZIPLIST = [m for r in REGULON_SETS for m in METHODS_ACTIVITY[r] for e in EVAL_DATASETS]

##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        ## run
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{method_activity}","{regulon_set}-{dataset}-{omic_type}-shadow_no-two_tailed.tsv.gz"), zip, regulon_set=REGULON_SETS_ZIPLIST, dataset=EVAL_DATASETS_ZIPLIST, method_activity=METHODS_ACTIVITY_ZIPLIST, omic_type=OMIC_TYPE_ZIPLIST),
        ## merge
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # make figures
        #os.path.join(RESULTS_DIR,"figures","regulon_evaluation"),
        #os.path.join(RESULTS_DIR,"figures","regulon_inference")
        
        
rule make_evaluation_labels:
    input:
        metadatas = METADATA_FILES
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"))
    run:
        import pandas as pd
        
        perts_oi = ["KNOCKDOWN","KNOCKOUT","OVEREXPRESSION"]
        
        for f in input.metadatas:
            # load
            metadata = pd.read_table(f)
            
            os.makedirs(output.output_dir, exist_ok=True)
            if "ENCORE" in f:
                for cell_line in metadata["cell_line"].unique():
                    # make labels
                    labels = metadata.loc[
                        metadata["cell_line"]==cell_line, 
                        ["cell_line","PERT_ENSEMBL","PERT_TYPE"]
                    ].drop_duplicates().copy()
                    
                    
                    dataset = "ENCOREKD" if "KD" in os.path.basename(f) else "ENCOREKO"
                    labels["study_accession"] = dataset
                    
                    labels["PERT_ID"] = labels[
                        ["study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"]
                    ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                    
                    # only simple perturbations
                    labels = labels.loc[labels["PERT_TYPE"].isin(perts_oi)]
                    
                    # save
                    labels.dropna().to_csv(os.path.join(output.output_dir,"%s_%s.tsv.gz") % (dataset, cell_line), **SAVE_PARAMS)
                
            elif "ENASFS" in f:
                # prepare labels
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL","PERT_TYPE"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                
                labels = metadata[["PERT_ID","PERT_GENE","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()

                # only simple perturbations
                labels = labels.loc[labels["PERT_TYPE"].isin(perts_oi)]

                # save
                labels.dropna().to_csv(os.path.join(output.output_dir,"ENASFS.tsv.gz"), **SAVE_PARAMS)
                
                
        print("Done!")
        

rule evaluate_regulons:
    input:
        signature = lambda wildcards: PERT_FILES[wildcards.omic_type][wildcards.dataset],
        regulons = os.path.join(RESULTS_DIR,"files","{regulon_set}-{omic_type}"),
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{method_activity}","{regulon_set}-{dataset}-{omic_type}-shadow_no-two_tailed.tsv.gz")
    params:
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz"),
        script_dir = SRC_DIR,
        shadow = "no",
        n_tails = "two",
        method_activity = "{method_activity}"
    conda: "publication_viper_splicing"
    threads: 1
    resources:
        # runtime = 3600*6, # h in seconds
        runtime = 60*24, # h in minutes 
        memory = 300, # GB
    shell:
        """
        nice Rscript {params.script_dir}/evaluate_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons} \
                    --eval_labels_file={params.eval_labels} \
                    --output_file={output} \
                    --method_activity={params.method_activity} \
                    --shadow_correction={params.shadow} \
                    --n_tails={params.n_tails}
        """
        
        
rule combine_evaluations:
    input:
        evaluations = [os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{method_activity}","{regulon_set}-{dataset}-{omic_type}-shadow_no-two_tailed.tsv.gz").format(regulon_set=r, dataset=d, omic_type="{omic_type}", method_activity=m) for r in REGULON_SETS for d in EVAL_DATASETS for m in METHODS_ACTIVITY[r]]
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz")
    params:
        omic_type = "{omic_type}"
    run:
        import pandas as pd
        
        print("Reading...")
        cols_oi = [
            "PERT_ID","auc_roc","auc_pr","avg_pr","mean_rank_percentile","median_rank_percentile",
            "grouping_var","n_pos_class","n_total","curves_type","eval_direction","eval_type","regulator",
            "n_networks_per_regulator","regulon_set_id","signature_id","shadow_correction","n_tails",
            "method_activity"
        ]
        evaluation = pd.concat([
            pd.read_table(f, usecols=cols_oi, low_memory=False).drop_duplicates() 
            for f in input.evaluations
        ])
        evaluation["omic_type"] = params.omic_type
        
        print("Saving...")
        evaluation.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        
    
rule figures_regulon_evaluation:
    input:
        evaluation_ex = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz"),
        regulators_per_target_robustness = os.path.join(RESULTS_DIR,"files","regulon_properties","regulators_per_target-EX.tsv.gz"),
        targets_per_regulator_robustness = os.path.join(RESULTS_DIR,"files","regulon_properties","targets_per_regulator-EX.tsv.gz"),
        regulators_per_target_thresholds = os.path.join(RESULTS_DIR,"files","regulon_properties","dPSIthresh-regulators_per_target-EX.tsv.gz"),
        targets_per_regulator_thresholds = os.path.join(RESULTS_DIR,"files","regulon_properties","dPSIthresh-targets_per_regulator-EX.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","regulon_evaluation"))
    conda: "publication_viper_splicing"
    threads: 1
    resources:
        # runtime = 3600*6, # h in seconds
        runtime = 60*6, # h in minutes 
        memory = 50, # GB
    shell:
        """
        Rscript scripts/figures_regulon_evaluation.R \
                    --evaluation_ex_file={input.evaluation_ex} \
                    --regulators_per_target_robustness_file={input.regulators_per_target_robustness} \
                    --targets_per_regulator_robustness_file={input.targets_per_regulator_robustness} \
                    --regulators_per_target_thresholds_file={input.regulators_per_target_thresholds} \
                    --targets_per_regulator_thresholds_file={input.targets_per_regulator_thresholds} \
                    --figs_dir={output}
        """
        
        
rule figures_inference_troubleshooting:
    input:
        experimental_pruned_path = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-EX"),
        aracne_and_experimental_path = os.path.join(RESULTS_DIR,"files","aracne_and_experimental_regulons-EX"),
        mlr_and_experimental_path = os.path.join(RESULTS_DIR,"files","mlr_and_experimental_regulons-EX")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","regulon_inference"))
    conda: "publication_viper_splicing"
    threads: 1
    resources:
        # runtime = 3600*6, # h in seconds
        runtime = 60*6, # h in minutes 
        memory = 50, # GB
    shell:
        """
        Rscript scripts/figures_inference_troubleshooting.R \
                    --experimental_pruned_path={input.experimental_pruned_path} \
                    --aracne_and_experimental_path={input.aracne_and_experimental_path} \
                    --mlr_and_experimental_path={input.mlr_and_experimental_path} \
                    --figs_dir={output}
        """