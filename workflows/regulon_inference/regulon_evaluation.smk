import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
BIN_DIR = os.path.join(ROOT,"bin")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr"] + EVENT_TYPES

PERT_SPLICING_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{omic_type}.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{omic_type}.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{omic_type}.tsv.gz')
}

PERT_GENEXPR_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'log2_fold_change_tpm.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'log2_fold_change_tpm.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','log2_fold_change_tpm.tsv.gz')
}

EVAL_DATASETS = PERT_GENEXPR_FILES.keys()

PERT_FILES = {
    "EX": PERT_SPLICING_FILES,
    "genexpr": PERT_GENEXPR_FILES
}

METADATA_FILES = [
    os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz")
]

REGULON_SETS = [
    "aracne_regulons_development",
    "experimentally_derived_regulons_pruned"
]


##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        ## run
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{regulon_set}-{dataset}-{omic_type}.tsv.gz"),
               zip, dataset=EVAL_DATASETS, omic_type=OMIC_TYPES, regulon_set=REGULON_SETS),
        ## merge
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # make figures
        expand(os.path.join(RESULTS_DIR,"figures","regulon_evaluation-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule make_evaluation_labels:
    input:
        metadatas = METADATA_FILES
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"))
    run:
        import pandas as pd
        
        for f in input.metadatas:
            # load
            metadata = pd.read_table(f)
            
            os.makedirs(output.output_dir, exist_ok=True)
            if "ENCORE" in f:
                for cell_line in metadata["cell_line"].unique():
                    # make labels
                    labels = metadata.loc[
                        metadata["cell_line"]==cell_line, ["PERT_GENE","PERT_ENSEMBL"]
                    ].drop_duplicates().copy()
                    labels["PERT_ID"] = metadata["PERT_ENSEMBL"]
                    labels["PERT_TYPE"] = "KNOCKDOWN" if "KD" in os.path.basename(f) else "KNOCKOUT"
                    
                    dataset = "ENCOREKD" if "KD" in os.path.basename(f) else "ENCOREKO"
                    
                    # save
                    labels.dropna().to_csv(os.path.join(output.output_dir,"%s_%s.tsv.gz") % (dataset, cell_line), **SAVE_PARAMS)
                
            elif "ENASFS" in f:
                # prepare labels
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                labels = metadata[["PERT_ID","PERT_GENE","PERT_ENSEMBL","PERT_TYPE"]].drop_duplicates()

                # save
                labels.dropna().to_csv(os.path.join(output.output_dir,"ENASFS.tsv.gz"), **SAVE_PARAMS)
                
                
        print("Done!")
        
        
rule evaluate_regulons:
    input:
        signature = lambda wildcards: PERT_FILES[wildcards.omic_type][wildcards.dataset],
        regulons = os.path.join(RESULTS_DIR,"files","{regulon_set}-{omic_type}"),
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{regulon_set}-{dataset}-{omic_type}.tsv.gz")
    params:
        script_dir = BIN_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons} \
                    --eval_labels_file={input.eval_labels} \
                    --output_file={output}
        """
        
        
rule combine_evaluations:
    input:
        evaluations = [os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{regulon_set}-{dataset}-{omic_type}.tsv.gz").format(regulon_set=r, dataset=d, omic_type="{omic_type}") for r in REGULON_SETS for d in PERT_FILES.keys()]
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz")
    run:
        import pandas as pd
    
        evaluation = pd.concat([pd.read_table(f) for f in input.evaluations])
        
        evaluation.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        
    
rule figures_regulon_evaluation:
    input:
        evaluation = os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","merged-{omic_type}.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","regulon_evaluation-{omic_type}"))
    shell:
        """
        Rscript scripts/figures_regulon_evaluation.R \
                    --evaluation_file={input.evaluation} \
                    --figs_dir={output}
        """