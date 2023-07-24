import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
BIN_DIR = os.path.join(ROOT,"bin")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]

PERT_FILES = {
    "ENCOREKD_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{event_type}.tsv.gz'),
    "ENCOREKD_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{event_type}.tsv.gz'),
    "ENCOREKO_HepG2": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{event_type}.tsv.gz'),
    "ENCOREKO_K562": os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{event_type}.tsv.gz'),
    "ENASFS": os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{event_type}.tsv.gz')
}

METADATA_FILES = [
    os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
    os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz")
]

##### RULES #####
rule all:
    input:
        # prepare evaluation labels
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels"),
        
        # evaluate regulons
        expand(os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{event_type}.tsv.gz"), dataset=PERT_FILES.keys(), event_type=EVENT_TYPES)
        
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
        signature = lambda wildcards: PERT_FILES[wildcards.dataset],
        regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}"),
        eval_labels = os.path.join(RESULTS_DIR,"files","regulon_evaluation_labels","{dataset}.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","regulon_evaluation_scores","{dataset}-{event_type}.tsv.gz")
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