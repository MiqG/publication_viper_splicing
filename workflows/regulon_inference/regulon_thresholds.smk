import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = EVENT_TYPES

THRESH_DPSI = [5,10,15,20,25,30,35,40,45]

##### RULES #####
rule all:
    input:
        # prune regulons with different thresholds
        expand(os.path.join(RESULTS_DIR,"files","dPSIthresh{thresh}_experimentally_derived_regulons_pruned-{omic_type}"), omic_type=EVENT_TYPES, thresh=THRESH_DPSI),
        
        # get properties
        expand(os.path.join(RESULTS_DIR,"files","regulon_properties","dPSIthresh-regulators_per_target-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","regulon_properties","dPSIthresh-targets_per_regulator-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES)
        
        
rule prune_regulons:
    input:
        regulons_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","dPSIthresh{thresh}_experimentally_derived_regulons_pruned-{omic_type}"))
    params:
        thresh = "{thresh}"
    run:
        import os
        import pandas as pd
        
        regulon_files = [os.path.join(input.regulons_dir,f) for f in os.listdir(input.regulons_dir) if f.endswith(".tsv.gz")]
        thresh = int(params.thresh)
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in regulon_files:
            print(regulon_file)
            
            regulon = pd.read_table(regulon_file)
            
            regulon = regulon.loc[regulon["likelihood"]>=thresh].copy()
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")            
     
    
rule properties_regulons:
    input:
        regulons_dirs = [os.path.join(RESULTS_DIR,"files","dPSIthresh{thresh}_experimentally_derived_regulons_pruned-{omic_type}").format(thresh=thresh, omic_type="{omic_type}") for thresh in THRESH_DPSI]
    output:
        regulators_per_target = os.path.join(RESULTS_DIR,"files","regulon_properties","dPSIthresh-regulators_per_target-{omic_type}.tsv.gz"),
        targets_per_regulator = os.path.join(RESULTS_DIR,"files","regulon_properties","dPSIthresh-targets_per_regulator-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        regulons_dirs = input.regulons_dirs
        regulators_per_targets = []
        targets_per_regulators = []
        for regulons_dir in regulons_dirs:
            regulons_files = [os.path.join(regulons_dir,f) for f in os.listdir(regulons_dir) if f.endswith(".tsv.gz")]
            regulons = pd.concat([pd.read_table(f) for f in regulons_files])        
            
            edges = regulons[["regulator","target"]].drop_duplicates()
            targets_per_regulator = edges.value_counts("regulator").reset_index().rename(columns={0:"n_targets"})
            regulators_per_target = edges.value_counts("target").reset_index().rename(columns={0:"n_regulators"})
            targets_per_regulator["regulon_set_id"] = os.path.basename(regulons_dir)
            regulators_per_target["regulon_set_id"] = os.path.basename(regulons_dir)
            
            regulators_per_targets.append(regulators_per_target)
            targets_per_regulators.append(targets_per_regulator)
            
        regulators_per_targets = pd.concat(regulators_per_targets)
        targets_per_regulators = pd.concat(targets_per_regulators)
            
        # save
        regulators_per_targets.to_csv(output.regulators_per_target, **SAVE_PARAMS)
        targets_per_regulators.to_csv(output.targets_per_regulator, **SAVE_PARAMS)
        
        print("Done!")