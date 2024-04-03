import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr"] + EVENT_TYPES

TOP_N = [100, 90, 80, 70, 60, 50, 40]

##### RULES #####
rule all:
    input:
        # derive topN regulons
        expand(os.path.join(RESULTS_DIR,"files","top{N}_experimentally_derived_regulons_pruned-{omic_type}"), omic_type=OMIC_TYPES, N=TOP_N),
        
        # get properties
        expand(os.path.join(RESULTS_DIR,"files","regulon_properties","regulators_per_target-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","regulon_properties","targets_per_regulator-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES)
        
        
rule make_topn_regulons:
    input:
        experimental_regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","top{N}_experimentally_derived_regulons_pruned-{omic_type}"))
    params:
        top_n = "{N}"
    run:
        import pandas as pd
        
        experimental_regulons_files = [os.path.join(input.experimental_regulons,f) for f in os.listdir(input.experimental_regulons) if f.endswith(".tsv.gz")]
        top_n = int(params.top_n)
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in experimental_regulons_files:
            print(regulon_file)
            
            regulon = pd.read_table(regulon_file)
            
            regulon = regulon.sort_values("likelihood", ascending=False).groupby("regulator").head(top_n)
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")              
     
    
rule properties_regulons:
    input:
        regulons_dirs = [os.path.join(RESULTS_DIR,"files","top{N}_experimentally_derived_regulons_pruned-{omic_type}").format(N=n, omic_type="{omic_type}") for n in TOP_N] + [os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")]
    output:
        regulators_per_target = os.path.join(RESULTS_DIR,"files","regulon_properties","regulators_per_target-{omic_type}.tsv.gz"),
        targets_per_regulator = os.path.join(RESULTS_DIR,"files","regulon_properties","targets_per_regulator-{omic_type}.tsv.gz")
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