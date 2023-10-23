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
        
