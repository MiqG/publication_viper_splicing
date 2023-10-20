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

##### RULES #####
rule all:
    input:
        # combine mlr and experimental regulons
        expand(os.path.join(RESULTS_DIR,"files","mlr_and_experimental_regulons-{omic_type}"), omic_type=OMIC_TYPES),

        
rule make_regulons:
    input:
        mlr_regulons = os.path.join(RESULTS_DIR,"files","mlr_regulons_development-{omic_type}"),
        experimental_regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","mlr_and_experimental_regulons-{omic_type}"))
    run:
        import numpy as np
        import pandas as pd
        
        # load mlr regulons
        mlr_regulons_files = [os.path.join(input.mlr_regulons,f) for f in os.listdir(input.mlr_regulons) if f.endswith(".tsv.gz")]
        assert len(mlr_regulons_files)==1
        mlr_regulons = pd.read_table(mlr_regulons_files[0]).rename(columns={"likelihood":"likelihood_mlr", "tfmode":"tfmode_mlr"})
        
        # change experimental likelihoods for mlr likelihoods
        experimental_regulons_files = [os.path.join(input.experimental_regulons,f) for f in os.listdir(input.experimental_regulons) if f.endswith(".tsv.gz")]
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in experimental_regulons_files:
            print(regulon_file)
            
            regulon = pd.read_table(regulon_file).rename(columns={"likelihood":"likelihood_experimental", "tfmode":"tfmode_experimental"})
            
            common_cols = list(set(regulon.columns).intersection(mlr_regulons.columns))
            regulon = pd.merge(regulon, mlr_regulons, how="left", on=common_cols)
            
            regulon["likelihood"] = regulon["likelihood_mlr"]
            regulon["tfmode"] = regulon["tfmode_experimental"]
            
            #idx_na = regulon["likelihood"].isnull()
            #regulon.loc[idx_na, "likelihood"] = regulon.loc[idx_na, "likelihood_experimental"]
            regulon = regulon.dropna()
            regulon = regulon.loc[np.abs(regulon["tfmode"]) > 0]
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")        