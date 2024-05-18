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

##### RULES #####
rule all:
    input:
        # combine aracne and experimental regulons
        expand(os.path.join(RESULTS_DIR,"files","aracne_and_experimental_regulons-{omic_type}"), omic_type=OMIC_TYPES),
        
        # combine mlr and experimental regulons
        expand(os.path.join(RESULTS_DIR,"files","mlr_and_experimental_regulons-{omic_type}"), omic_type=OMIC_TYPES),
        
        # combine aracne and mlr regulons
        expand(os.path.join(RESULTS_DIR,"files","aracne_and_mlr_regulons-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule make_regulons_aracne_and_experimental:
    input:
        aracne_regulons = os.path.join(RESULTS_DIR,"files","aracne_regulons_development-{omic_type}"),
        experimental_regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","aracne_and_experimental_regulons-{omic_type}"))
    run:
        import numpy as np
        import pandas as pd
        
        # load aracne regulons
        aracne_regulons_files = [os.path.join(input.aracne_regulons,f) for f in os.listdir(input.aracne_regulons) if f.endswith(".tsv.gz")]
        assert len(aracne_regulons_files)==1
        aracne_regulons = pd.read_table(aracne_regulons_files[0]).rename(columns={"likelihood":"likelihood_aracne", "tfmode":"tfmode_aracne"})
        
        # change experimental likelihoods for aracne likelihoods
        experimental_regulons_files = [os.path.join(input.experimental_regulons,f) for f in os.listdir(input.experimental_regulons) if f.endswith(".tsv.gz")]
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in experimental_regulons_files:
            print(regulon_file)
            
            regulon = pd.read_table(regulon_file).rename(columns={"likelihood":"likelihood_experimental", "tfmode":"tfmode_experimental"})
            
            common_cols = list(set(regulon.columns).intersection(aracne_regulons.columns))
            regulon = pd.merge(regulon, aracne_regulons, how="left", on=common_cols)
            
            regulon["likelihood"] = regulon["likelihood_aracne"]
            regulon["tfmode"] = regulon["tfmode_experimental"]
            
            idx = (np.abs(regulon["tfmode"]) > 0) & (~regulon["likelihood"].isnull())
            regulon = regulon.loc[idx]
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")              
        
        
        
rule make_regulons_mlr_and_experimental:
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
            
            idx = (np.abs(regulon["tfmode"]) > 0) & (~regulon["likelihood"].isnull())
            regulon = regulon.loc[idx]
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")        
        
        
rule make_regulons_aracne_and_mlr:
    input:
        aracne_regulons = os.path.join(RESULTS_DIR,"files","aracne_regulons_development-{omic_type}"),
        mlr_regulons = os.path.join(RESULTS_DIR,"files","mlr_regulons_development-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","aracne_and_mlr_regulons-{omic_type}"))
    run:
        import numpy as np
        import pandas as pd
        
        # load aracne regulons
        aracne_regulons_files = [os.path.join(input.aracne_regulons,f) for f in os.listdir(input.aracne_regulons) if f.endswith(".tsv.gz")]
        assert len(aracne_regulons_files)==1
        aracne_regulons = pd.read_table(aracne_regulons_files[0]).rename(columns={"likelihood":"likelihood_aracne", "tfmode":"tfmode_aracne"})

        # load mlr regulons
        mlr_regulons_files = [os.path.join(input.mlr_regulons,f) for f in os.listdir(input.mlr_regulons) if f.endswith(".tsv.gz")]
        assert len(mlr_regulons_files)==1
        mlr_regulons = pd.read_table(mlr_regulons_files[0]).rename(columns={"likelihood":"likelihood_mlr", "tfmode":"tfmode_mlr"})
        
        # combine regulons: keep ARACNe likelihoods, keep ARACNe negative MoR and MLR positive MoR
        os.makedirs(output.output_dir, exist_ok=True)
        
        common_cols = list(set(aracne_regulons.columns).intersection(mlr_regulons.columns))
        regulon = pd.merge(aracne_regulons, mlr_regulons, how="left", on=common_cols)
        
        regulon["likelihood"] = regulon["likelihood_aracne"]
        regulon["tfmode"] = np.nan
        regulon.loc[regulon["tfmode_aracne"]<0, "tfmode"] = -1 # negative
        regulon.loc[regulon["tfmode_mlr"]>0, "tfmode"] = 1 # positive
        regulon = regulon.loc[~regulon["tfmode"].isnull()]
        
        output_file = os.path.join(output.output_dir, os.path.basename(aracne_regulons_files[0]))
        regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")     