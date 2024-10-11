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
        # intersect postar3 and experimental regulons
        expand(os.path.join(RESULTS_DIR,"files","postar3_and_experimental_regulons-{omic_type}"), omic_type=OMIC_TYPES),
        
        # diff postar3 and experimental regulons
        expand(os.path.join(RESULTS_DIR,"files","experimental_without_postar3_regulons-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule make_postar3_and_experimental:
    input:
        clip_regulons = os.path.join(RESULTS_DIR,"files","postar3_clip_regulons-{omic_type}"),
        experimental_regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","postar3_and_experimental_regulons-{omic_type}"))
    run:
        import numpy as np
        import pandas as pd
        
        # get CLIP regulon files
        clip_regulons_files = [os.path.join(input.clip_regulons,f) for f in os.listdir(input.clip_regulons) if f.endswith(".tsv.gz")]
        assert len(clip_regulons_files)>1
        clip_regulons = pd.concat([pd.read_table(f) for f in clip_regulons_files])
        clip_interactions = clip_regulons[["regulator","target"]].drop_duplicates()
        
        # subset experimental interactions to only CLIP interactions
        experimental_regulons_files = [os.path.join(input.experimental_regulons,f) for f in os.listdir(input.experimental_regulons) if f.endswith(".tsv.gz")]
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in experimental_regulons_files:
            print(regulon_file)
            
            # load experimental regulon
            regulon = pd.read_table(regulon_file)
            
            # subset experimental regulon
            idx = (regulon["regulator"].isin(clip_interactions["regulator"]) & regulon["target"].isin(clip_interactions["target"]))
            regulon = regulon.loc[idx]
            
            if regulon.shape[0]>0:
                output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
                regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")
        
        
rule make_experimental_without_postar3:
    input:
        clip_regulons = os.path.join(RESULTS_DIR,"files","postar3_clip_regulons-{omic_type}"),
        experimental_regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimental_without_postar3_regulons-{omic_type}"))
    run:
        import numpy as np
        import pandas as pd
        
        # get CLIP regulon files
        clip_regulons_files = [os.path.join(input.clip_regulons,f) for f in os.listdir(input.clip_regulons) if f.endswith(".tsv.gz")]
        assert len(clip_regulons_files)>1
        clip_regulons = pd.concat([pd.read_table(f) for f in clip_regulons_files])
        clip_interactions = clip_regulons[["regulator","target"]].drop_duplicates()
        
        # subset experimental interactions to only CLIP interactions
        experimental_regulons_files = [os.path.join(input.experimental_regulons,f) for f in os.listdir(input.experimental_regulons) if f.endswith(".tsv.gz")]
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in experimental_regulons_files:
            print(regulon_file)
            
            # load experimental regulon
            regulon = pd.read_table(regulon_file)
            
            # subset experimental regulon for interactions not found
            idx = (
                (~regulon["regulator"].isin(clip_interactions["regulator"])) & 
                (~regulon["target"].isin(clip_interactions["target"]))
            )
            regulon = regulon.loc[idx]
            
            if regulon.shape[0]>0:
                output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
                regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")
