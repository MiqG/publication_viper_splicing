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
        # make regulons
        expand(os.path.join(RESULTS_DIR,"files","postar3_clip_regulons-{omic_type}"), omic_type=EVENT_TYPES),
        
        
rule make_regulons:
    input:
        mapped_clip = os.path.join(PREP_DIR,"clip_peaks_mapped","POSTAR3.tsv.gz"),
        regulators = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","postar3_clip_regulons-{omic_type}"))
    params:
        omic_type = "{omic_type}"
    run:
        import os
        import pandas as pd
        import numpy as np
        
        # load
        mapped_clip = pd.read_table(input.mapped_clip)
        regulators = pd.read_table(input.regulators)
        dataset = "POSTAR3"
        value_name = "clip_peak"
        
        # keep only splicing factors (we keep 97 out of 221)
        perts = mapped_clip.loc[mapped_clip["RBP"].isin(regulators["GENE"])].copy()
        
        # add SF info
        perts["PERT_GENE"] = perts["RBP"]
        perts = pd.merge(
            perts, 
            regulators[["GENE","ENSEMBL"]].rename(columns={"ENSEMBL":"PERT_ENSEMBL"}), 
            left_on="PERT_GENE", 
            right_on="GENE", 
            how="left"
        )
        
        # make standard network
        perts["regulator"] = perts["PERT_ENSEMBL"]
        perts["target"] = perts["EVENT"]
        perts["likelihood"] = 1
        perts["tfmode"] = 1
        
        # add pert id
        perts["PERT_ID"] = perts[
            ["experiment_id","experiment_model","PERT_ENSEMBL","experiment_method"]
        ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
        
        # create metaexperiments with a perturbation in each splicing factor
        # try to put perturbations from the same project together
        cols_oi = ["PERT_ENSEMBL","experiment_model","experiment_id"]
        gene_study = perts[cols_oi].drop_duplicates().groupby(
            cols_oi
        ).size().reset_index().rename(columns={0:"n"}).sort_values(["n","experiment_model"])
        gene_study["index"] = np.arange(0,len(gene_study))

        metaexperiments = {}
        it = 0
        while len(gene_study)>0:
            to_keep = gene_study["index"].isin(gene_study.groupby('PERT_ENSEMBL')["index"].min())
            metaexperiments["metaexperiment%s" % it] = gene_study.loc[to_keep].sort_values("PERT_ENSEMBL")
            gene_study = gene_study.loc[~to_keep].sort_values("PERT_ENSEMBL").copy()
            it = it + 1

        # save
        os.makedirs(output.output_dir, exist_ok=True)
        for metaexperiment_oi in metaexperiments.keys():
            metaexperiment = metaexperiments[metaexperiment_oi][cols_oi]
            perts_oi = pd.merge(metaexperiment, perts, on=cols_oi, how="left")

            # save
            if len(metaexperiment)>1:
                output_file = os.path.join(output.output_dir,"%s-%s-%s.tsv.gz") % (dataset, metaexperiment_oi, value_name)
                print("Saving %s..." % output_file)
                perts_oi.to_csv(output_file, **SAVE_PARAMS)

            
        print("Done!")

