import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

PERT_FILES = [
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"HepG2",'delta_psi-{event_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKD',"K562",'delta_psi-{event_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"HepG2",'delta_psi-{event_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENCOREKO',"K562",'delta_psi-{event_type}.tsv.gz'),
    os.path.join(PREP_DIR,'ground_truth_pert','ENASFS','delta_psi-{event_type}.tsv.gz')
]

EVENT_TYPES = ["EX"]

##### RULES #####
rule all:
    input:
        # make regulons
        expand(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{event_type}"), event_type=EVENT_TYPES),
        
        # prune regulons
        expand(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}"), event_type=EVENT_TYPES),
        
        
rule make_regulons:
    input:
        perts = PERT_FILES,
        metadata = os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"), # only for ENASFS
        regulators = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{event_type}"))
    run:
        import os
        import pandas as pd
        import numpy as np
        
        regulators = pd.read_table(input.regulators)

        # prep regulators
        regulators = regulators[["GENE","ENSEMBL"]]

        for f in input.perts:
            print("Loading %s..." % f)
            # load
            perts = pd.read_table(f, index_col=0)
            
            # get info
            if "ENCORE" in f:
                dataset = os.path.basename(os.path.dirname(os.path.dirname(f)))
                cell_line = os.path.basename(os.path.dirname(f))
            
            elif "ENASFS" in f:
                dataset = "ENASFS"
                metadata = pd.read_table(input.metadata)
                
            # prep perturbations
            perts.index.name = "EVENT"
            perts = perts.melt(
                ignore_index=False, var_name="ENSEMBL", value_name="delta_psi"
            ).dropna().reset_index().copy()

            # format
            perts["regulator"] = perts["ENSEMBL"]
            perts["target"] = perts["EVENT"]
            perts["likelihood"] = np.abs(perts["delta_psi"])
            perts["tfmode"] = (-1)*np.sign(perts["delta_psi"]) # they come from KD or KO, decrease activity
            
            os.makedirs(output.output_dir, exist_ok=True)
            if "ENCORE" in dataset:
                # subset
                perts = perts.loc[perts["ENSEMBL"].isin(regulators["ENSEMBL"])].copy()

                # add gene symbols
                perts = pd.merge(perts, regulators, on="ENSEMBL", how="left")

                # save
                output_file = os.path.join(output.output_dir,"%s-%s-delta_psi.tsv.gz") % (dataset, cell_line)
                print("Saving %s..." % output_file)
                perts.to_csv(output_file, **SAVE_PARAMS)
            
            elif dataset=="ENASFS":
                # correct ENSEMBL column
                X = perts["ENSEMBL"].str.split("___", expand=True)
                X.columns = ["study_accession","cell_line_name","ENSEMBL"]
                perts["PERT_ID"] = perts["ENSEMBL"]
                perts[["study_accession","cell_line_name","ENSEMBL"]] = X
                perts["regulator"] = perts["ENSEMBL"]
                
                # subset
                perts = perts.loc[perts["ENSEMBL"].isin(regulators["ENSEMBL"])].copy()
                
                # add gene symbols
                perts = pd.merge(perts, regulators, on="ENSEMBL", how="left")
                
                # add pert type
                metadata = metadata.loc[~metadata["PERT_ENSEMBL"].isnull()].copy()
                metadata["PERT_ID"] = metadata[
                    ["study_accession","cell_line_name","PERT_ENSEMBL"]
                ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
                perts = pd.merge(perts, metadata[["PERT_ID","PERT_TYPE"]].drop_duplicates(), on="PERT_ID", how="left")
                
                # subset pert types
                pert_types_oi = ["KNOCKDOWN","KNOCKOUT","OVEREXPRESSION"]
                perts = perts.loc[perts["PERT_TYPE"].isin(pert_types_oi)].copy()
                
                # correct overexpression sign as increase in activity
                idx = perts["PERT_TYPE"]=="OVEREXPRESSION"
                perts.loc[idx,"tfmode"] = -perts.loc[idx,"tfmode"]
                
                # create metacells with a perturbation in each splicing factor
                # try to put perturbations from the same project together
                cols_oi = ["ENSEMBL","cell_line_name","study_accession"]
                gene_study = perts[cols_oi].drop_duplicates().groupby(
                    cols_oi
                ).size().reset_index().rename(columns={0:"n"}).sort_values(["n","cell_line_name"])
                gene_study["index"] = np.arange(0,len(gene_study))
                
                metacells = {}
                it = 0
                while len(gene_study)>0:
                    to_keep = gene_study["index"].isin(gene_study.groupby('ENSEMBL')["index"].min())
                    metacells["metacell%s" % it] = gene_study.loc[to_keep].sort_values("ENSEMBL")
                    gene_study = gene_study.loc[~to_keep].sort_values("ENSEMBL").copy()
                    it = it + 1
                
                # save
                for metacell_oi in metacells.keys():
                    metacell = metacells[metacell_oi][cols_oi]
                    perts_oi = pd.merge(metacell, perts, on=cols_oi, how="left")

                    # save
                    if len(metacell)>1:
                        output_file = os.path.join(output.output_dir,"%s-%s-delta_psi.tsv.gz") % (dataset, metacell_oi)
                        print("Saving %s..." % output_file)
                        perts_oi.to_csv(output_file, **SAVE_PARAMS)

            
        print("Done!")

        
rule prune_regulons:
    input:
        regulons_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_raw-{event_type}")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}"))
    params:
        thresh = 15
    run:
        import os
        import pandas as pd
        
        regulon_files = [os.path.join(input.regulons_dir,f) for f in os.listdir(input.regulons_dir) if f.endswith(".tsv.gz")]
        thresh = params.thresh
        
        os.makedirs(output.output_dir, exist_ok=True)
        for regulon_file in regulon_files:
            print(regulon_file)
            
            regulon = pd.read_table(regulon_file)
            
            regulon = regulon.loc[regulon["likelihood"]>=thresh].copy()
            
            output_file = os.path.join(output.output_dir, os.path.basename(regulon_file))
            regulon.to_csv(output_file, **SAVE_PARAMS)
            
        print("Done!")