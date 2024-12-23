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
        # consider overfitting regulons
        os.path.join(RESULTS_DIR,"files","splicinglore_regulons-benchmarkable.tsv.gz")
        
        # make regulons
        # expand(os.path.join(RESULTS_DIR,"files","splicinglore_regulons-{omic_type}"), omic_type=EVENT_TYPES),
        

        
        
rule make_regulons:
    input:
        perts = os.path.join(PREP_DIR,"ground_truth_pert","SplicingLore","delta_psi-{omic_type}.tsv.gz"),
        regulators = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
    output:
        output_dir = directory(os.path.join(RESULTS_DIR,"files","splicinglore_regulons-{omic_type}"))
    params:
        omic_type = "{omic_type}"
    run:
        import os
        import pandas as pd
        import numpy as np
        
        # load
        perts = pd.read_table(input.perts)
        regulators = pd.read_table(input.regulators)
        dataset = "SplicingLore"
        value_name = "delta_psi"
        
        # keep only perturbations with splicing factors (we keep 62 out of 79)
        perts = perts.loc[perts["PERT_GENE"].isin(regulators["GENE"])].copy()
        
        # add info for regulator splicing factors
        perts = pd.merge(
            perts, 
            regulators[["GENE","ENSEMBL"]].rename(columns={"ENSEMBL":"PERT_ENSEMBL"}), 
            left_on="PERT_GENE", 
            right_on="GENE",
            how="left"
        )
        perts = perts.drop(columns=["GENE_y"]).rename(columns={"GENE_x":"GENE"})
        
        # make standard network
        perts["regulator"] = perts["PERT_ENSEMBL"]
        perts["target"] = perts["EVENT"]
        perts["likelihood"] = np.abs(perts["DeltaPSI"])
        perts["tfmode"] = (-1)*np.sign(perts["DeltaPSI"]) # DeltaPSI = KD - CTRL, we reverse for activity estimation

        # create metaexperiments with a perturbation in each splicing factor
        # try to put perturbations from the same project together
        cols_oi = ["PERT_ENSEMBL","cell_line","study_accession"]
        gene_study = perts[cols_oi].drop_duplicates().groupby(
            cols_oi
        ).size().reset_index().rename(columns={0:"n"}).sort_values(["n","cell_line"])
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

        
rule label_benchmark_regulons:
    input:
        metadata_ena = os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"),
        metadata_encorekd = os.path.join(PREP_DIR,"metadata","ENCOREKD.tsv.gz"),
        metadata_encoreko = os.path.join(PREP_DIR,"metadata","ENCOREKO.tsv.gz"),
        metadata_splicinglore = os.path.join(PREP_DIR,"metadata","SplicingLore.tsv.gz")
    output:
        metadata_splicinglore = os.path.join(RESULTS_DIR,"files","splicinglore_regulons-benchmarkable.tsv.gz")
    run:
        import pandas as pd
        
        # load
        metadata_ena = pd.read_table(input.metadata_ena)
        metadata_encorekd = pd.read_table(input.metadata_encorekd)
        metadata_encoreko = pd.read_table(input.metadata_encoreko)
        metadata_splicinglore = pd.read_table(input.metadata_splicinglore)
        
        # find samples in benchmark datasets in splicinglore (we should not consider them)
        metadata_splicinglore["in_ena"] = metadata_splicinglore["study_accession"].isin(metadata_ena["study_alias"])
        metadata_splicinglore["in_encorekd"] = metadata_splicinglore["study_accession"].isin(metadata_encorekd["experiment"])
        metadata_splicinglore["in_encoreko"] = metadata_splicinglore["study_accession"].isin(metadata_encoreko["experiment"])
        
        # add benchmark pert id
        id_cols = ["study_accession","cell_line","PERT_ENSEMBL","PERT_TYPE"]
        ## ENA
        metadata_ena["splicinglore_match"] = metadata_ena["study_alias"]
        metadata_ena = metadata_ena[["splicinglore_match"] + id_cols].drop_duplicates()
        metadata_ena["PERT_ID_BENCHMARK"] = metadata_ena[
            id_cols
        ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)        

        ## ENCOREKD
        metadata_encorekd["splicinglore_match"] = metadata_encorekd["experiment"]
        metadata_encorekd["study_accession"] = "ENCOREKD"
        metadata_encorekd = metadata_encorekd[["splicinglore_match"] + id_cols].drop_duplicates()
        metadata_encorekd["PERT_ID_BENCHMARK"] = metadata_encorekd[
                id_cols
            ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
        
        ## ENCOREKO
        metadata_encoreko["splicinglore_match"] = metadata_encoreko["experiment"]
        metadata_encoreko["study_accession"] = "ENCOREKO"
        metadata_encoreko = metadata_encoreko[["splicinglore_match"] + id_cols].drop_duplicates()
        metadata_encoreko["PERT_ID_BENCHMARK"] = metadata_encoreko[
                id_cols
            ].apply(lambda row: '___'.join(row.values.astype(str)), axis=1)
        
        ## merge benchmark metadata
        metadata_benchmarks = pd.concat([
            metadata_ena, metadata_encorekd, metadata_encoreko
        ])[["splicinglore_match","PERT_ID_BENCHMARK"]].drop_duplicates()
        
        ## join to splicinglore
        metadata_splicinglore = pd.merge(metadata_splicinglore, metadata_benchmarks, left_on="study_accession", right_on="splicinglore_match")
        
        # save
        metadata_splicinglore.to_csv(output.metadata_splicinglore, **SAVE_PARAMS)
        
        print("Done!")