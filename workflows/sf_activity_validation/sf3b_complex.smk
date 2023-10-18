"""
Author: Miquel Anglada Girotto
Contact: miquelangladagirotto [at] gmail [dot] com

Outline
-------
"""

import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
BIN_DIR = os.path.join(ROOT,"bin")
RESULTS_DIR = os.path.join(ROOT,"results","sf_activity_validation")
REGULONS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = EVENT_TYPES

##### RULES #####
rule all:
    input:
        # calculate signature
        expand(os.path.join(RESULTS_DIR,"files","signatures","sf3b_complex-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # compute viper SF activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","sf3b_complex-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # shortest paths to SF3B complex
        os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_sf3b_complex.tsv.gz'),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","validation_sf3b_complex")
        
        
rule compute_signatures:
    input:
        metadata_enasfs = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        splicing_enasfs = os.path.join(PREP_DIR,'event_psi','ENASFS-{omic_type}.tsv.gz'),
        metadata_sfdrugs = os.path.join(PREP_DIR,"metadata","sf_drugs.tsv.gz"),
        splicing_sfdrugs = os.path.join(PREP_DIR,"event_psi","sf_drugs-{omic_type}.tsv.gz"),
    output:
        metadatas = os.path.join(RESULTS_DIR,"files","metadata","sf3b_complex-{omic_type}.tsv.gz"),
        signatures = os.path.join(RESULTS_DIR,"files","signatures","sf3b_complex-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # ----- ENASFS -----
        # load
        metadata_enasfs = pd.read_table(input.metadata_enasfs)
        splicing_enasfs = pd.read_table(input.splicing_enasfs, index_col=0)
        
        # select perturbations
        perts_oi = ["E7107","MUT_PHF5A_Y36C","MUT_PHF5A_Y36C_AND_E7107","MUT_SUGP1_P636L","MUT_SF3B1_K700E"]
        metadata_enasfs = metadata_enasfs.loc[metadata_enasfs["condition"].isin(perts_oi)]
        metadata_enasfs = metadata_enasfs.loc[metadata_enasfs["study_accession"]!="PRJNA380104"]
        
        # delta PSI as the difference between conditions and the mean of the conditions
        signatures_enasfs = {}
        for sample_oi in metadata_enasfs["sampleID"]:
            # get the controls of the sample
            ctls = metadata_enasfs.loc[metadata_enasfs["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                psi_ctls = splicing_enasfs[ctls].mean(axis=1)
                
                # compute delta PSI
                dpsi = splicing_enasfs[sample_oi] - psi_ctls
                
                signatures_enasfs[sample_oi] = dpsi
                
                del dpsi, psi_ctls, ctls

        signatures_enasfs = pd.DataFrame(signatures_enasfs)
        
        # store
        metadata_enasfs["dataset"] = "ENASFS"
        metadatas.append(metadata_enasfs)
        signatures.append(signatures_enasfs)
        
        # delete
        del metadata_enasfs, signatures_enasfs
        
        # ----- SF drugs -----
        
        # load
        splicing_sfdrugs = pd.read_table(input.splicing_sfdrugs, index_col=0)
        metadata_sfdrugs = pd.read_table(input.metadata_sfdrugs)
        
        # subset
        perts_oi = [
            "DMSO","METHANOL",
            "ISOGINKGETIN","SPLICEOSTATIN_A","PLADIENOLIDE_B",
            "H3B-8800","H3B-8800_AND_SF3B1_K700E_MUTATION",
            "E7107","E7107_AND_SF3B1_K700E_MUTATION","E7107_AND_PHF5A_Y36C_MUTATION",
            "DMSO_AND_SF3B1_K700E_MUTATION","DMSO_AND_PHF5A_Y36C_MUTATION",
            "KD_CONTROL","KD_SF3B2"
        ]
        metadata_sfdrugs = metadata_sfdrugs.loc[metadata_sfdrugs["condition"].isin(perts_oi)]
        
        # subtract median DMSO to corresponding cell type and study
        signatures_sfdrugs = []
        for study in metadata_sfdrugs["study_accession"].unique():
            
            metadata_oi = metadata_sfdrugs.loc[metadata_sfdrugs["study_accession"]==study]
            signature = []
            for cell_line_oi in metadata_oi["cell_line_name"].unique():
                samples_cell_line = metadata_oi.loc[
                    metadata_oi["cell_line_name"] == cell_line_oi,
                    "sampleID"
                ]
                
                # DMSO by default
                samples_ctl = metadata_oi.loc[
                    (metadata_oi["cell_line_name"] == cell_line_oi) &\
                    (metadata_oi["condition"] == "DMSO"),
                    "sampleID"
                ]

                signature = (
                    splicing_sfdrugs[samples_cell_line] - \
                    splicing_sfdrugs[samples_ctl].median(axis=1).values.reshape(-1,1)
                )
                
                # spliceostatin A is diluted in methanol
                if any(metadata_oi["condition"]=="SPLICEOSTATIN_A"):
                    samples_ctl = metadata_oi.loc[
                        (metadata_oi["cell_line_name"] == cell_line_oi) &\
                        (metadata_oi["condition"] == "METHANOL"),
                        "sampleID"
                    ].tolist()
                    
                    samples_spliceostatin = metadata_oi.loc[
                        (metadata_oi["cell_line_name"] == cell_line_oi) &\
                        (metadata_oi["condition"] == "SPLICEOSTATIN_A"),
                        "sampleID"
                    ].tolist()
                    
                    samples_oi = samples_ctl + samples_spliceostatin
                    
                    signature.loc[:,samples_oi] = (
                        splicing_sfdrugs[samples_oi] - \
                        splicing_sfdrugs[samples_ctl].median(axis=1).values.reshape(-1,1)
                    )
                    
                signatures_sfdrugs.append(signature)
                
        signatures_sfdrugs = pd.concat(signatures_sfdrugs, axis=1)
        
        # store
        metadata_sfdrugs["dataset"] = "SF_DRUGS"
        metadatas.append(metadata_sfdrugs)
        signatures.append(signatures_sfdrugs)
        
        # delete
        del metadata_sfdrugs, signatures_sfdrugs
        
        # save
        pd.concat(metadatas).to_csv(output.metadatas, **SAVE_PARAMS)
        pd.concat(signatures, axis=1).reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")

        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","sf3b_complex-{omic_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","sf3b_complex-{omic_type}.tsv.gz")
    params:
        script_dir = BIN_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """
        

rule shortest_paths_stringdb:
    input:
        ppi = os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz'),
        sources = os.path.join(SUPPORT_DIR,"splicing_factors","sf3b_complex-symbol.txt"),
        targets = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-symbol.txt")
    output:
        os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_sf3b_complex.tsv.gz')
    threads: 16
    shell:
        """
        nice python scripts/ppi_path_lengths.py \
                    --ppi_file={input.ppi} \
                    --sources_file={input.sources} \
                    --targets_file={input.targets} \
                    --output_file={output} \
                    --n_jobs={threads}
        """ 
        
        
rule make_figures:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","sf3b_complex-EX.tsv.gz"),
        metadata = os.path.join(RESULTS_DIR,"files","metadata","sf3b_complex-EX.tsv.gz"),
        splicing_factors = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv"),
        shortest_paths = os.path.join(RESULTS_DIR,'files','ppi','shortest_path_lengths_to_sf3b_complex.tsv.gz')
    output:
        directory(os.path.join(RESULTS_DIR,"figures","validation_sf3b_complex"))
    shell:
        """
        
        Rscript scripts/figures_sf3b_complex.R \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --splicing_factors_file={input.splicing_factors} \
                    --shortest_paths_file={input.shortest_paths} \
                    --figs_dir={output}
        """