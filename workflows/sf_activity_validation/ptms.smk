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
        expand(os.path.join(RESULTS_DIR,"files","signatures","ptms-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # compute viper SF activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","ptms-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # figures
        #os.path.join(RESULTS_DIR,"figures","validation_ptms")
        
        
rule compute_signatures:
    input:
        metadata_enasfs = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        splicing_enasfs = os.path.join(PREP_DIR,'event_psi','ENASFS-{omic_type}.tsv.gz'),
        metadata_Lu2021 = os.path.join(PREP_DIR,"metadata","Lu2021.tsv.gz"),
        splicing_Lu2021 = os.path.join(PREP_DIR,"event_psi","Lu2021-{omic_type}.tsv.gz"),
        metadata_sfdrugs = os.path.join(PREP_DIR,"metadata","sf_drugs.tsv.gz"),
        splicing_sfdrugs = os.path.join(PREP_DIR,"event_psi","sf_drugs-{omic_type}.tsv.gz"),
        metadata_sfptms = os.path.join(PREP_DIR,"metadata","sf_ptms.tsv.gz"),
        splicing_sfptms = os.path.join(PREP_DIR,"event_psi","sf_ptms-{omic_type}.tsv.gz"),
    output:
        metadatas = os.path.join(RESULTS_DIR,"files","metadata","ptms-{omic_type}.tsv.gz"),
        signatures = os.path.join(RESULTS_DIR,"files","signatures","ptms-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        metadatas = []
        signatures = []
        
        # ----- ENASFS -----
        # load
        metadata_enasfs = pd.read_table(input.metadata_enasfs)
        splicing_enasfs = pd.read_table(input.splicing_enasfs, index_col=0)
        
        # select perturbations
        idx = (
            metadata_enasfs["condition"].isin(["T3","KH-CB19"]) | 
            metadata_enasfs["condition"].str.contains("CLK")
        )
        metadata_enasfs = metadata_enasfs.loc[idx]
        
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
        
        # ----- Lu2021 -----
        splicing_Lu2021 = pd.read_table(input.splicing_Lu2021, index_col=0)
        metadata_Lu2021 = pd.read_table(input.metadata_Lu2021)
        
        # subset
        metadata_Lu2021 = metadata_Lu2021.loc[metadata_Lu2021["condition"].isin(["DMSO","MS023"])]
        
        # subtract median DMSO to corresponding cell type
        signatures_Lu2021 = []
        for cell_line_oi in metadata_Lu2021["cell_line_name"].unique():
            samples_cell_line = metadata_Lu2021.loc[
                metadata_Lu2021["cell_line_name"] == cell_line_oi,
                "sampleID"
            ]
            samples_ctl = metadata_Lu2021.loc[
                (metadata_Lu2021["cell_line_name"] == cell_line_oi) &\
                (metadata_Lu2021["condition"] == "DMSO"),
                "sampleID"
            ]
            
            signatures_Lu2021.append(
                splicing_Lu2021[samples_cell_line] - \
                splicing_Lu2021[samples_ctl].median(axis=1).values.reshape(-1,1)
            )
        signatures_Lu2021 = pd.concat(signatures_Lu2021, axis=1)        
        
        # store
        metadata_Lu2021["dataset"] = "Lu2021"
        metadatas.append(metadata_Lu2021)
        signatures.append(signatures_Lu2021)
        
        # delete
        del metadata_Lu2021, signatures_Lu2021
        
        # ----- SF drugs -----
        
        # load
        splicing_sfdrugs = pd.read_table(input.splicing_sfdrugs, index_col=0)
        metadata_sfdrugs = pd.read_table(input.metadata_sfdrugs)
        
        # subset
        perts_oi = ["DMSO","GSK3326595","PALBOCICLIB"]
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
        
        # ----- SF PTMs -----
        
        # load
        splicing_sfptms = pd.read_table(input.splicing_sfptms, index_col=0)
        metadata_sfptms = pd.read_table(input.metadata_sfptms)
        
        # subset
        metadata_sfptms = metadata_sfptms.loc[~metadata_sfptms["condition"].isnull()]
        
        # subtract median DMSO to corresponding cell type and study
        signatures_sfptms = []
        for study in metadata_sfptms["study_accession"].unique():
            
            metadata_oi = metadata_sfptms.loc[metadata_sfptms["study_accession"]==study]
            signature = []
            for cell_line_oi in metadata_oi["cell_line_name"].unique():
                samples_cell_line = metadata_oi.loc[
                    metadata_oi["cell_line_name"] == cell_line_oi,
                    "sampleID"
                ]
                
                # DMSO by default
                samples_ctl = metadata_oi.loc[
                    (metadata_oi["cell_line_name"] == cell_line_oi) &\
                    (metadata_oi["condition"].isin(["DMSO","CONTROL"])),
                    "sampleID"
                ]

                signature = (
                    splicing_sfptms[samples_cell_line] - \
                    splicing_sfptms[samples_ctl].median(axis=1).values.reshape(-1,1)
                )
                
                signatures_sfptms.append(signature)
                
        signatures_sfptms = pd.concat(signatures_sfptms, axis=1)
        
        # store
        metadata_sfptms["dataset"] = "SF_PTMs"
        metadatas.append(metadata_sfptms)
        signatures.append(signatures_sfptms)
        
        # delete
        del metadata_sfptms, signatures_sfptms
        
        # save
        pd.concat(metadatas).to_csv(output.metadatas, **SAVE_PARAMS)
        pd.concat(signatures, axis=1).reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")

        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","ptms-{omic_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","ptms-{omic_type}.tsv.gz")
    params:
        script_dir = BIN_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """
        
        
# rule make_figures:
#     input:
#         genexpr = os.path.join(PREP_DIR,"genexpr_tpm","sf_ptms.tsv.gz"),
#         protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","sf_ptms.tsv.gz"),
#         metadata = os.path.join(PREP_DIR,"metadata","sf_ptms.tsv.gz"),
#         annotation = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
#     output:
#         directory(os.path.join(RESULTS_DIR,"figures","validation_ptms"))
#     shell:
#         """
#         Rscript scripts/figures_ptms.R \
#                     --genexpr_file={input.genexpr} \
#                     --protein_activity_file={input.protein_activity} \
#                     --metadata_file={input.metadata} \
#                     --annotation_file={input.annotation} \
#                     --figs_dir={output}
#         """