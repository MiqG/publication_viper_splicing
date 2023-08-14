"""
Author: Miquel Anglada Girotto
Contact: miquelangladagirotto [at] gmail [dot] com

Outline
-------
- study activities of cancer-driver splicing program 
    - during tumorigenesis in brain organoids
    - upon cancer-driver mutations
"""

import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
BIN_DIR = os.path.join(ROOT,"bin")
RESULTS_DIR = os.path.join(ROOT,"results","cancer_splicing_program")
REGULONS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = EVENT_TYPES

DATASETS = ["Bian2018","driver_mutations"]

###### RULES ######
rule all:
    input:
        # compute signatures
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        # compute protein activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # figures
        #os.path.join(RESULTS_DIR,"figures","tumor_initiation-{omic_type}")
        
        
rule compute_signature_within:
    input:
        splicing = os.path.join(PREP_DIR,"event_psi","Bian2018-{omic_type}.tsv.gz"),
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","Bian2018-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        splicing = pd.read_table(input.splicing, index_col=0)
        
        # subtract median PT
        signature = splicing
        signature = signature - splicing.median(axis=1).values.reshape(-1,1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
        
# rule compute_signatures_organoids:
#     input:
#         metadata = os.path.join(PREP_DIR,"metadata","Bian2018.tsv.gz"),
#         splicing = os.path.join(PREP_DIR,"event_psi","Bian2018-{omic_type}.tsv.gz")
#     output:
#         signatures = os.path.join(RESULTS_DIR,"files","signatures","Bian2018-{omic_type}.tsv.gz")
#     run:
#         import pandas as pd
        
#         # load
#         metadata = pd.read_table(input.metadata)
#         splicing = pd.read_table(input.splicing, index_col=0)
        
#         # delta PSI as the difference between conditions and the mean of the conditions
#         signatures = {}
#         for sample_oi in metadata["sampleID"]:
#             # get the controls of the sample
#             ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
#             # controls will be np.nan
#             if isinstance(ctls,str):
#                 ctls = ctls.split(",")
#                 psi_ctls = splicing[ctls].mean(axis=1)
                
#                 # compute delta PSI
#                 dpsi = splicing[sample_oi] - psi_ctls
                
#                 signatures[sample_oi] = dpsi
                
#                 del dpsi, psi_ctls, ctls
                
        
#         signatures = pd.DataFrame(signatures)
        
#         # save
#         signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
#         print("Done!")
        
        
rule compute_signatures_mutations:
    input:
        metadata_enasfs = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        splicing_enasfs = os.path.join(PREP_DIR,'event_psi','ENASFS-{omic_type}.tsv.gz'),
        metadata_sfptms = os.path.join(PREP_DIR,"metadata","sf_ptms.tsv.gz"),
        splicing_sfptms = os.path.join(PREP_DIR,"event_psi","sf_ptms-{omic_type}.tsv.gz"),
    output:
        metadatas = os.path.join(RESULTS_DIR,"files","metadata","driver_mutations-{omic_type}.tsv.gz"),
        signatures = os.path.join(RESULTS_DIR,"files","signatures","driver_mutations-{omic_type}.tsv.gz")
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
            metadata_enasfs["PERT_TYPE"]=="MUTATION"
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
        
        # ----- SF PTMs -----
        
        # load
        splicing_sfptms = pd.read_table(input.splicing_sfptms, index_col=0)
        metadata_sfptms = pd.read_table(input.metadata_sfptms)
        
        # subset
        ## mutation PHF5A
        metadata_sfptms = metadata_sfptms.loc[metadata_sfptms["study_accession"]=="PRJNA506256"]
        
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
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz")
    params:
        script_dir = BIN_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """

        
rule make_figures:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","Bian2018-{omic_type}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Bian2018.tsv.gz"),
    output:
        directory(os.path.join(RESULTS_DIR,"figures","tumor_initiation-{omic_type}"))
    shell:
        """
        Rscript scripts/figures_sf3b_complex.R \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --figs_dir={output}
        """