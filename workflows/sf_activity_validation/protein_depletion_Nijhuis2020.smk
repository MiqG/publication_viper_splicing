"""
Author: Miquel Anglada Girotto
Contact: miquelangladagirotto [at] gmail [dot] com
"""

import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
SRC_DIR = os.path.join(ROOT,"src")
RESULTS_DIR = os.path.join(ROOT,"results","sf_activity_validation")
REGULONS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]

##### RULES #####
rule all:
    input:
        # calculate signature
        expand(os.path.join(RESULTS_DIR,"files","signatures","protein_depletion-Nijhuis2020-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        # compute viper SF activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","protein_depletion-Nijhuis2020-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        # figures
        expand(os.path.join(RESULTS_DIR,"figures","protein_depletion-Nijhuis2020-{event_type}"), event_type=EVENT_TYPES)
        
        
rule compute_signature:
    input:
        splicing = os.path.join(PREP_DIR,"event_psi","Nijhuis2020-{event_type}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Nijhuis2020.tsv.gz")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","protein_depletion-Nijhuis2020-{event_type}.tsv.gz")
    run:
        import pandas as pd
        
        splicing = pd.read_table(input.splicing, index_col=0)
        metadata = pd.read_table(input.metadata)
        
        # delta PSI as the difference between conditions and the mean of the conditions
        signature = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                psi_ctls = splicing[ctls].median(axis=1)
                
                # compute delta PSI
                dpsi = splicing[sample_oi] - psi_ctls
                
                signature[sample_oi] = dpsi
                
                del dpsi, psi_ctls, ctls
                
        
        signature = pd.DataFrame(signature)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","protein_depletion-Nijhuis2020-{event_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","protein_depletion-Nijhuis2020-{event_type}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """


rule figures_protein_depletion:
    input:
        proteomics = os.path.join(RAW_DIR,"articles","Nijhuis2020","supplementary_data","proteomics_lfq_intensity.tsv.gz"),
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","Nijhuis2020.tsv.gz"),
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","protein_depletion-Nijhuis2020-{event_type}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Nijhuis2020.tsv.gz"),
        gene_info = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","protein_depletion-Nijhuis2020-{event_type}"))
    shell:
        """
        Rscript scripts/figures_protein_depletion_Nijhuis2020.R \
                    --proteomics_file={input.proteomics} \
                    --genexpr_file={input.genexpr} \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --gene_info_file={input.gene_info} \
                    --figs_dir={output}
        """