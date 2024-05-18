"""
Author: Miquel Anglada Girotto
Contact: miquelangladagirotto [at] gmail [dot] com

Outline
-------
- proliferation
    - MKI67 expression vs cancer-driver SFs in STN TCGA
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

# TCGA
metadata = pd.read_table(os.path.join(RAW_DIR,"TCGA","metadata","PANCAN.tsv.gz"))
metadata = metadata.groupby(
    ["sample_type_clean","cancer_type"]
).size().reset_index().rename(columns={0:"n"})
metadata = metadata.loc[metadata["n"]>=3]

DATASETS = ["CCLE"]

##### RULES ######
rule all:
    input:
        # compute signatures within
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"), zip, dataset=DATASETS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # compute viper SF activities        
        expand(os.path.join(RESULTS_DIR,"files","program_enrichment_scores","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # figures
        expand(os.path.join(RESULTS_DIR,"figures","proliferation-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule compute_signature_within:
    input:
        splicing_pt = os.path.join(PREP_DIR,"event_psi","{dataset}-{omic_type}.tsv.gz"),
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        splicing_pt = pd.read_table(input.splicing_pt, index_col=0)
        
        # subtract median PT
        signature = splicing_pt
        signature = signature - splicing_pt.median(axis=1).values.reshape(-1,1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
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
        protein_activity_ccle = os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-{omic_type}.tsv.gz"),
        genexpr_ccle = os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz"),
        metadata_ccle = os.path.join(PREP_DIR,"metadata","CCLE.tsv.gz"),
        driver_types = os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
    output:
        directory(os.path.join(RESULTS_DIR,"figures","proliferation-{omic_type}"))
    shell:
        """
        Rscript scripts/figures_proliferation.R \
                    --protein_activity_ccle_file={input.protein_activity_ccle} \
                    --genexpr_ccle_file={input.genexpr_ccle} \
                    --metadata_ccle_file={input.metadata_ccle} \
                    --driver_types_file={input.driver_types} \
                    --figs_dir={output}
        """