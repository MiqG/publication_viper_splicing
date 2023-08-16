"""
Author: Miquel Anglada Girotto
Contact: miquelangladagirotto [at] gmail [dot] com

Outline
-------
- proliferation
    - MKI67 expression vs cell line doubling vs cancer-driver SFs in CCLE
    - MKI67 expression vs cancer-driver SFs in STN TCGA
- tumor evasion
    - markers vs cancer-driver SFs
- metastasis
    - MetMap scores vs cancer-driver SFs in CCLE
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

CANCER_TYPES_STN = metadata.loc[metadata["sample_type_clean"]=="SolidTissueNormal","cancer_type"].values

DATASETS = ["PANCAN-SolidTissueNormal","CCLE"]

##### RULES ######
rule all:
    input:
        # make PANCAN for STN
        expand(os.path.join(PREP_DIR,"event_psi","PANCAN-{sample}-{omic_type}.tsv.gz"), zip, sample=["SolidTissueNormal"], omic_type=OMIC_TYPES),
        expand(os.path.join(PREP_DIR,"genexpr_tpm","PANCAN-{sample}.tsv.gz"), zip, sample=["SolidTissueNormal"]),

        # compute signatures within
        expand(os.path.join(RESULTS_DIR,"files","signatures","{dataset}-{omic_type}.tsv.gz"), zip, dataset=DATASETS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # compute viper SF activities        
        expand(os.path.join(RESULTS_DIR,"files","program_enrichment_scores","{dataset}-{omic_type}.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # compute cancer-driver NES
        
        
        # figures
        #expand(os.path.join(RESULTS_DIR,"figures","cancer_hallmarks-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule merge_event_psi_by_sample_type:
    input:
        psi = [os.path.join(PREP_DIR,"event_psi","{cancer}-{sample}.tsv.gz").format(cancer=c, sample="{sample}") for c in CANCER_TYPES_STN]
    output:
        psi = os.path.join(PREP_DIR,"event_psi","PANCAN-{sample}.tsv.gz")
    run:
        import pandas as pd
        
        psi = pd.concat([pd.read_table(f, index_col=0) for f in input.psi], axis=1)
        psi.reset_index().to_csv(output.psi, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule merge_genexpr_tpm_by_sample_type:
    input:
        genexpr = [os.path.join(PREP_DIR,"genexpr_tpm","{cancer}-{sample}.tsv.gz").format(cancer=c, sample="{sample}") for c in CANCER_TYPES_STN]
    output:
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","PANCAN-{sample}.tsv.gz")
    run:
        import pandas as pd
        
        genexpr = pd.concat([pd.read_table(f, index_col=0) for f in input.genexpr], axis=1)
        genexpr.reset_index().to_csv(output.genexpr, **SAVE_PARAMS)
        
        print("Done!")
        
        
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
        
        
rule compute_program_enrichment_score:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{dataset}-{omic_type}.tsv.gz"),
        gene_sets = os.path.join("scripts","driver_types-gene_sets.tsv")
    output:
        os.path.join(RESULTS_DIR,"files","program_enrichment_scores","{dataset}-{omic_type}.tsv.gz")
    params:
        script_dir = BIN_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_gene_sets_enrichment_score.R \
                    --omic_table_file={input.protein_activity} \
                    --gene_sets_file={input.gene_sets} \
                    --output_file={output}
        """
    
    

rule make_figures:
    input:
        protein_activity_tcga_stn = os.path.join(RESULTS_DIR,"files","protein_activity","PANCAN-SolidTissueNormal-{omic_type}.tsv.gz"),
        genexpr_tcga_stn = os.path.join(PREP_DIR,"genexpr_tpm","PANCAN-SolidTissueNormal.tsv.gz"),
        metadata_tcga = os.path.join(PREP_DIR,"metadata","PANCAN.tsv.gz"),
        protein_activity_ccle = os.path.join(RESULTS_DIR,"files","protein_activity","CCLE-{omic_type}.tsv.gz"),
        genexpr_ccle = os.path.join(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz"),
        metadata_ccle = os.path.join(PREP_DIR,"metadata","CCLE.tsv.gz"),
        doublings_ccle = os.path.join(PREP_DIR,"doubling_times","CCLE.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","cancer_hallmarks-{omic_type}"))
    shell:
        """
        Rscript scripts/figures_sf3b_complex.R \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --figs_dir={output}
        """