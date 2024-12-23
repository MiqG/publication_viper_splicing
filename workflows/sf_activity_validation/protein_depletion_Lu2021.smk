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
        expand(os.path.join(RESULTS_DIR,"files","signatures","protein_depletion-Lu2021-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        # compute viper SF activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","protein_depletion-Lu2021-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        # figures
        expand(os.path.join(RESULTS_DIR,"figures","protein_depletion-Lu2021-{event_type}"), event_type=EVENT_TYPES)
        
        
rule compute_signature:
    input:
        splicing = os.path.join(PREP_DIR,"event_psi","Lu2021-{event_type}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Lu2021.tsv.gz")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","protein_depletion-Lu2021-{event_type}.tsv.gz")
    run:
        import pandas as pd
        
        splicing = pd.read_table(input.splicing, index_col=0)
        metadata = pd.read_table(input.metadata)
        
        # subtract median DMSO to corresponding cell type
        signature = []
        for cell_line_oi in metadata["cell_line_name"].unique():
            samples_cell_line = metadata.loc[
                metadata["cell_line_name"] == cell_line_oi,
                "sampleID"
            ]
            samples_ctl = metadata.loc[
                (metadata["cell_line_name"] == cell_line_oi) &\
                (metadata["condition"] == "DMSO"),
                "sampleID"
            ]
            
            signature.append(
                splicing[samples_cell_line] - \
                splicing[samples_ctl].median(axis=1).values.reshape(-1,1)
            )
        signature = pd.concat(signature, axis=1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","protein_depletion-Lu2021-{event_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","protein_depletion-Lu2021-{event_type}.tsv.gz")
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
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","Lu2021.tsv.gz"),
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","protein_depletion-Lu2021-{event_type}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Lu2021.tsv.gz")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","protein_depletion-Lu2021-{event_type}"))
    shell:
        """
        Rscript scripts/figures_protein_depletion_Lu2021.R \
                    --genexpr_file={input.genexpr} \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --figs_dir={output}
        """