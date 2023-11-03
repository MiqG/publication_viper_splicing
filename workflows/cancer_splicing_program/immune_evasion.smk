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

CONDITIONS = ["PRE_PD1","PRE_CombPD1_CTLA4","ON_PD1","ON_CombPD1_CTLA4"]
CONDITIONS = ["PRE","ON"]

##### RULES ######
rule all:
    input:
        # make datasets
        expand(os.path.join(PREP_DIR,'event_psi','Riaz2017-{condition}-{omic_type}.tsv.gz'), condition=CONDITIONS, omic_type=OMIC_TYPES),
        expand(os.path.join(PREP_DIR,'genexpr_tpm','Riaz2017-{condition}.tsv.gz'), condition=CONDITIONS),
        
        # compute signatures within
        expand(os.path.join(RESULTS_DIR,"files","signatures","Riaz2017-{condition}-{omic_type}.tsv.gz"), condition=CONDITIONS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","Riaz2017-{condition}-{omic_type}.tsv.gz"), condition=CONDITIONS, omic_type=OMIC_TYPES),
        #expand(os.path.join(RESULTS_DIR,"files","protein_activity","Riaz2017-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # figures
        
        
rule split_psi_by_condition:
    input:
        metadata = os.path.join(PREP_DIR,"metadata","Riaz2017.tsv.gz"),
        psi = os.path.join(PREP_DIR,'event_psi','Riaz2017-{omic_type}.tsv.gz'),
    output:
        psi = os.path.join(PREP_DIR,'event_psi','Riaz2017-{condition}-{omic_type}.tsv.gz')
    params:
        condition = "{condition}"
    run:
        import pandas as pd
        
        metadata = pd.read_table(input.metadata)
        psi = pd.read_table(input.psi, index_col=0)
        
        idx = metadata["treatment_status"] == params.condition
        samples_oi = list(set(metadata.loc[idx,"sampleID"]).intersection(psi.columns))
        
        psi[samples_oi].reset_index().to_csv(output.psi, **SAVE_PARAMS)
        
        print("Done!")
        

rule split_genexpr_by_condition:
    input:
        metadata = os.path.join(PREP_DIR,"metadata","Riaz2017.tsv.gz"),
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','Riaz2017.tsv.gz')
    output:
        genexpr = os.path.join(PREP_DIR,'genexpr_tpm','Riaz2017-{condition}.tsv.gz')
    params:
        condition = "{condition}"
    run:
        import pandas as pd
        
        metadata = pd.read_table(input.metadata)
        genexpr = pd.read_table(input.genexpr, index_col=0)
        
        idx = metadata["treatment_status"] == params.condition
        samples_oi = list(set(metadata.loc[idx,"sampleID"]).intersection(genexpr.columns))
        
        genexpr[samples_oi].reset_index().to_csv(output.genexpr, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_signature_within:
    input:
        splicing = os.path.join(PREP_DIR,'event_psi','Riaz2017-{condition}-{omic_type}.tsv.gz')
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","Riaz2017-{condition}-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        splicing = pd.read_table(input.splicing, index_col=0)
        
        # subtract median PT
        signature = splicing
        signature = signature - splicing.median(axis=1).values.reshape(-1,1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")        


rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","Riaz2017-{condition}-{omic_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","Riaz2017-{condition}-{omic_type}.tsv.gz")
    params:
        script_dir = BIN_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """
        
        
# rule combine_protein_activities:
#     input:
#         signatures = [os.path.join(RESULTS_DIR,"files","protein_activity","Riaz2017-{condition}-{omic_type}.tsv.gz").format(condition=c, omic_type="{omic_type}") for c in CONDITIONS]
#     output:
#         signatures = os.path.join(RESULTS_DIR,"files","protein_activity","Riaz2017-{omic_type}.tsv.gz")
#     run:
#         import pandas as pd
        
#         signatures = pd.concat([pd.read_table(signature, index_col=0) for signature in input.signatures], axis=1)
#         signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
#         print("Done!")
