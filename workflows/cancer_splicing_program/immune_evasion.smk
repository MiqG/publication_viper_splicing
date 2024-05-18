import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","cancer_splicing_program")
REGULONS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = EVENT_TYPES

CONDITIONS = ["PRE","ON"]

##### RULES ######
rule all:
    input:
        # make datasets
        expand(os.path.join(PREP_DIR,'event_psi','Riaz2017-{condition}-{omic_type}.tsv.gz'), condition=CONDITIONS, omic_type=OMIC_TYPES),
        
        # survival analysis
        expand(os.path.join(RESULTS_DIR,'files',"survival_analysis",'splicing-{omic_type}-Riaz2017-{condition}-surv.tsv.gz'), condition=CONDITIONS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,'files',"survival_analysis",'splicing-{omic_type}-Riaz2017-{condition}-cat.tsv.gz'), condition=CONDITIONS, omic_type=OMIC_TYPES),
        
        # figures
        expand(os.path.join(RESULTS_DIR,"figures","immune_evasion-{omic_type}"), omic_type=OMIC_TYPES)        
        
        
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
        

rule survival_analysis_splicing:
    input:
        splicing = os.path.join(PREP_DIR,'event_psi','Riaz2017-{condition}-{omic_type}.tsv.gz'),
        metadata = os.path.join(PREP_DIR,"metadata","Riaz2017.tsv.gz")
    output:
        surv = os.path.join(RESULTS_DIR,'files',"survival_analysis",'splicing-{omic_type}-Riaz2017-{condition}-surv.tsv.gz'),
        cat = os.path.join(RESULTS_DIR,'files',"survival_analysis",'splicing-{omic_type}-Riaz2017-{condition}-cat.tsv.gz')
    params:
        sample_col = "sampleID",
        surv_event_col = "OS_event",
        surv_time_col = "OS_time"
    shell:
        """
        Rscript scripts/survival_analysis.R \
                    --data_matrix_file={input.splicing} \
                    --metadata_file={input.metadata} \
                    --sample_col={params.sample_col} \
                    --surv_event_col={params.surv_event_col} \
                    --surv_time_col={params.surv_time_col} \
                    --output_surv_file={output.surv} \
                    --output_cat_file={output.cat}
        """
        
        
rule figures_immune_evasion:
    input:
        annotation = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz'),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}"),
        splicing = os.path.join(PREP_DIR,"event_psi","Riaz2017-PRE-{omic_type}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","Riaz2017.tsv.gz"),
        enrichments_reactome = os.path.join(RESULTS_DIR,"figures","cancer_program","figdata","cancer_program","enrichments_reactome.tsv.gz"),
        immune_screen = os.path.join(SUPPORT_DIR,"supplementary_tables_literature","Dubrot2022-suptabs-41590_2022_1315_MOESM2_ESM.xlsx"), # Sup. Tab. 13
        human2mouse = os.path.join(RAW_DIR,"BIOMART","human2mouse.tsv"),
        survival_analysis = os.path.join(RESULTS_DIR,'files',"survival_analysis",'splicing-{omic_type}-Riaz2017-PRE-surv.tsv.gz'),
        protein_impact = os.path.join(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
    output:
        directory(os.path.join(RESULTS_DIR,"figures","immune_evasion-{omic_type}"))
    shell:
        """
        Rscript scripts/figures_immune_evasion.R \
                    --annotation_file={input.annotation} \
                    --regulons_path={input.regulons_path} \
                    --splicing_file={input.splicing} \
                    --metadata_file={input.metadata} \
                    --enrichments_reactome_file={input.enrichments_reactome} \
                    --immune_screen_file={input.immune_screen} \
                    --human2mouse_file={input.human2mouse} \
                    --survival_analysis_file={input.survival_analysis} \
                    --protein_impact_file={input.protein_impact} \
                    --figs_dir={output}
        """