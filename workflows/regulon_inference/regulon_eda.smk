import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
BIN_DIR = os.path.join(ROOT,"bin")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]

##### RULES #####
rule all:
    input:
        # run gene set enrichment analysis
        expand(os.path.join(RESULTS_DIR,"regulons_eda_gsea","experimentally_derived_regulons_pruned-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        
rule run_gsea:
    input:
        regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}"),
        annotation = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv'),
        msigdb_dir = os.path.join(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs'),
        protein_impact = os.path.join(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
    output:
        os.path.join(RESULTS_DIR,"regulons_eda_gsea","experimentally_derived_regulons_pruned-{event_type}.tsv.gz")
    shell:
        """
        Rscript scripts/run_gsea.R \
                    --regulons_path={input.regulons} \
                    --annotation_file={input.annotation} \
                    --msigdb_dir={input.msigdb_dir} \
                    --protein_impact_file={input.protein_impact} \
                    --output_file={output}
        """