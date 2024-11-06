import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]

##### RULES #####
rule all:
    input:
        # run gene set enrichment analysis
        expand(os.path.join(RESULTS_DIR,"files","regulons_eda_gsea","experimentally_derived_regulons_pruned-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        # jaccard distances
        expand(os.path.join(RESULTS_DIR,"files","regulons_eda_jaccard","experimentally_derived_regulons_pruned-{event_type}.tsv.gz"), event_type=EVENT_TYPES),
        
        # figures
        expand(os.path.join(RESULTS_DIR,'figures','eda_regulons-{event_type}'), event_type=EVENT_TYPES)
        
        
rule run_gsea:
    input:
        regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}"),
        annotation = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv'),
        msigdb_dir = os.path.join(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs'),
        protein_impact = os.path.join(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
    output:
        os.path.join(RESULTS_DIR,"files","regulons_eda_gsea","experimentally_derived_regulons_pruned-{event_type}.tsv.gz")
    shell:
        """
        Rscript scripts/run_gsea.R \
                    --regulons_path={input.regulons} \
                    --annotation_file={input.annotation} \
                    --msigdb_dir={input.msigdb_dir} \
                    --protein_impact_file={input.protein_impact} \
                    --output_file={output}
        """

        
rule compute_jaccard_distances:
    input:
        regulons_dir = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}")
    output:
        dist = os.path.join(RESULTS_DIR,"files","regulons_eda_jaccard","experimentally_derived_regulons_pruned-{event_type}.tsv.gz")
    run:
        import os
        import pandas as pd
        import numpy as np
        from sklearn.metrics import pairwise_distances
        
        regulons = pd.concat([pd.read_table(os.path.join(input.regulons_dir,f)) for f in os.listdir(input.regulons_dir) if f.endswith(".tsv.gz")])
        regulons = regulons[["regulator","target"]].drop_duplicates().copy()
        regulons["edge"] = True
        regulons_mat = pd.pivot_table(
            regulons, values="edge", index="target", columns="regulator", fill_value=False
        ).astype(bool)
        
        dist = pairwise_distances(regulons_mat.T.values, metric="jaccard")
        dist = pd.DataFrame(dist, index=regulons_mat.columns, columns=regulons_mat.columns)
        dist = dist.where(np.triu(np.ones(dist.shape), k=1).astype(bool))
        dist = dist.stack()
        dist.name = "distance"
        dist.index.names = ["regulator_a","regulator_b"]
        dist = dist.reset_index()
        dist["method"] = "jaccard"

        dist.to_csv(output.dist, **SAVE_PARAMS)
        
        print("Done!")
        

rule figures_eda_regulons:
    input:
        regulons = os.path.join(RESULTS_DIR,"files","experimentally_derived_regulons_pruned-{event_type}"),
        regulons_clip = os.path.join(RESULTS_DIR,"files","postar3_clip_regulons-EX"),
        protein_impact = os.path.join(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz'),
        annotation = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv'),
        splicing_factors = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,'figures','eda_regulons-{event_type}'))
    shell:
        """
        Rscript scripts/figures_eda_regulons.R \
                    --regulons_dir={input.regulons} \
                    --regulons_clip_dir={input.regulons_clip} \
                    --protein_impact_file={input.protein_impact} \
                    --annotation_file={input.annotation} \
                    --splicing_factors_file={input.splicing_factors} \
                    --figs_dir={output}
        """