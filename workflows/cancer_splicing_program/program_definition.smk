import os
import pandas as pd

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
SRC_DIR = os.path.join(ROOT,"src")
RESULTS_DIR = os.path.join(ROOT,"results","cancer_splicing_program")
REGULONS_DIR = os.path.join(ROOT,"results","regulon_inference")
PACT_CCLE_DIR = os.path.join(ROOT,"results","sf_activity_ccle")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}
PADJ_METHOD = 'fdr_bh'

metadata = pd.read_table(os.path.join(RAW_DIR,"TCGA","metadata","PANCAN.tsv.gz"))
metadata = metadata.groupby(
    ["sample_type_clean","cancer_type"]
).size().reset_index().rename(columns={0:"n"})
metadata = metadata.loc[metadata["n"]>=3]

CANCER_TYPES = metadata["cancer_type"].values
SAMPLE_TYPES = metadata["sample_type_clean"].values

CANCER_TYPES_PTSTN = [
    'BLCA',
    'BRCA',
    'COAD',
    'HNSC',
    'KICH',
    'KIRC',
    'KIRP',
    'LIHC',
    'LUAD',
    'LUSC',
    'PRAD',
    'STAD',
    'THCA',
    'UCEC'
]
CANCER_TYPES_METPT = ["BRCA","SKCM","THCA"]
CANCER_TYPES_PT = metadata.loc[
    metadata["sample_type_clean"].isin(["PrimaryTumor","PrimaryBloodDerivedCancerPeripheralBlood"]),"cancer_type"
].values
CANCER_TYPES_STN = metadata.loc[metadata["sample_type_clean"]=="SolidTissueNormal","cancer_type"].values

DIFF_CONDITIONS = {
    "PrimaryTumor_vs_SolidTissueNormal":{
        "a":"PrimaryTumor",
        "b":"SolidTissueNormal"
    },
    "Metastatic_vs_PrimaryTumor":{
        "a":"Metastatic",
        "b":"PrimaryTumor"
    }
}

DIFF_CANCER_TYPES = {
    "PrimaryTumor_vs_SolidTissueNormal": CANCER_TYPES_PTSTN,
    "Metastatic_vs_PrimaryTumor": CANCER_TYPES_METPT
}

CANCER_TYPES_BY_SAMPLETYPE = {
    s: metadata.loc[metadata["sample_type_clean"]==s,"cancer_type"].values for s in metadata["sample_type_clean"].unique()
}

CANCER_TYPES_SURV = metadata.loc[
    metadata["sample_type_clean"].isin(["PrimaryTumor","PrimaryBloodDerivedCancerPeripheralBlood"]),
    ["cancer_type","sample_type_clean"]
]
OMICS = ["protein_activity","genexpr_tpm"]
##### RULES #####
rule all:
    input:
        # calculate signatures
        expand(os.path.join(RESULTS_DIR,"files","signatures","{cancer}-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz"), cancer=CANCER_TYPES_PTSTN),
        expand(os.path.join(RESULTS_DIR,"files","signatures","{cancer}-Metastatic_vs_PrimaryTumor-EX.tsv.gz"), cancer=CANCER_TYPES_METPT),
        expand(os.path.join(RESULTS_DIR,"files","signatures","{cancer}-{sample}-EX.tsv.gz"), zip, cancer=CANCER_TYPES, sample=SAMPLE_TYPES),

        # compute viper SF activities
        ## PT vs STN
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz"), cancer=CANCER_TYPES_PTSTN, sample=["PrimaryTumor_vs_SolidTissueNormal"]),
        ## MET vs PT
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz"), cancer=CANCER_TYPES_METPT, sample=["Metastatic_vs_PrimaryTumor"]),
        ## within
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz"), zip, cancer=CANCER_TYPES, sample=SAMPLE_TYPES),

        # differential analyses
        ## SF activities
        expand(os.path.join(RESULTS_DIR,'files',"diff_protein_activity",'{cancer}-{comparison}.tsv.gz'), comparison=["PrimaryTumor_vs_SolidTissueNormal"], cancer=CANCER_TYPES_PTSTN),
        expand(os.path.join(RESULTS_DIR,'files',"diff_protein_activity",'{cancer}-{comparison}.tsv.gz'), comparison=["Metastatic_vs_PrimaryTumor"], cancer=CANCER_TYPES_METPT),
        ## gene expression
        expand(os.path.join(RESULTS_DIR,'files',"diff_genexpr_tpm",'{cancer}-{comparison}.tsv.gz'), comparison=["PrimaryTumor_vs_SolidTissueNormal"], cancer=CANCER_TYPES_PTSTN),
        expand(os.path.join(RESULTS_DIR,'files',"diff_genexpr_tpm",'{cancer}-{comparison}.tsv.gz'), comparison=["Metastatic_vs_PrimaryTumor"], cancer=CANCER_TYPES_METPT),
        ## merge
        expand(os.path.join(RESULTS_DIR,'files','PANCAN','{omic}-mannwhitneyu-{comparison}.tsv.gz'), comparison=DIFF_CANCER_TYPES.keys(), omic=OMICS),
        ## define cancer program
        os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz'),
        
        # survival
        ## SF activity
        expand(os.path.join(RESULTS_DIR,'files',"survival_analysis",'{omic}-{cancer}-{sample_type}-surv.tsv.gz'), zip, cancer=CANCER_TYPES_SURV["cancer_type"].values, sample_type=CANCER_TYPES_SURV["sample_type_clean"].values, omic=["protein_activity"]),
        expand(os.path.join(RESULTS_DIR,'files',"survival_analysis",'{omic}-{cancer}-{sample_type}-cat.tsv.gz'), zip, cancer=CANCER_TYPES_SURV["cancer_type"].values, sample_type=CANCER_TYPES_SURV["sample_type_clean"].values, omic=["protein_activity"]),
        ## gene expression
        expand(os.path.join(RESULTS_DIR,'files',"survival_analysis",'{omic}-{cancer}-{sample_type}-surv.tsv.gz'), zip, cancer=CANCER_TYPES_SURV["cancer_type"].values, sample_type=CANCER_TYPES_SURV["sample_type_clean"].values, omic=["genexpr_tpm"]),
        expand(os.path.join(RESULTS_DIR,'files',"survival_analysis",'{omic}-{cancer}-{sample_type}-cat.tsv.gz'), zip, cancer=CANCER_TYPES_SURV["cancer_type"].values, sample_type=CANCER_TYPES_SURV["sample_type_clean"].values, omic=["genexpr_tpm"]),
        ## merge
        expand(os.path.join(RESULTS_DIR,'files','PANCAN',"{omic}-survival_analysis-{surv_type}.tsv.gz"), surv_type=["surv","cat"], omic=OMICS),

        # cross regulation
        expand(os.path.join(RESULTS_DIR,'files',"sf_cross_regulation",'{omic}-{cancer}-{sample_type}.tsv.gz'), zip, cancer=CANCER_TYPES_SURV["cancer_type"].values, sample_type=CANCER_TYPES_SURV["sample_type_clean"].values, omic=OMICS),
        expand(os.path.join(RESULTS_DIR,'files','PANCAN',"{omic}-sf_cross_regulation.tsv.gz"), omic=OMICS),

        # correl genexpr vs activity
        expand(os.path.join(RESULTS_DIR,'files',"sf_activity_regulation",'genexpr_tpm_vs_activity-{cancer}-{sample_type}.tsv.gz'), zip, cancer=CANCER_TYPES_SURV["cancer_type"].values, sample_type=CANCER_TYPES_SURV["sample_type_clean"].values),
        os.path.join(RESULTS_DIR,'files','PANCAN',"genexpr_tpm_vs_activity.tsv.gz"),

        # figures
        os.path.join(RESULTS_DIR,"figures","cancer_program")
        

rule compute_signature_pt_vs_stn:
    input:
        splicing_pt = os.path.join(PREP_DIR,"event_psi","{cancer}-PrimaryTumor-EX.tsv.gz"),
        splicing_stn = os.path.join(PREP_DIR,"event_psi","{cancer}-SolidTissueNormal-EX.tsv.gz")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{cancer}-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz")
    run:
        import pandas as pd
        
        splicing_pt = pd.read_table(input.splicing_pt, index_col=0)
        splicing_stn = pd.read_table(input.splicing_stn, index_col=0)
        
        # subtract median from the other group
        signature = pd.concat([
            splicing_pt - splicing_stn.median(axis=1).values.reshape(-1,1), 
            splicing_stn - splicing_pt.median(axis=1).values.reshape(-1,1)
        ], axis=1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
        
        
rule compute_signature_met_vs_pt:
    input:
        splicing_pt = os.path.join(PREP_DIR,"event_psi","{cancer}-PrimaryTumor-EX.tsv.gz"),
        splicing_met = os.path.join(PREP_DIR,"event_psi","{cancer}-Metastatic-EX.tsv.gz")
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{cancer}-Metastatic_vs_PrimaryTumor-EX.tsv.gz")
    run:
        import pandas as pd
        
        splicing_pt = pd.read_table(input.splicing_pt, index_col=0)
        splicing_met = pd.read_table(input.splicing_met, index_col=0)
        
        # subtract median PT
        signature = pd.concat([
            splicing_pt - splicing_met.median(axis=1).values.reshape(-1,1), 
            splicing_met - splicing_pt.median(axis=1).values.reshape(-1,1)
        ], axis=1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")
        

ruleorder: compute_signature_pt_vs_stn > compute_signature_within
ruleorder: compute_signature_met_vs_pt > compute_signature_within
rule compute_signature_within:
    input:
        splicing_pt = os.path.join(PREP_DIR,"event_psi","{cancer}-{sample}-EX.tsv.gz"),
    output:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{cancer}-{sample}-EX.tsv.gz")
    run:
        import pandas as pd
        
        splicing_pt = pd.read_table(input.splicing_pt, index_col=0)
        
        # subtract median PT
        signature = splicing_pt
        signature = signature - splicing_pt.median(axis=1).values.reshape(-1,1)
        
        # save
        signature.reset_index().to_csv(output.signature, **SAVE_PARAMS)
        
        print("Done!")

#ruleorder: combine_protein_activity > compute_protein_activity      
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","{cancer}-{sample}-EX.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-EX")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample}-EX.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
        """

        
rule compute_differential_protein_activity:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{comparison}-EX.tsv.gz"),
        metadata = os.path.join(PREP_DIR,'metadata','{cancer}.tsv.gz')
    output:
        os.path.join(RESULTS_DIR,'files',"diff_protein_activity",'{cancer}-{comparison}.tsv.gz')
    params:
        script_dir = SRC_DIR,
        padj_method = PADJ_METHOD,
        condition_a = lambda wildcards: DIFF_CONDITIONS[wildcards.comparison]["a"],
        condition_b = lambda wildcards: DIFF_CONDITIONS[wildcards.comparison]["b"],
        comparison_col = "sample_type_clean",
        sample_col = "sampleID"
    shell:
        """
        python {params.script_dir}/MannWhitneyU.py \
                    --data_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --sample_col={params.sample_col} \
                    --comparison_col={params.comparison_col} \
                    --condition_a={params.condition_a} \
                    --condition_b={params.condition_b} \
                    --output_file={output} \
                    --padj_method={params.padj_method} 
        """ 
        
        
rule compute_differential_genexpr:
    input:
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","{cancer}.tsv.gz"),
        metadata = os.path.join(PREP_DIR,'metadata','{cancer}.tsv.gz')
    output:
        os.path.join(RESULTS_DIR,'files',"diff_genexpr_tpm",'{cancer}-{comparison}.tsv.gz')
    params:
        script_dir = SRC_DIR,
        padj_method = PADJ_METHOD,
        condition_a = lambda wildcards: DIFF_CONDITIONS[wildcards.comparison]["a"],
        condition_b = lambda wildcards: DIFF_CONDITIONS[wildcards.comparison]["b"],
        comparison_col = "sample_type_clean",
        sample_col = "sampleID"
    shell:
        """
        python {params.script_dir}/MannWhitneyU.py \
                    --data_file={input.genexpr} \
                    --metadata_file={input.metadata} \
                    --sample_col={params.sample_col} \
                    --comparison_col={params.comparison_col} \
                    --condition_a={params.condition_a} \
                    --condition_b={params.condition_b} \
                    --output_file={output} \
                    --padj_method={params.padj_method} 
        """    
        
        
rule combine_differential_results:
    input:
        diff_files = lambda wildcards: [os.path.join(RESULTS_DIR,'files',"diff_{omic}",'{cancer}-{comparison}.tsv.gz').format(cancer=cancer, comparison=wildcards.comparison, omic="{omic}") for cancer in DIFF_CANCER_TYPES[wildcards.comparison]],
    params:
        comparison = '{comparison}'
    output:
        os.path.join(RESULTS_DIR,'files','PANCAN','{omic}-mannwhitneyu-{comparison}.tsv.gz')
    run:
        import os
        import pandas as pd
        
        dfs = []
        for diff_file in input.diff_files:
            # combine
            df = pd.read_table(diff_file)
            
            # add cancer type
            cancer_type = os.path.basename(diff_file).replace(".tsv.gz","").split("-")[0]
            print(cancer_type)
            df['cancer_type'] = cancer_type
            
            dfs.append(df)
            
            del df
            
        dfs = pd.concat(dfs)  
        dfs.to_csv(output[0], sep='\t', index=False, compression='gzip')
        
        print("Done!")
        

rule define_cancer_program:
    input:
        diff_activity = os.path.join(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz'),
        gene_annotation = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
    output:
        cancer_program = os.path.join(RESULTS_DIR,'files','PANCAN','cancer_program.tsv.gz')
    params:
        thresh_fdr = 0.05,
        thresh_n_sum = 5
    run:
        import pandas as pd
        import numpy as np
        
        # load data
        diff_activity = pd.read_table(input.diff_activity)
        gene_annotation = pd.read_table(input.gene_annotation)
        
        # oncogenic and tumor suppressor splicing factors are those that preferentially
        # activate or inactivate recurrently across cancer types
        
        ## keep significant
        diff_activity = diff_activity.loc[diff_activity["padj"] < params.thresh_fdr].copy()
        
        ## add driver type
        diff_activity["driver_type"] = diff_activity["condition_a-median"].apply(
            lambda x: "Oncogenic" if x>0 else "Tumor suppressor"
        )
        
        ## count recurrence
        cancer_program = diff_activity.groupby(
            ["regulator","driver_type"]
        ).size().reset_index().rename(columns={0:"n"})
        cancer_program["n_sign"] = [
            n if driver_type=="Oncogenic" else -n 
            for n, driver_type in cancer_program[["n","driver_type"]].values
        ]
        n_sum = cancer_program.groupby("regulator")["n_sign"].sum().reset_index().rename(
            columns={"n_sign":"n_sum"}
        )
        cancer_program = pd.merge(cancer_program, n_sum, on="regulator", how="left")
        
        ## classify
        cancer_program = cancer_program.loc[
            np.abs(cancer_program["n_sum"]) > params.thresh_n_sum
        ].copy()
        cancer_program = cancer_program.loc[
            cancer_program.groupby("regulator")["n"].idxmax()
        ].copy()
        
        # add gene annotations
        gene_annotation = gene_annotation.rename(
            columns={"Approved symbol":"GENE", "Ensembl gene ID":"ENSEMBL"}
        )
        cancer_program = cancer_program.rename(
            columns={"regulator":"ENSEMBL"}
        )
        cancer_program = pd.merge(
            cancer_program, gene_annotation[["GENE","ENSEMBL"]], on="ENSEMBL", how="left"
        )
        
        # save
        cancer_program.to_csv(output.cancer_program, **SAVE_PARAMS)
        
        print("Done!")

        
rule survival_analysis_protein_activity:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample_type}-EX.tsv.gz"),
        metadata = os.path.join(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz')
    output:
        surv = os.path.join(RESULTS_DIR,'files',"survival_analysis",'protein_activity-{cancer}-{sample_type}-surv.tsv.gz'),
        cat = os.path.join(RESULTS_DIR,'files',"survival_analysis",'protein_activity-{cancer}-{sample_type}-cat.tsv.gz')
    params:
        sample_col = "sample",
        surv_event_col = "OS",
        surv_time_col = "OS.time"
    shell:
        """
        Rscript scripts/survival_analysis.R \
                    --data_matrix_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --sample_col={params.sample_col} \
                    --surv_event_col={params.surv_event_col} \
                    --surv_time_col={params.surv_time_col} \
                    --output_surv_file={output.surv} \
                    --output_cat_file={output.cat}
        """
        
        
rule survival_analysis_genexpr_tpm:
    input:
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","{cancer}-{sample_type}.tsv.gz"),
        metadata = os.path.join(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz'),
        sfs = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-ensembl.txt")
    output:
        surv = os.path.join(RESULTS_DIR,'files',"survival_analysis",'genexpr_tpm-{cancer}-{sample_type}-surv.tsv.gz'),
        cat = os.path.join(RESULTS_DIR,'files',"survival_analysis",'genexpr_tpm-{cancer}-{sample_type}-cat.tsv.gz')
    params:
        sample_col = "sample",
        surv_event_col = "OS",
        surv_time_col = "OS.time"
    shell:
        """
        Rscript scripts/survival_analysis.R \
                    --data_matrix_file={input.genexpr} \
                    --metadata_file={input.metadata} \
                    --sample_col={params.sample_col} \
                    --surv_event_col={params.surv_event_col} \
                    --surv_time_col={params.surv_time_col} \
                    --features_oi_file={input.sfs} \
                    --output_surv_file={output.surv} \
                    --output_cat_file={output.cat}
        """
    
    
rule merge_survival_analysis:
    input:
        surv_files = [os.path.join(RESULTS_DIR,'files',"survival_analysis",'{omic}-{cancer}-{sample_type}-{surv_type}.tsv.gz').format(cancer=cancer, sample_type=sample_type, surv_type="{surv_type}", omic="{omic}") for cancer, sample_type in CANCER_TYPES_SURV.values]
    output:
        os.path.join(RESULTS_DIR,'files','PANCAN',"{omic}-survival_analysis-{surv_type}.tsv.gz")
    run:
        import os
        import pandas as pd
        
        # surv
        dfs = []
        for f in input.surv_files:
            cancer_type = os.path.basename(f).split("-")[1]
            sample_type = os.path.basename(f).split("-")[2]
            df = pd.read_table(f)
            df["cancer_type"] = cancer_type
            df["sample_type"] = sample_type
            dfs.append(df)
        
        dfs = pd.concat(dfs)
        dfs.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        

rule sf_cross_regulation_protein_activity:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample_type}-EX.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,'files',"sf_cross_regulation",'protein_activity-{cancer}-{sample_type}.tsv.gz')
    run:
        import pandas as pd
        import numpy as np
        
        # load
        protein_activity = pd.read_table(input.protein_activity, index_col=0)
        
        # compute correlations
        correlations = protein_activity.T.corr(method="pearson")
        correlations = correlations.where(np.triu(np.ones(correlations.shape), k=1).astype(bool))
        correlations = correlations.stack()
        correlations.name = "correlation"
        correlations.index.names = ["regulator_a","regulator_b"]
        correlations = correlations.reset_index()
        correlations["method"] = "pearson"
        
        # save
        correlations.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        
        
rule sf_cross_regulation_genexpr_tpm:
    input:
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","{cancer}.tsv.gz"),
        sfs = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-ensembl.txt")
    output:
        os.path.join(RESULTS_DIR,'files',"sf_cross_regulation",'genexpr_tpm-{cancer}-{sample_type}.tsv.gz')
    run:
        import pandas as pd
        import numpy as np
        
        # load
        genexpr = pd.read_table(input.genexpr, index_col=0)
        sfs = pd.read_table(input.sfs, header=None)[0].tolist()
        
        # prep
        genexpr = genexpr.loc[genexpr.index.isin(sfs)].copy()
        
        # compute correlations
        correlations = genexpr.T.corr(method="pearson")
        correlations = correlations.where(np.triu(np.ones(correlations.shape), k=1).astype(bool))
        correlations = correlations.stack()
        correlations.name = "correlation"
        correlations.index.names = ["regulator_a","regulator_b"]
        correlations = correlations.reset_index()
        correlations["method"] = "pearson"
        
        # save
        correlations.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
    
    
rule merge_sf_cross_regulation:
    input:
        corr_files = [os.path.join(RESULTS_DIR,'files',"sf_cross_regulation",'{omic}-{cancer}-{sample_type}.tsv.gz').format(cancer=cancer, sample_type=sample_type, omic="{omic}") for cancer, sample_type in CANCER_TYPES_SURV.values]
    output:
        os.path.join(RESULTS_DIR,'files','PANCAN',"{omic}-sf_cross_regulation.tsv.gz")
    run:
        import os
        import pandas as pd
        
        # surv
        dfs = []
        for f in input.corr_files:
            cancer_type = os.path.basename(f).split("-")[1]
            sample_type = os.path.basename(f).split("-")[2].replace(".tsv.gz","")
            df = pd.read_table(f)
            df["cancer_type"] = cancer_type
            df["sample_type"] = sample_type
            dfs.append(df)
        
        dfs = pd.concat(dfs)
        dfs.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
    

rule sf_activity_regulation_genexpr:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","{cancer}-{sample_type}-EX.tsv.gz"),
        genexpr = os.path.join(PREP_DIR,"genexpr_tpm","{cancer}.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,'files',"sf_activity_regulation",'genexpr_tpm_vs_activity-{cancer}-{sample_type}.tsv.gz')
    run:
        import pandas as pd
        from scipy import stats
        
        # load
        protein_activity = pd.read_table(input.protein_activity, index_col=0)
        genexpr = pd.read_table(input.genexpr, index_col=0)
        
        # prep
        common_sfs = set(protein_activity.index).intersection(genexpr.index)
        common_samples = set(protein_activity.columns).intersection(genexpr.columns)
        protein_activity = protein_activity.loc[common_sfs,common_samples].copy()
        genexpr = genexpr.loc[common_sfs,common_samples].copy()
        
        # compute correlations
        correlations = genexpr.apply(
            lambda row: protein_activity.corrwith(row, axis=1, method="spearman"), 
            axis=1
        )
        correlations = correlations.melt(ignore_index=False).reset_index()
        correlations.columns = ["sf_genexpr","sf_activity","correlation"]
        correlations["method"] = "spearman"
        
        # save
        correlations.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")
        
        
rule merge_sf_activity_regulation_genexpr:
    input:
        corr_files = [os.path.join(RESULTS_DIR,'files',"sf_activity_regulation",'genexpr_tpm_vs_activity-{cancer}-{sample_type}.tsv.gz').format(cancer=cancer, sample_type=sample_type) for cancer, sample_type in CANCER_TYPES_SURV.values]
    output:
        os.path.join(RESULTS_DIR,'files','PANCAN',"genexpr_tpm_vs_activity.tsv.gz")
    run:
        import os
        import pandas as pd
        
        # surv
        dfs = []
        for f in input.corr_files:
            cancer_type = os.path.basename(f).split("-")[1]
            sample_type = os.path.basename(f).split("-")[2].replace(".tsv.gz","")
            df = pd.read_table(f)
            df["cancer_type"] = cancer_type
            df["sample_type"] = sample_type
            dfs.append(df)
        
        dfs = pd.concat(dfs)
        dfs.to_csv(output[0], **SAVE_PARAMS)
        
        print("Done!")

        
rule figures_cancer_program:
    input:
        diff_activity = os.path.join(RESULTS_DIR,'files','PANCAN','protein_activity-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz'),
        diff_genexpr = os.path.join(RESULTS_DIR,'files','PANCAN','genexpr_tpm-mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal.tsv.gz'),
        survival_activity = os.path.join(RESULTS_DIR,'files','PANCAN',"protein_activity-survival_analysis-surv.tsv.gz"),
        survival_genexpr = os.path.join(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-survival_analysis-surv.tsv.gz"),
        sf_crossreg_activity = os.path.join(RESULTS_DIR,'files','PANCAN',"protein_activity-sf_cross_regulation.tsv.gz"),
        sf_crossreg_genexpr = os.path.join(RESULTS_DIR,'files','PANCAN',"genexpr_tpm-sf_cross_regulation.tsv.gz"),
        assocs_gene_dependency = os.path.join(PACT_CCLE_DIR,"files","protein_activity_vs_demeter2","CCLE.tsv.gz"),
        sf_activity_vs_genexpr = os.path.join(RESULTS_DIR,'files','PANCAN',"genexpr_tpm_vs_activity.tsv.gz"),
        ontology_chea = os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","PANCAN.tsv.gz"),
        annotation = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv'),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-EX"),
        regulons_jaccard = os.path.join(REGULONS_DIR,"files","regulons_eda_jaccard","experimentally_derived_regulons_pruned-EX.tsv.gz"),
        gene_annotation = os.path.join(RAW_DIR,"HGNC","gene_annotations.tsv.gz"),
        msigdb_dir = os.path.join(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs'),
        immune_screen = os.path.join(SUPPORT_DIR,"supplementary_tables_literature","Dubrot2022-suptabs-41590_2022_1315_MOESM2_ESM.xlsx"),
        human2mouse = os.path.join(RAW_DIR,"BIOMART","human2mouse.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","cancer_program"))
    shell:
        """
        Rscript scripts/figures_cancer_program.R \
                    --diff_activity_file={input.diff_activity} \
                    --diff_genexpr_file={input.diff_genexpr} \
                    --survival_activity_file={input.survival_activity} \
                    --survival_genexpr_file={input.survival_genexpr} \
                    --sf_crossreg_activity_file={input.sf_crossreg_activity} \
                    --sf_crossreg_genexpr_file={input.sf_crossreg_genexpr} \
                    --assocs_gene_dependency_file={input.assocs_gene_dependency} \
                    --sf_activity_vs_genexpr_file={input.sf_activity_vs_genexpr} \
                    --metadata_file={input.metadata} \
                    --regulons_path={input.regulons_path} \
                    --annotation_file={input.annotation} \
                    --regulons_jaccard_file={input.regulons_jaccard} \
                    --gene_annotation_file={input.gene_annotation} \
                    --ontology_chea_file={input.ontology_chea} \
                    --msigdb_dir={input.msigdb_dir} \
                    --immune_screen_file={input.immune_screen} \
                    --human2mouse_file={input.human2mouse} \
                    --figs_dir={output}
        """