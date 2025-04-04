import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = EVENT_TYPES # ["genexpr"] + 

GENEXPR_FILES = {
    "CardosoMoreira2020": os.path.join(PREP_DIR,'genexpr_tpm','CardosoMoreira2020.tsv.gz'),
    "PANCAN_STN": os.path.join(PREP_DIR,'genexpr_tpm','PANCAN-SolidTissueNormal.tsv.gz'),
    "PANCAN_PT": os.path.join(PREP_DIR,'genexpr_tpm','PANCAN-PrimaryTumor.tsv.gz')
}

SPLICING_FILES = {
    "CardosoMoreira2020": os.path.join(PREP_DIR,'event_psi_imputed','CardosoMoreira2020-{omic_type}.tsv.gz'),
    "PANCAN_STN": os.path.join(PREP_DIR,'event_psi_imputed','PANCAN-SolidTissueNormal-{omic_type}.tsv.gz'),
    "PANCAN_PT": os.path.join(PREP_DIR,'event_psi_imputed','PANCAN-PrimaryTumor-{omic_type}.tsv.gz')
}

DATASETS = list(SPLICING_FILES.keys())


##### RULES #####
rule all:
    input:
        # make regulons
        expand(os.path.join(RESULTS_DIR,"files","mlr_regulons","{dataset}-{omic_type}","pruned","regulons.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        
        # make regulon set
        expand(os.path.join(RESULTS_DIR,"files","mlr_regulons_{dataset}-{omic_type}"), omic_type=OMIC_TYPES, dataset=DATASETS),
        
        # combine regulon sets
        expand(os.path.join(RESULTS_DIR,"files","mlr_regulons_combined-{omic_type}"), omic_type=OMIC_TYPES)
        
        
rule regulon_inference_mlr:
    input:
        regulators = lambda wildcards: GENEXPR_FILES[wildcards.dataset],
        targets = lambda wildcards: SPLICING_FILES[wildcards.dataset] if wildcards.omic_type!="genexpr" else GENEXPR_FILES[wildcards.dataset],
        regulators_oi = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-ensembl.txt"),
    output:
        os.path.join(RESULTS_DIR,"files","mlr_regulons","{dataset}-{omic_type}","pruned","regulons.tsv.gz")
    params:
        thresh_pvalue = 0.01
    threads: 24
    shell:
        """
        nice python scripts/run_mlr.py \
                    --regulators_file={input.regulators} \
                    --regulators_oi_file={input.regulators_oi} \
                    --targets_file={input.targets} \
                    --thresh_pvalue={params.thresh_pvalue} \
                    --n_jobs={threads} \
                    --output_file={output}
        """

        
rule make_regulon_sets:
    input:
        regulons = os.path.join(RESULTS_DIR,"files","mlr_regulons","{dataset}-{omic_type}","pruned","regulons.tsv.gz")
    output:
        regulons_dir = directory(os.path.join(RESULTS_DIR,"files","mlr_regulons_{dataset}-{omic_type}"))
    run:
        import os
        import shutil
        
        os.makedirs(output.regulons_dir, exist_ok=True)

        f = input.regulons
        dataset = os.path.basename(os.path.dirname(os.path.dirname(f))).split("-")[0]
        filename = os.path.join(output.regulons_dir,"%s.tsv.gz") % dataset
        shutil.copy(f, filename)
        print("Copied", filename)
            
        print("Done!")
        
        
rule combine_regulon_sets:
    input:
        regulons_dir = [os.path.join(RESULTS_DIR,"files","mlr_regulons_{dataset}-{omic_type}").format(dataset=d, omic_type="{omic_type}") for d in DATASETS]
    output:
        regulons_dir = directory(os.path.join(RESULTS_DIR,"files","mlr_regulons_combined-{omic_type}"))
    run:
        import os
        import shutil
        
        os.makedirs(output.regulons_dir, exist_ok=True)

        for regulons_dir in input.regulons_dir:
            # get all regulon files
            regulon_files = [os.path.join(regulons_dir,f) for f in os.listdir(regulons_dir) if f.endswith("tsv.gz")]
            
            # copy each regulon file
            for f in regulon_files:
                filename = os.path.join(output.regulons_dir, os.path.basename(f))
                shutil.copy(f, filename)
                print("Copied", filename)
        
        print("Done!")    