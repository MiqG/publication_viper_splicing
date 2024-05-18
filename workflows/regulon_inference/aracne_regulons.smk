import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
PREP_DIR = os.path.join(ROOT,"data","prep")
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,"results","regulon_inference")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

EVENT_TYPES = ["EX"]
OMIC_TYPES = ["genexpr"] + EVENT_TYPES

GENEXPR_FILES = {
    "CardosoMoreira2020": os.path.join(PREP_DIR,'genexpr_tpm','CardosoMoreira2020.tsv.gz'),
}

SPLICING_FILES = {
    "CardosoMoreira2020": os.path.join(PREP_DIR,'event_psi_imputed','CardosoMoreira2020-{omic_type}.tsv.gz'),
}

DATASETS = list(SPLICING_FILES.keys())

# params
PARAMS = {
    "ARACNE_MI_THRESH_SEED": 1,
    "ARACNE_MI_THRESH_PVALUE": 0.00000001,
    "ARACNE_N_BOOTSTRAPS": 100,
    ## pruning
    "ARACNE_BOOTSTRAP_MAX_TARGETS": 1000,
    ## consolidation
    "ARACNE_CONSOLIDATE_PVALUE": 0.05,
}

N_BOOTSTRAPS = PARAMS["ARACNE_N_BOOTSTRAPS"]
BOOTSTRAPS = list(range(N_BOOTSTRAPS))
##### RULES #####
rule all:
    input:
        # make regulons
        expand(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_threshold"), dataset=DATASETS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_bootstrap_{boot_i}"), dataset=DATASETS, omic_type=OMIC_TYPES, boot_i=BOOTSTRAPS),
        expand(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_pruning"), dataset=DATASETS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned","network.txt"), dataset=DATASETS, omic_type=OMIC_TYPES),
        expand(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned","regulons.tsv.gz"), dataset=DATASETS, omic_type=OMIC_TYPES),
        # make regulon set
        expand(os.path.join(RESULTS_DIR,"files","aracne_regulons_development-{omic_type}"), omic_type=OMIC_TYPES)
        
               
# ----- ARACNe network reverse engineering -----
rule decompress_inputs:
    input:
        regulators = lambda wildcards: GENEXPR_FILES[wildcards.dataset],
        targets = lambda wildcards: SPLICING_FILES[wildcards.dataset] if wildcards.omic_type!="genexpr" else GENEXPR_FILES[wildcards.dataset],
    output:
        regulators = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","regulators.tsv"),
        targets = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","targets.tsv")
    shell:
        """
        set -eo pipefail
        
        zcat {input.regulators} > {output.regulators}
        zcat {input.targets} > {output.targets}
        
        echo "Done!"
        """

rule regulon_inference_aracne_java_threshold:
    input:
        regulators = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","regulators.tsv"),
        targets = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","targets.tsv"),
        regulators_oi = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-ensembl.txt")
    output:
        touch(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_threshold"))
    params:
        src = "scripts/aracne/ARACNe-AP",
        output_dir = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}"),
        random_seed = PARAMS["ARACNE_MI_THRESH_SEED"],
        mi_pvalue_thresh = str(PARAMS["ARACNE_MI_THRESH_PVALUE"]).replace("e","E"),
    threads: 12
    resources:
        runtime = 3600*12, # 12h 
        memory = 5, # 5GB
    shell:
        """
        set -eo pipefail
        
        java -Xmx{resources.memory}G -jar {params.src}/dist/aracne.jar \
                --expfile_upstream {input.regulators} \
                --tfs {input.regulators_oi} \
                --expfile_downstream {input.targets} \
                --output {params.output_dir} \
                --pvalue {params.mi_pvalue_thresh} \
                --seed {params.random_seed} \
                --threads {threads} \
                --calculateThreshold
                
        echo "Done!"
        """
        
        
rule regulon_inference_aracne_java_bootstrap:
    input:
        os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_threshold"),
        regulators = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","regulators.tsv"),
        targets = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","targets.tsv"),
        regulators_oi = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors-ensembl.txt")
    output:
        touch(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_bootstrap_{boot_i}"))
    params:
        src = "scripts/aracne/ARACNe-AP",
        output_dir = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}"),
        random_seed = "{boot_i}",
        mi_pvalue_thresh = str(PARAMS["ARACNE_MI_THRESH_PVALUE"]).replace("e","E")
    threads: 6
    resources:
        runtime = 3600*12, # 12h 
        memory = 5, # 5GB
    shell:
        """
        set -eo pipefail

        java -Xmx{resources.memory}G -jar {params.src}/dist/aracne.jar \
                --expfile_upstream {input.regulators} \
                --tfs {input.regulators_oi} \
                --expfile_downstream {input.targets} \
                --output {params.output_dir} \
                --pvalue {params.mi_pvalue_thresh} \
                --seed {params.random_seed} \
                --threads {threads} \
                --nodpi
                
        echo "Done!"
        """

        
rule regulon_inference_aracne_prune_bootstraps:
    input:
        os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_threshold"),
        [os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_bootstrap_{boot_i}").format(dataset="{dataset}", omic_type="{omic_type}", boot_i=boot_i) for boot_i in BOOTSTRAPS]
    output:
        touch(os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_pruning"))
    params:
        output_dir = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}"),
        max_targets = PARAMS["ARACNE_BOOTSTRAP_MAX_TARGETS"]
    threads: 1
    resources:
        runtime = 3600*1, # 1 h
        memory = 10 # 10 GB
    shell:
        """
        set -eo pipefail

        perl scripts/aracne/filter_arachne_bootstraps.pl \
                    --regulon-size {params.max_targets} \
                    {params.output_dir} \
                    {params.output_dir}/pruned

        echo "Done!"
        """


rule regulon_inference_aracne_java_consolidation:
    input:
        os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_threshold"),
        [os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_bootstrap_{boot_i}").format(dataset="{dataset}", omic_type="{omic_type}", boot_i=boot_i) for boot_i in BOOTSTRAPS],
        os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}",".done","done_pruning")
    output:
        os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned","network.txt")
    params:
        src = "scripts/aracne/ARACNe-AP",
        output_dir = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned"),
        consolidate_pvalue = PARAMS["ARACNE_CONSOLIDATE_PVALUE"]
    threads: 12
    resources:
        runtime = 3600*12, # 12h 
        memory = 20, # 15GB
    shell:
        """
        set -eo pipefail
        
        java -Xmx{resources.memory}G -jar {params.src}/dist/aracne.jar \
                --output {params.output_dir} \
                --threads {threads} \
                --consolidatepvalue {params.consolidate_pvalue} \
                --consolidate
                
        echo "Done!"
        """
        

rule prepare_regulons:
    input:
        regulators = lambda wildcards: GENEXPR_FILES[wildcards.dataset],
        targets = lambda wildcards: SPLICING_FILES[wildcards.dataset] if wildcards.omic_type!="genexpr" else GENEXPR_FILES[wildcards.dataset],
        aracne_network = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned","network.txt")
    output:
        regulons = os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned","regulons.tsv.gz")
    threads: 1
    resources:
        runtime = 3600*1, # 1 h
        memory = 10 # GB
    shell:
        """
        set -eo pipefail
        
        Rscript scripts/aracne/estimate_mor.R \
                    --genexpr_file={input.regulators} \
                    --splicing_file={input.targets} \
                    --aracne_network_file={input.aracne_network} \
                    --output_as_edgelist=TRUE \
                    --output_file={output.regulons}
                    
        echo "Done!"
        """
        
        
rule make_regulon_sets:
    input:
        regulons = [os.path.join(RESULTS_DIR,"files","aracne_regulons","{dataset}-{omic_type}","pruned","regulons.tsv.gz").format(dataset=d, omic_type="{omic_type}") for d in DATASETS]
    output:
        regulons_dir = directory(os.path.join(RESULTS_DIR,"files","aracne_regulons_development-{omic_type}"))
    run:
        import os
        import shutil
        
        os.makedirs(output.regulons_dir, exist_ok=True)
        for f in input.regulons:
            dataset = os.path.basename(os.path.dirname(os.path.dirname(f))).split("-")[0]
            filename = os.path.join(output.regulons_dir,"%s.tsv.gz") % dataset
            shutil.copy(f, filename)
            print("Copied", filename)
            
        print("Done!")
        