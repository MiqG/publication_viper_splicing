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
OMIC_TYPES = EVENT_TYPES

##### RULES #####
rule all:
    input:
        # calculate signature
        expand(os.path.join(RESULTS_DIR,"files","signatures","combinatorial_perturbations-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # compute viper SF activities
        expand(os.path.join(RESULTS_DIR,"files","protein_activity","combinatorial_perturbations-{omic_type}.tsv.gz"), omic_type=OMIC_TYPES),
        
        # figures
        os.path.join(RESULTS_DIR,"figures","validation_combinatorial_perturbations")
        
        
rule compute_signatures:
    input:
        metadata = os.path.join(PREP_DIR,'metadata','ENASFS.tsv.gz'),
        splicing = os.path.join(PREP_DIR,'event_psi','ENASFS-{omic_type}.tsv.gz')
    output:
        signatures = os.path.join(RESULTS_DIR,"files","signatures","combinatorial_perturbations-{omic_type}.tsv.gz")
    run:
        import pandas as pd
        
        # load
        metadata = pd.read_table(input.metadata)
        splicing = pd.read_table(input.splicing, index_col=0)
        
        # select combinatorial perturbations
        perts_oi = ["KNOCKDOWN_DOUBLE","KNOCKOUT_DOUBLE","KNOCKDOWN_TRIPLE","KNOCKDOWN_QUADRUPLE"]
        metadata = metadata.loc[metadata["PERT_TYPE"].isin(perts_oi)]
        
        # delta PSI as the difference between conditions and the mean of the conditions
        signatures = {}
        for sample_oi in metadata["sampleID"]:
            # get the controls of the sample
            ctls = metadata.loc[metadata["sampleID"]==sample_oi, "control_samples"].values[0]
            
            # controls will be np.nan
            if isinstance(ctls,str):
                ctls = ctls.split("||")
                psi_ctls = splicing[ctls].mean(axis=1)
                
                # compute delta PSI
                dpsi = splicing[sample_oi] - psi_ctls
                
                signatures[sample_oi] = dpsi
                
                del dpsi, psi_ctls, ctls
                
        
        signatures = pd.DataFrame(signatures)
        
        # save
        signatures.reset_index().to_csv(output.signatures, **SAVE_PARAMS)
        
        print("Done!")

        
rule compute_protein_activity:
    input:
        signature = os.path.join(RESULTS_DIR,"files","signatures","combinatorial_perturbations-{omic_type}.tsv.gz"),
        regulons_path = os.path.join(REGULONS_DIR,"files","experimentally_derived_regulons_pruned-{omic_type}")
    output:
        os.path.join(RESULTS_DIR,"files","protein_activity","combinatorial_perturbations-{omic_type}.tsv.gz")
    params:
        script_dir = SRC_DIR
    shell:
        """
        Rscript {params.script_dir}/compute_protein_activity.R \
                    --signature_file={input.signature} \
                    --regulons_path={input.regulons_path} \
                    --output_file={output}
       """
        
        
rule make_figures:
    input:
        protein_activity = os.path.join(RESULTS_DIR,"files","protein_activity","combinatorial_perturbations-EX.tsv.gz"),
        metadata = os.path.join(PREP_DIR,"metadata","ENASFS.tsv.gz"),
        splicing_factors = os.path.join(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,"figures","validation_combinatorial_perturbations"))
    shell:
        """
        Rscript scripts/figures_combinatorial_perturbations.R \
                    --protein_activity_file={input.protein_activity} \
                    --metadata_file={input.metadata} \
                    --splicing_factors_file={input.splicing_factors} \
                    --figs_dir={output}
        """
