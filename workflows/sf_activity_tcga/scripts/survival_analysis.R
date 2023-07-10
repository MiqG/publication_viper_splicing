#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - Validate whether protein activity inferrence indicates that indisulam targets RBM39
# - Investigate the effect of MS023, inhibitor of Type I PRMT enzymes
#     - no cell growth effect in vitro, but strong suppression of tumor growth in vivo

require(optparse)
require(tidyverse)
require(survival)
require(survminer)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_tcga")
# protein_activity_file = file.path(RESULTS_DIR,"files","protein_activity","LAML-PrimaryBloodDerivedCancerPeripheralBlood-EX.tsv.gz")
# metadata_file = file.path(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz')
# surv_event_col = "OS"
# surv_time_col = "OS.time"
# sample_col = "sample"


##### FUNCTIONS #####
run_cutpoint = function(df, surv_time_col, surv_event_col){
    result_cut = surv_cutpoint(df, time=surv_time_col, event=surv_event_col, variables="activity")
    result_cat = surv_categorize(result_cut)
    fit_survdiff = survival::survdiff(
        as.formula(sprintf("Surv(%s, %s) ~ activity", surv_time_col, surv_event_col)), 
        data=result_cat
    )
    pvalue_survdiff = stats::pchisq(fit_survdiff[["chisq"]], length(fit_survdiff[["n"]]) - 1, lower.tail = FALSE)
    fit_coxph = survival::coxph(
        as.formula(sprintf("Surv(%s, %s) ~ activity", surv_time_col, surv_event_col)), 
        data=df
    )
    result_surv = data.frame(
        regulator = df[["regulator"]] %>% unique(),
        activity_cutoff = result_cut[["activity"]][["estimate"]],
        survdiff_pvalue = pvalue_survdiff,
        survdiff_method = "Log-rank",
        coxph_coef = fit_coxph[["coefficients"]],
        coxph_pvalue = summary(fit_coxph)[["waldtest"]][["pvalue"]]
    )
    result_cat[["regulator"]] = df[["regulator"]] %>% unique()
    result_cat[["sampleID"]] = df[["sampleID"]]
    
    result = list(
        "surv" = result_surv,
        "cat" = result_cat
    )
    return(result)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--protein_activity_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--sample_col", type="character"),
        make_option("--surv_event_col", type="character"),
        make_option("--surv_time_col", type="character"),
        make_option("--output_surv_file", type="character"),
        make_option("--output_cat_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    protein_activity_file = args[["protein_activity_file"]]
    metadata_file = args[["metadata_file"]]
    sample_col = args[["sample_col"]]
    surv_event_col = args[["surv_event_col"]]
    surv_time_col = args[["surv_time_col"]]
    output_surv_file = args[["output_surv_file"]]
    output_cat_file = args[["output_cat_file"]]
    
    # load
    protein_activity = read_tsv(protein_activity_file)
    metadata = read_tsv(metadata_file)
    
    # compute KM into two groups
    X = protein_activity %>%
        pivot_longer(-regulator, names_to="sampleID", values_to="activity") %>%
        left_join(
            metadata[,c(sample_col, surv_event_col, surv_time_col)] %>% distinct(),
            by=c("sampleID"=sample_col)
        ) %>%
        drop_na()
    
    result = lapply(
            X %>% pull(regulator) %>% unique(), function(regulator_oi){
            df = X %>% filter(regulator==regulator_oi)
            res = run_cutpoint(df, surv_time_col, surv_event_col)
        })
    
    result_surv = lapply(result, function(x){x[["surv"]]}) %>% do.call(rbind,.)
    result_surv[["survdiff_fdr"]] = p.adjust(result_surv[["survdiff_pvalue"]], method="fdr")
    result_surv[["coxph_fdr"]] = p.adjust(result_surv[["coxph_pvalue"]], method="fdr")
    result_cat = lapply(result, function(x){x[["cat"]]}) %>% do.call(rbind,.)
    
    # save
    write_tsv(result_surv, output_surv_file)
    write_tsv(result_cat, output_cat_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}