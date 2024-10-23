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

# variables
CONFOUNDER_VARS = c(
    "age_at_initial_pathologic_diagnosis",
    "gender",
    "ajcc_pathologic_tumor_stage"
)

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","sf_activity_tcga")
# data_matrix_file = file.path(RESULTS_DIR,"files","protein_activity","LAML-PrimaryBloodDerivedCancerPeripheralBlood-EX.tsv.gz")
# data_matrix_file = file.path(PREP_DIR,"genexpr_tpm","LAML-PrimaryBloodDerivedCancerPeripheralBlood.tsv.gz")
# data_matrix_file = file.path(PREP_DIR,"genexpr_tpm","CHOL-PrimaryTumor.tsv.gz")
# metadata_file = file.path(RAW_DIR,'UCSCXena','TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz')
# surv_event_col = "OS"
# surv_time_col = "OS.time"
# sample_col = "sample"
# consider_confounders = "True"
# features_oi_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors-ensembl.txt")

##### FUNCTIONS #####
run_cutpoint = function(df, surv_time_col, surv_event_col, consider_confounders){
    result = tryCatch(
        {
            
            if (consider_confounders=="True"){
                variables = c("feature_value", CONFOUNDER_VARS)
                vars_formula = paste(variables, collapse=" + ")
            } else {
                variables = "feature_value"
                vars_formula = "feature_value"
            }
            
            result_cut = surv_cutpoint(df, time=surv_time_col, event=surv_event_col, variables=variables)
            result_cat = surv_categorize(result_cut) %>% data.frame()
            fit_survdiff = survival::survdiff(
                as.formula(sprintf("Surv(%s, %s) ~ %s", surv_time_col, surv_event_col, vars_formula)), 
                data=result_cat
            )
            pvalue_survdiff = stats::pchisq(fit_survdiff[["chisq"]], length(fit_survdiff[["n"]]) - 1, lower.tail = FALSE)
            fit_coxph = survival::coxph(
                as.formula(sprintf("Surv(%s, %s) ~ %s", surv_time_col, surv_event_col, vars_formula)), 
                data=df
            )
            result_surv = data.frame(
                feature = df[["feature"]] %>% unique(),
                feature_value_cutoff = result_cut[["feature_value"]][["estimate"]],
                survdiff_pvalue = pvalue_survdiff,
                survdiff_method = "Log-rank",
                coxph_coef = fit_coxph[["coefficients"]][["feature_value"]],
                coxph_pvalue = summary(fit_coxph)[["waldtest"]][["pvalue"]]
            )
            result_cat[["feature"]] = df[["feature"]] %>% unique()
            result_cat[["sampleID"]] = df[["sampleID"]]

            result = list(
                "surv" = result_surv,
                "cat" = result_cat
            )
            return(result)
        }, error=function(cond){
            result_surv = data.frame(
                    feature = df[["feature"]] %>% unique(),
                    feature_value_cutoff = NA,
                    survdiff_pvalue = NA,
                    survdiff_method = NA,
                    coxph_coef = NA,
                    coxph_pvalue = NA
                )
            
            result_cat = data.frame(
                    df[[surv_time_col]],
                    df[[surv_event_col]],
                    NA,
                    df[["feature"]] %>% unique(),
                    df[["sampleID"]]
                )
            colnames(result_cat) = c(surv_time_col, surv_event_col, "feature_value", "feature", "sampleID")
            
            
            result = list(
                "surv" = result_surv,
                "cat" = result_cat
            )
        }
    )
    return(result)
}


parseargs = function(){
    
    option_list = list( 
        make_option("--data_matrix_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--sample_col", type="character"),
        make_option("--surv_event_col", type="character"),
        make_option("--surv_time_col", type="character"),
        make_option("--features_oi_file", type="character", default=NULL),
        make_option("--consider_confounders", type="character", default="False"),
        make_option("--output_surv_file", type="character"),
        make_option("--output_cat_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    data_matrix_file = args[["data_matrix_file"]]
    metadata_file = args[["metadata_file"]]
    sample_col = args[["sample_col"]]
    surv_event_col = args[["surv_event_col"]]
    surv_time_col = args[["surv_time_col"]]
    features_oi_file = args[["features_oi_file"]]
    consider_confounders = args[["consider_confounders"]]
    output_surv_file = args[["output_surv_file"]]
    output_cat_file = args[["output_cat_file"]]
    
    # load
    data_matrix = read_tsv(data_matrix_file)
    metadata = read_tsv(metadata_file)
    
    # prep
    ## rename first column
    data_matrix = data_matrix %>%
        dplyr::rename("feature"=colnames(data_matrix)[1])
    
    ## filter
    if (!is.null(features_oi_file)){
        print("Filtering features...")
        features_oi = readLines(features_oi_file)
        
        data_matrix = data_matrix %>%
            filter(feature %in% features_oi)
    }
    
    # compute KM into two groups
    metadata_cols = c(sample_col, surv_event_col, surv_time_col)

    if (consider_confounders=="True") {
        metadata_cols = c(metadata_cols, CONFOUNDER_VARS)
    }
    
    X = data_matrix %>%
        pivot_longer(-feature, names_to="sampleID", values_to="feature_value") %>%
        left_join(
            metadata[,metadata_cols] %>% distinct(),
            by=c("sampleID"=sample_col)
        ) %>%
        drop_na()
    
    result = lapply(
            X %>% pull(feature) %>% unique(), function(feature_oi){
            df = X %>% filter(feature==feature_oi)
            res = run_cutpoint(df, surv_time_col, surv_event_col, consider_confounders)
        })
    
    result_surv = lapply(result, function(x){x[["surv"]]}) %>% bind_rows()
    result_surv[["survdiff_fdr"]] = p.adjust(result_surv[["survdiff_pvalue"]], method="fdr")
    result_surv[["coxph_fdr"]] = p.adjust(result_surv[["coxph_pvalue"]], method="fdr")
    result_cat = lapply(result, function(x){x[["cat"]]}) %>% bind_rows()
    
    # save
    write_tsv(result_surv, output_surv_file)
    write_tsv(result_cat, output_cat_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}