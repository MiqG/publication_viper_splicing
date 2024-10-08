#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables
REGULON_SETS = c(
    'aracne_regulons_CardosoMoreira2020',
    'mlr_regulons_CardosoMoreira2020',
    'aracne_regulons_PANCAN_STN',
    'mlr_regulons_PANCAN_STN',
    'experimentally_derived_regulons_pruned',
    'aracne_and_experimental_regulons',
    'mlr_and_experimental_regulons',
    'aracne_and_mlr_regulons',
    'top100_experimentally_derived_regulons_pruned',
    'top90_experimentally_derived_regulons_pruned',
    'top80_experimentally_derived_regulons_pruned',
    'top70_experimentally_derived_regulons_pruned',
    'top60_experimentally_derived_regulons_pruned',
    'top50_experimentally_derived_regulons_pruned',
    'top40_experimentally_derived_regulons_pruned'
)

SETS_MAIN = c(
    'aracne_regulons_CardosoMoreira2020',
    'mlr_regulons_CardosoMoreira2020',
    'aracne_regulons_PANCAN_STN',
    'mlr_regulons_PANCAN_STN',
    'splicinglore_regulons',
    'experimentally_derived_regulons_pruned'
)

SETS_ROBUSTNESS = c(
    'top40_experimentally_derived_regulons_pruned',
    'top50_experimentally_derived_regulons_pruned',
    'top60_experimentally_derived_regulons_pruned',
    'top70_experimentally_derived_regulons_pruned',
    'top80_experimentally_derived_regulons_pruned',
    'top90_experimentally_derived_regulons_pruned',
    'top100_experimentally_derived_regulons_pruned',
    'experimentally_derived_regulons_pruned'
)

SETS_THRESHOLDS = c(
    'dPSIthresh5_experimentally_derived_regulons_pruned',
    'dPSIthresh10_experimentally_derived_regulons_pruned',
    'dPSIthresh15_experimentally_derived_regulons_pruned',
    'dPSIthresh20_experimentally_derived_regulons_pruned',
    'dPSIthresh25_experimentally_derived_regulons_pruned',
    'dPSIthresh30_experimentally_derived_regulons_pruned',
    'dPSIthresh35_experimentally_derived_regulons_pruned',
    'dPSIthresh40_experimentally_derived_regulons_pruned',
    'dPSIthresh45_experimentally_derived_regulons_pruned'
)

SETS_LIKELIHOOD = c(
    'aracne_and_experimental_regulons',
    'mlr_and_experimental_regulons',
    'experimentally_derived_regulons_pruned'
)

SETS_MOR = c(
    'aracne_and_mlr_regulons'
)

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "#616161"
PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)

METHODS_ACTIVITY = c(
    "gsea",
    "correlation_spearman",
    "correlation_pearson",
    "viper"
)

SF_CLASS = c("Core", "RBP", "Other")

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,"support")
# RESULTS_DIR = file.path(ROOT,"results","regulon_inference")
# evaluation_ex_file = file.path(RESULTS_DIR,"files","regulon_evaluation_scores","merged-EX.tsv.gz")
# regulators_per_target_robustness_file = file.path(RESULTS_DIR,"files","regulon_properties","regulators_per_target-EX.tsv.gz")
# targets_per_regulator_robustness_file = file.path(RESULTS_DIR,"files","regulon_properties","targets_per_regulator-EX.tsv.gz")
# regulators_per_target_thresholds_file = file.path(RESULTS_DIR,"files","regulon_properties","dPSIthresh-regulators_per_target-EX.tsv.gz")
# targets_per_regulator_thresholds_file = file.path(RESULTS_DIR,"files","regulon_properties","dPSIthresh-targets_per_regulator-EX.tsv.gz")
# splicing_factors_file = file.path(SUPPORT_DIR,"splicing_factors","splicing_factors.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","regulon_evaluation")

##### FUNCTIONS #####
plot_evaluation = function(evaluation){
    plts = list()
    
    X = evaluation
    
    # general evaluation 
    ## across perturbations
    plts[["evaluation-general-perturbations-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        filter(curves_type=="combined" & eval_direction=="perturbations" & n_tails=="two") %>%
        ggplot(aes(x=method_activity, y=auc_roc, group=interaction(method_activity, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y = 0.3, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="ROC AUC | Perturbations", fill="Inference Type")
    
    
    ## across regulators
    plts[["evaluation-general-regulators-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        filter(curves_type=="combined" & eval_direction=="regulators" & n_tails=="two") %>%
        ggplot(aes(x=method_activity, y=auc_roc, group=interaction(method_activity, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y = 0.3, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="ROC AUC | Regulators", fill="Inference Type")
    
    # evaluation by single profile
    ## perturbations
    plts[["evaluation-single-perturbations-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        filter(curves_type=="by_group" & eval_direction=="perturbations" & n_tails=="two") %>%
        ggplot(aes(x=method_activity, y=auc_roc, group=interaction(method_activity, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y=-0.1, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="AUC ROC | Perturbation", fill="Inference Type")    
    
    ## regulators
    plts[["evaluation-single-regulators-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        filter(curves_type=="by_group" & eval_direction=="regulators" & n_tails=="two") %>%
        ggplot(aes(x=method_activity, y=auc_roc, group=interaction(method_activity, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y=-0.1, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="AUC ROC | Regulator", fill="Inference Type")        
    
    # batch effects
    ## between benchmark datasets
    plts[["evaluation-batch_effect_datasets-perturbations-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(curves_type=="by_group" & eval_direction=="perturbations" & n_tails=="two") %>%
        filter(regulon_set=="experimentally_derived_regulons_pruned") %>%
        ggplot(aes(x=signature_id, y=auc_roc, group=interaction(signature_id, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y=-0.1, label=label), 
            . %>% 
            count(signature_id, regulon_set, method_activity, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set+method_activity) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Benchmark Signature", y="AUC ROC | Perturbation", fill="Inference Type")
    
    plts[["evaluation-batch_effect_datasets-regulators-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(curves_type=="by_group" & eval_direction=="regulators" & n_tails=="two") %>%
        filter(regulon_set=="experimentally_derived_regulons_pruned") %>%
        ggplot(aes(x=signature_id, y=auc_roc, group=interaction(signature_id, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y=-0.1, label=label), 
            . %>% 
            count(signature_id, regulon_set, method_activity, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set+method_activity) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="AUC ROC | Regulator", fill="Inference Type")
    
    ## performance considering the number of networks per regulator
    plts[["evaluation-batch_effect_n_networks-regulators-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(curves_type=="by_group" & eval_direction=="regulators" & n_tails=="two") %>%
        filter(regulon_set=="experimentally_derived_regulons_pruned" & method_activity=="viper") %>%
        ggplot(aes(x=as.factor(n_networks_per_regulator), y=auc_roc, group=interaction(n_networks_per_regulator, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(color=eval_type), size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        color_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y=-0.1, label=label), 
            . %>% 
            count(n_networks_per_regulator, regulon_set, method_activity, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 0) +
        facet_wrap(~regulon_set+method_activity) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="N. Independent Signatures per Regulator", y="AUC ROC | Regulator", fill="Inference Type", color="Inference Type")
    
    ## performance of different regulators according to their general class
    plts[["evaluation-sf_class-regulators-auc_roc-violin"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY),
            sf_class = factor(sf_class, levels=SF_CLASS)
        ) %>%
        filter(curves_type=="by_group" & eval_direction=="regulators" & n_tails=="two") %>%
        filter(method_activity=="viper") %>%
        ggplot(aes(x=sf_class, y=auc_roc, group=interaction(sf_class, eval_type))) +
        geom_violin(aes(fill=eval_type), color=NA, trim=TRUE) +
        geom_point(aes(fill=eval_type), color=PAL_DARK, size=0.1, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        stat_summary(fun=median, geom="crossbar", linewidth=0.1, width=0.5, color="black", position=position_dodge(0.9)) + 
        fill_palette(PAL_EVAL_TYPE) + 
        geom_text(
            aes(y=-0.1, label=label), 
            . %>% 
            count(sf_class, regulon_set, method_activity, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~regulon_set+method_activity) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="AUC ROC | Regulator", fill="Inference Type") 
    
    ## (TODO) explore performance on diverse cell types through ENASFS
    X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(curves_type=="by_group" & eval_direction=="regulators" & n_tails=="two") %>%
        filter(regulon_set=="experimentally_derived_regulons_pruned" & signature_id=="ENASFS")
    
    ## (TODO) combined performance empirical networks with different thresholds
    
    ## (TODO) combined performance empirical netowkrs with different ablations
    
    ####### OLD ######
    
    X = evaluation %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, n_tails, regulon_set_id, pert_type_lab, regulator, n_targets_median, n_total_regulators, n_total_targets) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        filter(omic_type==omic_type_oi)
    
    # main networks
    plts[["evaluation-ranking_perc_vs_regulon_set_vs_pert_type-main-box"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="two") %>%
        ggplot(aes(x=pert_type_lab, y=ranking_perc, 
                   group=interaction(pert_type_lab, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~eval_direction+regulon_set_id, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(omic_type, pert_type_lab, regulon_set_id, eval_direction, eval_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Validation Perturbation", y="Evaluation Score", fill="Inference Type")

    
    plts[["evaluation-ranking_perc_vs_regulon_set-main-box"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MAIN)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    
    # dPSI thresholds networks
    plts[["evaluation-ranking_perc_vs_regulon_set-dpsi_thresh-box"]] = X %>%
        filter(regulon_set %in% SETS_THRESHOLDS) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator, n_total_regulators, n_total_targets) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_THRESHOLDS)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=0.1, label=n_total_regulators), 
            . %>% distinct(omic_type, eval_direction, eval_type, regulon_set, n_total_regulators), 
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=0.05, label=n_total_targets), 
            . %>% distinct(omic_type, eval_direction, eval_type, regulon_set, n_total_targets), 
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # robustness networks
    plts[["evaluation-ranking_perc_vs_regulon_set-robustness-box"]] = X %>%
        filter(regulon_set %in% SETS_ROBUSTNESS) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator, n_targets_median) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_ROBUSTNESS)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        geom_text(
            aes(y=0.1, label=n_targets_median), 
            . %>% distinct(omic_type, eval_direction, eval_type, regulon_set, n_targets_median), 
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # likelihood networks
    plts[["evaluation-ranking_perc_vs_regulon_set-likelihood-box"]] = X %>%
        filter(regulon_set %in% SETS_LIKELIHOOD) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_LIKELIHOOD)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # mor networks
    plts[["evaluation-ranking_perc_vs_regulon_set-mor-box"]] = X %>%
        filter(regulon_set %in% SETS_MOR) %>%
        filter(n_tails=="two") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MOR)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    # one-tailed networks
    plts[["evaluation-ranking_perc_vs_regulon_set-main_one_tailed-box"]] = X %>%
        filter(regulon_set %in% SETS_MAIN) %>%
        filter(n_tails=="one") %>%
        group_by(omic_type, eval_direction, eval_type, regulon_set, regulator) %>%
        summarize(ranking_perc = median(ranking_perc, na.rm=TRUE)) %>%
        ungroup() %>%
        mutate(regulon_set = factor(regulon_set, levels=SETS_MAIN)) %>%
        ggplot(aes(x=regulon_set, y=ranking_perc, 
                   group=interaction(regulon_set, eval_type))) +
        geom_boxplot(aes(fill=eval_type), width=0.5, outlier.size=0.1, 
                     position=position_dodge(0.5)) +
        fill_palette(PAL_EVAL_TYPE) + 
        theme_pubr() +
        facet_wrap(~omic_type+eval_direction, ncol=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(regulon_set, eval_direction, eval_type, omic_type) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(x="Regulon Set", y="Evaluation Score", fill="Inference Type")
    
    names(plts) = sprintf("%s-%s", omic_type_oi, names(plts))
    
    return(plts)
}


make_plots = function(evaluation){
    plts = list(
        plot_evaluation(evaluation, "EX")
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(evaluation){
    figdata = list(
        "regulon_evaluation" = list(
            "evaluation" = evaluation
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)   
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    omic_types = c("EX")
    for (omic_type_oi in omic_types){
        # main
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set_vs_pert_type-main-box", omic_type_oi), '.pdf', figs_dir, width=6.5, height=12)
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-main-box", omic_type_oi), '.pdf', figs_dir, width=7, height=7)
        # robustness
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-robustness-box", omic_type_oi), '.pdf', figs_dir, width=12, height=7)
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-dpsi_thresh-box", omic_type_oi), '.pdf', figs_dir, width=12, height=7)
        # likelihood
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-likelihood-box", omic_type_oi), '.pdf', figs_dir, width=7, height=7)
        # mor
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-mor-box", omic_type_oi), '.pdf', figs_dir, width=3.5, height=7)
        # one-tailed
        save_plt(plts, sprintf("%s-evaluation-ranking_perc_vs_regulon_set-main_one_tailed-box", omic_type_oi), '.pdf', figs_dir, width=7, height=7)
    }
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        d = file.path(dir,'figdata',x)
        dir.create(d, recursive=TRUE)
        lapply(names(figdata[[x]]), function(nm){
            df = figdata[[x]][[nm]]
            filename = file.path(d, paste0(nm,'.tsv.gz'))
            write_tsv(df, filename)
            
            print(filename)
        })
    })
}


parseargs = function(){
    
    option_list = list( 
        make_option("--evaluation_ex_file", type="character"),
        make_option("--regulators_per_target_robustness_file", type="character"),
        make_option("--targets_per_regulator_robustness_file", type="character"),
        make_option("--regulators_per_target_thresholds_file", type="character"),
        make_option("--targets_per_regulator_thresholds_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    evaluation_ex_file = args[["evaluation_ex_file"]]
    regulators_per_target_robustness_file = args[["regulators_per_target_robustness_file"]]
    targets_per_regulator_robustness_file = args[["targets_per_regulator_robustness_file"]]
    regulators_per_target_thresholds_file = args[["regulators_per_target_thresholds_file"]]
    targets_per_regulator_thresholds_file = args[["targets_per_regulator_thresholds_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    evaluation = read_tsv(evaluation_ex_file)
    splicing_factors = read_tsv(splicing_factors_file)
    
    targets_per_regulator = list(
        read_tsv(targets_per_regulator_robustness_file),
        read_tsv(targets_per_regulator_thresholds_file)
    ) %>% bind_rows()
    regulators_per_target = list(
        read_tsv(regulators_per_target_robustness_file),
        read_tsv(regulators_per_target_thresholds_file)
    ) %>% bind_rows()
    
    targets_per_regulator = targets_per_regulator %>%
        group_by(regulon_set_id) %>%
        summarize(
            n_targets_median = median(n_targets),
            n_total_regulators = n()
        ) %>%
        ungroup()
    n_targets_all = targets_per_regulator %>% 
        filter(regulon_set_id=="experimentally_derived_regulons_pruned-EX") %>% 
        pull(n_targets_median)
    targets_per_regulator = targets_per_regulator %>%
        mutate(perc_targets_median = round(100*n_targets_median/n_targets_all,1))
    
    regulators_per_target = regulators_per_target %>%
        group_by(regulon_set_id) %>%
        summarize(
            n_regulators_median = median(n_regulators),
            n_total_targets = n()
        ) %>%
        ungroup()
    n_regulators_all = regulators_per_target %>% 
        filter(regulon_set_id=="experimentally_derived_regulons_pruned-EX") %>% 
        pull(n_regulators_median)
    regulators_per_target = regulators_per_target %>%
        mutate(perc_regulators_median = round(100*n_regulators_median/n_regulators_all,1))
    
    # prep
    evaluation = evaluation %>%
        mutate(
            regulon_set = gsub("-.*","",regulon_set_id)
        ) %>%
        left_join(splicing_factors, by=c("regulator"="ENSEMBL")) %>%
        mutate(
            sf_class = case_when(
                !is.na(spliceosome_db_complex) ~ "Core",
                in_go_rbp ~ "RBP",
                TRUE ~ "Other"
            )
        )
    
    # plot
    plts = make_plots(evaluation)
    
    # make figdata
    figdata = make_figdata(evaluation)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}