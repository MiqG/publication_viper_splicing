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
    'aracne_regulons_PANCAN_STN',
    'mlr_regulons_PANCAN_STN',
    'experimentally_derived_regulons_pruned'
)

SETS_CLIP = c(
    "experimentally_derived_regulons_pruned",
    "postar3_clip_regulons"
    # postar3_and_experimental_regulons,
    # experimental_without_postar3_regulons
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

METHODS_ACTIVITY = c(
    "gsea",
    "correlation_spearman",
    "correlation_pearson",
    "viper"
)

SF_CLASS = c("Core", "RBP", "Other")

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "#616161"
PAL_EVAL_TYPE = c(
    "random" = "lightgrey",
    "real" = "orange"
)

PAL_METHODS_ACTIVITY = c(
    "gsea"="#383F51",
    "correlation_spearman"="#B9BAA3",
    "correlation_pearson"="#8E8DBE",
    "viper"="#A22C29"
)

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
    
    # evaluation by dataset
    x = X %>%
        group_by(regulon_set, method_activity, curves_type, eval_direction, eval_type, signature_id) %>%
        summarize(auc_roc = median(auc_roc, na.rm=TRUE)) %>%    
        ungroup()
    
    plts[["evaluation-general-median_auc_roc-box"]] = x %>%
        filter(regulon_set%in%SETS_MAIN & eval_type=="real") %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_MAIN),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        ggplot(aes(x=regulon_set, y=auc_roc, group=interaction(regulon_set, method_activity))) +
        geom_boxplot(aes(color=method_activity), fill=NA, outlier.shape=NA, position=position_dodge(0.9)) +
        geom_point(aes(color=method_activity, shape=signature_id), size=0.5, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        color_palette(PAL_METHODS_ACTIVITY) + 
        geom_text(
            aes(y = 0.3, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="median(ROC AUC)", color="Inference Type", shape="Held-out Dataset")
    
    # evaluation by SF class of empirical splicing networks
    plts[["evaluation-held_out_ds-raw_auc_roc-box"]] = X %>%
        filter(eval_type=="real" & regulon_set=="experimentally_derived_regulons_pruned") %>%
        mutate(
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        ggplot(aes(x=signature_id, y=auc_roc, group=interaction(signature_id, method_activity))) +
        geom_boxplot(aes(color=method_activity), fill=NA, outlier.shape=NA, position=position_dodge(0.9)) +
        geom_point(aes(color=method_activity, shape=signature_id), size=0.5, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        color_palette(PAL_METHODS_ACTIVITY) + 
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(signature_id, method_activity, regulon_set, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Held-out Dataset", y="ROC AUC", color="Inference Type", shape="Held-out Dataset")

    # evaluation by SF class of empirical splicing networks
    plts[["evaluation-sf_class-raw_auc_roc-box"]] = X %>%
        filter(eval_type=="real" & regulon_set=="experimentally_derived_regulons_pruned" & eval_direction=="regulators") %>%
        mutate(
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY),
            sf_class = factor(sf_class, levels=SF_CLASS)
        ) %>%
        ggplot(aes(x=sf_class, y=auc_roc, group=interaction(sf_class, method_activity))) +
        geom_boxplot(aes(color=method_activity), fill=NA, outlier.shape=NA, position=position_dodge(0.9)) +
        geom_point(aes(color=method_activity, shape=signature_id), size=0.5, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        color_palette(PAL_METHODS_ACTIVITY) + 
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(sf_class, method_activity, regulon_set, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="SF Class", y="ROC AUC", color="Inference Type", shape="Held-out Dataset")


    # evaluation by SF class and held-out datasets of empirical splicing networks
    plts[["evaluation-sf_class_vs_held_out_ds-raw_auc_roc-box"]] = X %>%
        filter(eval_type=="real" & regulon_set=="experimentally_derived_regulons_pruned" & eval_direction=="regulators") %>%
        mutate(
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY),
            sf_class = factor(sf_class, levels=SF_CLASS)
        ) %>%
        ggplot(aes(x=sf_class, y=auc_roc, group=interaction(sf_class, method_activity))) +
        geom_boxplot(aes(color=method_activity), fill=NA, outlier.shape=NA, position=position_dodge(0.9)) +
        geom_point(aes(color=method_activity, shape=signature_id), size=0.5, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        color_palette(PAL_METHODS_ACTIVITY) + 
        geom_text(
            aes(y = -0.1, label=label), 
            . %>% 
            count(sf_class, signature_id, method_activity, regulon_set, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction+signature_id) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="SF Class", y="ROC AUC", color="Inference Type", shape="Held-out Dataset")
    
    
    # (TODO) CLIP benchmark
    plts[["evaluation-clip-median_auc_roc-box"]] = x %>%
        filter(eval_type=="real" & regulon_set%in%SETS_CLIP) %>%
        mutate(
            regulon_set = factor(regulon_set, levels=SETS_CLIP),
            method_activity = factor(method_activity, levels=METHODS_ACTIVITY)
        ) %>%
        ggplot(aes(x=regulon_set, y=auc_roc, group=interaction(regulon_set, method_activity))) +
        geom_boxplot(aes(color=method_activity), fill=NA, outlier.shape=NA, 
                     position=position_dodge2(0.9, preserve="single")) +
        geom_point(aes(color=method_activity, shape=signature_id), size=0.5, 
                   position=position_jitterdodge(dodge.width=0.9, jitter.width=0.1)) +
        color_palette(PAL_METHODS_ACTIVITY) + 
        geom_text(
            aes(y = 0.3, label=label), 
            . %>% 
            count(method_activity, regulon_set, eval_type, eval_direction) %>% 
            mutate(label=paste0("n=",n)),
            position=position_dodge(0.9), size=FONT_SIZE, family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle = 45) +
        facet_wrap(~eval_direction) +  
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Method", y="median(ROC AUC)", color="Inference Type", shape="Held-out Dataset")
    
    return(plts)
}


make_plots = function(evaluation){
    plts = list(
        plot_evaluation(evaluation)
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
        ) %>% 
        # drop summarized evaluation
        filter(curves_type=="by_group" & n_tails=="two")
    
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