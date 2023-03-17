#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# make figures of evaluate target inference algorithms.
# median accuracy of 0.973 and 0.975 in HepG2 and K562 with threshold of 0.2

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(extrafont)

# variables

# formatting
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_DARK = "darkgreen"
PAL_DUAL = c("darkgreen","#7570B3") # '#1B9E77''#7570B3'

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,"results","sf_targets_inference")

#evaluation_file = file.path(RESULTS_DIR,"files","inference_evaluation","merged.tsv.gz")

# viper_result_file = file.path(RESULTS_DIR,"files","validations","aracne","event_psi","ENCORE","HepG2-regulon_CCLE.tsv.gz")
# regulons_file = file.path(RESULTS_DIR,"files","target_inference","aracne","event_psi","CCLE.tsv.gz")

# encore_logfc_file = file.path(PREP_DIR,'ground_truth_kd','ENCORE',"HepG2",'log2_fold_change_tpm.tsv.gz')

# figs_dir = file.path(RESULTS_DIR,'figures','validations')

##### FUNCTIONS #####
# plot_evaluation = function(evaluation){
#     plts = list()
    
#     X = evaluation
    
#     n_kds = X %>%
#         distinct(sf_target_inference_method, kd_cell_line, KD_ENSEMBL) %>%
#         count(sf_target_inference_method, kd_cell_line) %>%
#         mutate(label=sprintf("%s (n=%s)",kd_cell_line,n))
    
#     plts[["evaluation-thresh_vs_prop_correct-line"]] = X %>%
#         left_join(n_kds, by=c("sf_target_inference_method","kd_cell_line")) %>%
#         ggplot(aes(x=threshold_classification, y=prop_correct, group=KD_ENSEMBL)) +
#         geom_line(size=0.1, color="grey", alpha=0.5) +
#         geom_smooth(aes(color=kd_cell_line, fill=kd_cell_line, group=kd_cell_line), 
#                     se=FALSE, span=0.2, size=LINE_SIZE, linetype="dashed", alpha=0.5, method="loess") +
#         color_palette("Dark2") +
#         theme_pubr(legend="none") +
#         facet_wrap(~label+sf_target_inference_method) +
#         theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="Classification Threshold", y="Proportion Correct") +
#         lims(x=c(0,1), y=c(0,1))
    
#     x = X %>%
#         left_join(n_kds, by=c("kd_cell_line","sf_target_inference_method")) %>% 
#         group_by(sf_target_inference_method, label, kd_cell_line, threshold_classification) %>%
#         summarize(
#             med_tpr=median(tpr, na.rm=TRUE), 
#             med_fpr=median(fpr, na.rm=TRUE),
#             med_recall=median(recall, na.rm=TRUE),
#             med_precision=median(precision, na.rm=TRUE)
#         )
#     plts[["evaluation-fpr_vs_tpr-line"]] = x %>%
#         arrange(threshold_classification) %>%
#         ggplot(aes(x=med_fpr, y=med_tpr)) +
#         geom_line(aes(color=sf_target_inference_method), size=LINE_SIZE, linetype="dashed") +
#         geom_point(aes(color=sf_target_inference_method), size=1) +
#         color_palette("Dark2") +
#         facet_wrap(~label) +
#         theme_pubr() + 
#         theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="FPR", y="TPR", color="Inference Method")
    
#     plts[["evaluation-recall_vs_precision-line"]] = x %>%
#         ggplot(aes(x=med_recall, y=med_precision)) +
#         geom_line(aes(color=sf_target_inference_method), size=LINE_SIZE, linetype="dashed") +
#         geom_point(aes(color=sf_target_inference_method), size=1) +
#         color_palette("Dark2") +
#         facet_wrap(~label) +
#         theme_pubr() + 
#         theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
#         labs(x="Recall", y="Precision", color="Inference Method")
    
#     return(plts)
# }


plot_regulons = function(regulons){
    plts = list()
    
    X = regulons
    
    total_upreg = X[["upstream_regulator"]] %>% unique() %>% length()
    total_targets = X[["target"]] %>% unique() %>% length()
    
    # number of targets per upstream regulator?
    plts[["regulons-n_targets_per_reg-distr"]] = X %>%
        count(upstream_regulator) %>%
        gghistogram(x="n", fill=PAL_DARK, color=NA) +
        labs(x="# Targets per Upstream Regulator", y="Count")
    
    # pleiotropy? how many upstream regulators per target?
    plts[["regulons-pleiotropy-distr"]] = X %>%
        count(target) %>%
        gghistogram(x="n", fill=PAL_DARK, color=NA) +
        labs(x="# Upstream Regulators per Target", y="Count")
    
    return(plts)
}


plot_viper_activities = function(viper_result, encore_logfc){
    plts = list()
    
    X = viper_result %>%
        pivot_longer(
            -c(upstream_regulator), 
            names_to="KD", 
            values_to="viper_activity"
        ) %>%
        group_by(KD) %>%
        mutate(
            # rank viper activity
            viper_activity_ranking = rank(viper_activity),
            # mark if upstream regulator was experimentally validated
            is_validated = upstream_regulator == KD
        ) %>%
        ungroup() %>%
        left_join(
            encore_logfc %>%
            pivot_longer(-ID, names_to="KD", values_to="log2_fc"),
            by = c("upstream_regulator"="ID","KD")
        )
    
    # randomize selection of the same N
    x_validated = X %>% filter(is_validated)
    x_random_non_validated = X %>%
        filter(!is_validated) %>%
        sample_n(size=nrow(x_validated), seed=1234)
    x = rbind(x_validated, x_random_non_validated)
    
    # are the rankings of validated SFs higher than non-validated?
    plts[["viper_activities-activity_vs_is_validated-violin"]] = x %>%
        ggviolin(x="is_validated", y="viper_activity", 
                 fill="is_validated", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black") +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        fill_palette(PAL_DUAL) +
        geom_text(
            aes(label=label, y=-2),
            . %>% count(is_validated) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY,
        ) +
        guides(fill="none") +
        labs(x="Is Validated Experimentally", y="VIPER Activity")
    
    plts[["viper_activities-activity_ranking_vs_is_validated-violin"]] = x %>%
        ggviolin(x="is_validated", y="viper_activity_ranking", 
                 fill="is_validated", color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1, fill=NA, color="black") +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        fill_palette(PAL_DUAL) +
        geom_text(
            aes(label=label, y=-5),
            . %>% count(is_validated) %>% mutate(label=sprintf("n=%s",n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        guides(fill="none") +
        labs(x="Is Validated Experimentally", y="VIPER Activity Ranking")
    
    # does predicted activity correlate with experimental logFC?
    plts[["viper_activities-fc_kds_vs_activity-scatter"]] = X %>%
        ggplot(aes(x=log2_fc, y=viper_activity)) +
        geom_scattermore(
            aes(color=is_validated),
            pixels = c(1000,1000), pointsize=7, alpha=0.5) +
        color_palette(PAL_DUAL) +
        theme_pubr() +
        facet_wrap(~is_validated, scales="free_x") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(color="none") +
        labs(x="Log2FC KD vs WT", y="Viper Activity")
    
    plts[["viper_activities-fc_kds_vs_activity_ranking-scatter"]] = X %>%
        ggplot(aes(x=log2_fc, y=viper_activity_ranking)) +
        geom_scattermore(
            aes(color=is_validated),
            pixels = c(1000,1000), pointsize=7, alpha=0.5) +
        color_palette(PAL_DUAL) +
        theme_pubr() +
        facet_wrap(~is_validated, scales="free_x") +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(color="none") +
        labs(x="Log2FC KD vs WT", y="Viper Activity Ranking")
    
    return(plts)
}


make_plots = function(viper_result, encore_logfc, regulons){
    plts = list(
        plot_viper_activities(viper_result, encore_logfc),
        plot_regulons(regulons)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(viper_result, encore_logfc, regulons){
    figdata = list(
        "target_inference" = list(
            "viper_result" = viper_result
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
#     save_plt(plts, "evaluation-thresh_vs_prop_correct-line", '.pdf', figs_dir, width=8, height=8)
#     save_plt(plts, "evaluation-fpr_vs_tpr-line", '.pdf', figs_dir, width=8, height=6)
#     save_plt(plts, "evaluation-recall_vs_precision-line", '.pdf', figs_dir, width=8, height=6)
    
    save_plt(plts, "viper_activities-activity_vs_is_validated-violin", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "viper_activities-activity_ranking_vs_is_validated-violin", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "viper_activities-fc_kds_vs_activity-scatter", '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, "viper_activities-fc_kds_vs_activity_ranking-scatter", '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, "regulons-n_targets_per_reg-distr", '.pdf', figs_dir, width=4, height=4)
    save_plt(plts, "regulons-pleiotropy-distr", '.pdf', figs_dir, width=4, height=4)
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
        make_option("--viper_result_file", type="character"),
        make_option("--encore_logfc_file", type="character"),
        make_option("--regulons_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    viper_result_file = args[["viper_result_file"]]
    encore_logfc_file = args[["encore_logfc_file"]]
    regulons_file = args[["regulons_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    viper_result = read_tsv(viper_result_file)
    encore_logfc = read_tsv(encore_logfc_file)
    regulons = read_tsv(regulons_file)
    
    # plot
    plts = make_plots(viper_result, encore_logfc, regulons)
    
    # make figdata
    figdata = make_figdata(viper_result, encore_logfc, regulons)

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}